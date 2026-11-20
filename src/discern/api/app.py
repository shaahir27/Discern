"""Serve the fixed research archive and the local Jinja interface."""

from collections import OrderedDict
from collections.abc import AsyncIterator, Awaitable, Callable
from contextlib import asynccontextmanager
from dataclasses import dataclass
from pathlib import Path
from secrets import token_urlsafe
from urllib.parse import parse_qs, urlencode

import numpy as np
from fastapi import FastAPI, HTTPException, Query, Request
from fastapi.concurrency import run_in_threadpool
from fastapi.exceptions import RequestValidationError
from fastapi.responses import HTMLResponse, JSONResponse, RedirectResponse, Response
from fastapi.staticfiles import StaticFiles
from fastapi.templating import Jinja2Templates
from starlette.exceptions import HTTPException as StarletteHTTPException

from ..decision import Predictor
from ..documentation import model_card
from ..importers import detect_and_parse_upload
from ..label_shift import deduplicate_patient_probabilities, estimate_prevalence_mlls
from ..longitudinal import analyze_patient_timelines
from ..panels import PANELS, audit_sample_panel
from ..prevalence import adjust, archived_predictive_values
from ..quality import check_csv_quality
from ..research import (
    BANNER,
    MAX_CSV_BYTES,
    MODE,
    PROFILE_COLUMNS,
    REPORTS,
    Archive,
    UploadError,
    csv_text,
    read_batch,
    read_json,
    workspace,
)
from ..surveillance import batch_recurrence_audit

UI = Path(__file__).resolve().parents[1] / "ui"


def asset_url(name: str) -> str:
    version = (UI / "static" / name).stat().st_mtime_ns
    return f"/static/{name}?v={version}"


templates = Jinja2Templates(directory=UI / "templates")
templates.env.globals.update(banner=BANNER, mode=MODE, max_csv_bytes=MAX_CSV_BYTES, asset_url=asset_url)
UPLOAD_COOKIE = "discern_upload"
DEFAULT_PREVALENCE = 0.05


@dataclass
class UploadedBatch:
    id: str
    owner: str
    quality: dict
    records: dict[str, dict]
    timelines: dict[str, dict]
    prevalence: dict
    rows: list[dict]
    alerts: list[dict]
    illustrative: bool


@asynccontextmanager
async def lifespan(api: FastAPI) -> AsyncIterator[None]:
    root = workspace()
    api.state.archive = Archive(root)
    api.state.predictor = Predictor(root / "models/lr_model", root / "results", config_path=root / "config.yaml")
    api.state.uploads = OrderedDict()
    yield
    api.state.uploads.clear()


app = FastAPI(title="Discern Research Explorer", lifespan=lifespan, docs_url=None, redoc_url=None, openapi_url=None)
app.mount("/static", StaticFiles(directory=UI / "static"), name="static")


@app.middleware("http")
async def research_header(request: Request, call_next: Callable[[Request], Awaitable[Response]]) -> Response:
    response = await call_next(request)
    response.headers["X-Research-Use"] = BANNER
    response.headers["X-Content-Type-Options"] = "nosniff"
    response.headers["Content-Security-Policy"] = "default-src 'self'; object-src 'none'; frame-ancestors 'none'; form-action 'self'"
    if request.method == "POST" or "batch" in request.query_params or request.url.path.startswith(("/result", "/patient/", "/api/patient/")):
        response.headers["Cache-Control"] = "no-store"
    return response


@app.exception_handler(StarletteHTTPException)
async def http_error(_request: Request, error: StarletteHTTPException) -> JSONResponse:
    return JSONResponse(
        {"detail": error.detail, "banner": BANNER}, status_code=error.status_code, headers=error.headers,
    )


@app.exception_handler(UploadError)
async def upload_error(_request: Request, error: UploadError) -> JSONResponse:
    return JSONResponse({
        "detail": str(error), "banner": BANNER,
        "error": {"code": error.code, "line": error.line, "column": error.column},
    }, status_code=error.status)


@app.exception_handler(RequestValidationError)
async def invalid_request(_request: Request, _error: RequestValidationError) -> JSONResponse:
    return JSONResponse({"detail": "Check the request format.", "banner": BANNER}, status_code=422)


@app.exception_handler(Exception)
async def server_error(_request: Request, _error: Exception) -> JSONResponse:
    return JSONResponse({"detail": "The archive is temporarily unavailable.", "banner": BANNER}, status_code=500,
                        headers={"X-Research-Use": BANNER})


def archive(request: Request) -> Archive:
    return request.app.state.archive


def render(request: Request, name: str, *, status: int = 200, **context: object) -> Response:
    return templates.TemplateResponse(request=request, name=f"{name}.html", context=context, status_code=status)


def csv_response(text: str, filename: str) -> Response:
    return Response(text, media_type="text/csv", headers={"Content-Disposition": f'attachment; filename="{filename}"'})


def queue_csv(rows: list[dict]) -> Response:
    return csv_response(csv_text(["rank", "sample_id", "call", "call_label", "probability", "fold"], rows), "replay-queue.csv")


async def small_body(request: Request, limit: int = MAX_CSV_BYTES) -> bytes:
    payload = bytearray()
    async for chunk in request.stream():
        if len(payload) + len(chunk) > limit:
            raise UploadError("file_too_large", "Choose a smaller CSV of recorded sample IDs.", status=413)
        payload.extend(chunk)
    return bytes(payload)


@app.get("/health")
def health(request: Request) -> dict:
    return {"status": "ok", "versions": archive(request).versions, "mode": MODE, "banner": BANNER}


@app.get("/replay/{sample_id}")
def replay(request: Request, sample_id: str) -> dict:
    try:
        return archive(request).replay(sample_id)
    except KeyError:
        raise HTTPException(404, "That sample is outside the replay collection.") from None


@app.get("/replay-samples.csv")
def replay_samples(request: Request) -> Response:
    return csv_response(csv_text(["sample_id"], list(archive(request).roster.values())), "replay-samples.csv")


@app.get("/example.csv")
def example(request: Request) -> Response:
    data = archive(request)
    return csv_response(csv_text(PROFILE_COLUMNS, data.profile(data.example_id)), "discern-example.csv")


@app.get("/evidence/{gene}/{mutation}")
def evidence(request: Request, gene: str, mutation: str) -> dict:
    if len(gene) > 12 or len(mutation) > 512:
        raise HTTPException(422, "Use a gene and variant from an archived evidence card.")
    return {**archive(request).card(f"{gene}@{mutation}"), "banner": BANNER}


@app.post("/batch")
async def batch(request: Request) -> Response:
    if request.headers.get("content-type", "").split(";", 1)[0].strip().lower() != "text/csv":
        raise UploadError("unsupported_media_type", "Send a text/csv body containing recorded sample IDs.", status=415)
    data = archive(request)
    return queue_csv(data.queue(data.imported_ids(await small_body(request))))


@app.get("/", response_class=HTMLResponse)
def home(request: Request) -> Response:
    return sample_home(request)


def sample_home(request: Request, error: UploadError | None = None) -> Response:
    data = archive(request)
    return render(request, "home", title="Check a sample", active="home", roster=data.roster.values(),
                  example_id=data.example_id, error=str(error) if error else "", status=error.status if error else 200)


async def form_csv(request: Request) -> tuple[bytes, str]:
    if request.headers.get("content-type", "").split(";", 1)[0].strip().lower() != "application/x-www-form-urlencoded":
        raise UploadError("unsupported_media_type", "Use the CSV upload form.", status=415)
    # Each CSV byte can expand to three bytes when the browser percent-encodes the form.
    body = await small_body(request, limit=3 * MAX_CSV_BYTES + len("csv=") + len("&assay=icmr_nirt_minion") + 64)
    try:
        form = parse_qs(body.decode("utf-8"), max_num_fields=2, strict_parsing=True,
                        keep_blank_values=True, errors="strict")
    except ValueError:
        raise UploadError("invalid_form", "Open the CSV upload form and submit the file again.") from None
    if set(form) not in ({"csv"}, {"csv", "assay"}) or any(len(values) != 1 for values in form.values()):
        raise UploadError("invalid_form", "Open the CSV upload form and submit the file again.")
    raw_csv = form["csv"][0].encode("utf-8")
    assay = form.get("assay", ["wgs"])[0].strip() or "wgs"
    return raw_csv, assay


@app.post("/review", response_class=HTMLResponse)
async def review(request: Request) -> Response:
    return await review_upload(request, single_redirect=True)


@app.post("/queue", response_class=HTMLResponse)
async def queue_import(request: Request) -> Response:
    return await review_upload(request, single_redirect=False)


@app.post("/longitudinal-example")
async def open_longitudinal_example(request: Request) -> Response:
    return await review_upload(request, single_redirect=False, illustrative=True)


async def review_upload(request: Request, *, single_redirect: bool, illustrative: bool = False) -> Response:
    data = archive(request)
    try:
        if illustrative:
            payload, assay = longitudinal_example_path().read_bytes(), "wgs"
        else:
            payload, assay = await form_csv(request)
        parsed = await run_in_threadpool(parse_upload, payload, assay, data)
        rows, quality, archived = parsed
        if archived:
            if single_redirect and len(rows) == 1:
                destination = "/result?" + urlencode({"sample_id": next(iter(rows))})
                return RedirectResponse(destination, status_code=303)
            return render(request, "queue", title="Review queue", active="queue",
                          rows=data.queue(list(rows)), scope="import", upload_quality=quality)
        predictor = request.app.state.predictor
        records, queue_rows, timelines, prevalence, alerts = await run_in_threadpool(
            score_upload, rows, predictor,
        )
    except UploadError as error:
        if single_redirect:
            return sample_home(request, error)
        return render(request, "queue", status=error.status, title="Review queue", active="queue", scope="all",
                      rows=data.queue(list(data.roster)), error=str(error))

    uploads = request.app.state.uploads
    batch_id = token_urlsafe(18)
    previous = uploaded_batch(request)
    owner = previous.owner if previous else token_urlsafe(24)
    uploads[batch_id] = UploadedBatch(
        batch_id, owner, quality, records, timelines, prevalence, queue_rows, alerts, illustrative,
    )
    # shortcut: per-process cache for the local, single-worker explorer, use shared storage for multi-worker hosting.
    if len(uploads) > 16:
        uploads.popitem(last=False)
    if illustrative and timelines:
        patient_id = min(timelines)
        destination = request.url_for("patient_timeline_view", patient_id=patient_id)
        response = RedirectResponse(str(destination.include_query_params(batch=batch_id)), status_code=303)
    elif single_redirect and len(records) == 1:
        destination = "/result?" + urlencode({"sample_id": next(iter(records)), "batch": batch_id})
        response = RedirectResponse(destination, status_code=303)
    else:
        response = render_uploaded_queue(request, uploads[batch_id], batch_id)
    response.set_cookie(UPLOAD_COOKIE, owner, httponly=True, samesite="strict",
                        secure=request.url.scheme == "https")
    return response


def parse_upload(payload: bytes, assay: str, data: Archive) -> tuple[dict[str, dict], dict, bool]:
    parsed = detect_and_parse_upload(payload)
    if parsed["format"] in ("discern_csv", "empty"):
        headers, rows = read_batch(payload)
        quality = check_csv_quality(payload)
        if all(row["sample_id"] in data.roster for _, row in rows):
            return {sid: {} for sid in data.imported_ids(payload)}, quality, True
        for issue in quality["issues"]:
            if issue["type"] != "duplicate":
                raise UploadError("invalid_profile", "Check this field in the uploaded profile.",
                                  line=issue["line"], column=issue["field"])
        inferred_assay = assay
        sample_ids: list[str] = []
    else:
        sample_ids = parsed["sample_ids"]
        rows = list(enumerate(parsed["rows"], 2))
        headers = PROFILE_COLUMNS
        inferred_assay = parsed["inferred_assay"] if assay == "wgs" else assay
        quality = check_csv_quality(
            csv_text(PROFILE_COLUMNS, parsed["rows"]) if rows
            else csv_text(["sample_id"], [{"sample_id": sid} for sid in sample_ids]),
        )
        for issue in quality["issues"]:
            if issue["type"] != "duplicate":
                raise UploadError("invalid_profile", "Check this field in the uploaded profile.",
                                  line=issue["line"], column=issue["field"])
        quality["format"] = parsed["format"]
    samples: dict[str, dict] = {}
    for sample_id in sample_ids:
        samples[sample_id] = {"variants": [], "assay": inferred_assay}
    for line, row in rows:
        sample_id = row["sample_id"]
        sample = samples.setdefault(sample_id, {"variants": [], "assay": inferred_assay})
        for field in ("patient_id", "collection_date", "assay"):
            value = row.get(field, "")
            if not value:
                continue
            if field == "assay" and value not in PANELS:
                raise UploadError("invalid_assay", "Choose a supported assay.",
                                  line=line, column=field)
            if field == "patient_id" and (len(value) > 64 or "/" in value or value in {".", ".."}):
                raise UploadError("invalid_patient_id", "Use a patient identifier without path separators.",
                                  line=line, column=field)
            if field in sample and (field != "assay" or sample.get("assay_in_file")) and sample[field] != value:
                raise UploadError("inconsistent_metadata", "Use the same metadata on all rows for this sample.",
                                  line=line, column=field)
            sample[field] = value
            if field == "assay":
                sample["assay_in_file"] = True
        if headers != ["sample_id"]:
            sample["variants"].append({
                "gene": row["gene"], "mutation": row["mutation"], "frs": float(row["frs"]),
            })
    return samples, quality, False


def score_upload(samples: dict[str, dict], predictor: Predictor) -> tuple[dict, list, dict, dict, list]:
    records = {}
    queue_rows = []
    visits = []
    for sid, sample in samples.items():
        variants = sample["variants"]
        assay = sample["assay"]
        pred = predictor.predict(variants, prevalence=DEFAULT_PREVALENCE, sample_id=sid, assay=assay)
        model = pred["model"]
        fractions: dict[str, float] = {}
        for variant in variants:
            key = variant["gene"] + "@" + variant["mutation"].split("@")[-1]
            fractions[key] = max(fractions.get(key, 0.0), variant["frs"])
        basic_record = {
            "sample_id": sid,
            "call": pred["call"],
            "call_label": pred["call_label"],
            "probability": model["probability_at_prevalence"],
            "fold": "Live",
            "source": "Uploaded profile",
            "phenotype": "Live prediction",
            "is_live": True,
            "patient_id": sample.get("patient_id"),
        }
        queue_rows.append(basic_record)
        records[sid] = {
            **basic_record,
            "probability_train": model["probability_train_prevalence"],
            "score_raw": model["score_raw"],
            "evidence": [{**card, "frs": fractions[card["variant"]]} for card in pred["evidence"]],
            "tier_a": pred["tier_a"],
            "clinical_reflex": pred["clinical_reflex"],
            "panel_audit": pred["panel_audit"],
            "mode": "Live model evaluation (Discern Tier A + Tier B)",
            "banner": BANNER,
            "score_context": {
                "reference_prevalence": predictor.pi_tr,
                "uncertainty": None,
                "uncertainty_status": "live",
            },
            "warnings": pred["warnings"],
        }
        visits.append({
            **basic_record,
            "collection_date": sample.get("collection_date"),
            "probability_train": model["probability_train_prevalence"],
            "variants": variants,
        })

    timelines = analyze_patient_timelines(visits)
    for row in queue_rows:
        if timeline := timelines.get(row["patient_id"]):
            row["n_visits"] = timeline["n_visits"]
            row["has_longitudinal_alerts"] = timeline["has_alerts"]
            records[row["sample_id"]]["longitudinal_timeline"] = timeline
    alerts = batch_recurrence_audit(
        {sid: row["variants"] for sid, row in samples.items()}, predictor.ev,
        patient_map={sid: row.get("patient_id") or sid for sid, row in samples.items()},
    )

    probabilities = deduplicate_patient_probabilities(visits)
    prevalence = estimate_prevalence_mlls(probabilities, predictor.pi_tr)
    if prevalence["status"] == "estimated":
        metrics = predictor.cat_metrics
        prevalence["ppv_hat"], prevalence["npv_hat"] = archived_predictive_values(metrics, prevalence["pi_hat"])
        prevalence["ppv_low"], prevalence["npv_high"] = archived_predictive_values(metrics, prevalence["ci95_low"])
        prevalence["ppv_high"], prevalence["npv_low"] = archived_predictive_values(metrics, prevalence["ci95_high"])
        prevalence["curve_points"] = [
            {"pi": float(pi), "ppv": archived_predictive_values(metrics, float(pi))[0]}
            for pi in np.linspace(0.005, 0.25, 40)
        ]
    queue_rows.sort(key=lambda row: (row["call"] != "U", -row["probability"], row["sample_id"]))
    for rank, row in enumerate(queue_rows, 1):
        row["rank"] = rank if row["call"] == "U" else None
    return records, queue_rows, timelines, prevalence, alerts


def uploaded_batch(request: Request, batch_id: str | None = None) -> UploadedBatch | None:
    uploads = request.app.state.uploads
    owner = request.cookies.get(UPLOAD_COOKIE)
    if not owner:
        return None
    if batch_id is not None:
        batch = uploads.get(batch_id)
        return batch if batch and batch.owner == owner else None
    return next((batch for batch in reversed(uploads.copy().values()) if batch.owner == owner), None)


@app.get("/result", response_class=HTMLResponse)
def result(
    request: Request, sample_id: str, batch: str | None = None,
    prevalence: float | None = Query(None, ge=0.001, le=0.5),
) -> Response:
    data = archive(request)
    manual_prevalence = prevalence
    prevalence = DEFAULT_PREVALENCE if prevalence is None else prevalence
    if batch is not None:
        upload = uploaded_batch(request, batch)
        if upload is None or sample_id not in upload.records:
            return render(request, "error", status=404, title="Upload unavailable", active="home",
                          detail="Open the result from your upload or submit the file again.")
        predictor = request.app.state.predictor
        record = upload.records[sample_id].copy()
        estimated = upload.prevalence.get("pi_hat")
        prefilled = manual_prevalence is None and estimated is not None and upload.prevalence["converged"]
        if prefilled:
            prevalence = min(0.5, max(0.001, estimated))
        record["probability"] = float(adjust(record["probability_train"], predictor.pi_tr, prevalence))
        ppv, npv = archived_predictive_values(predictor.cat_metrics, prevalence)
        return render(request, "result", active="home", title="Bedaquiline result", record=record,
                      prevalence=prevalence, ppv=ppv, npv=npv, upload_quality=upload.quality,
                      prefilled_from_batch=prefilled, batch_id=batch)

    try:
        record = data.replay(sample_id)
    except KeyError:
        return render(request, "error", status=404, title="Recorded isolate unavailable", active="home",
                      detail="Choose an isolate from the replay collection.")
    record["panel_audit"] = audit_sample_panel(data.profile(sample_id), "wgs")
    ppv, npv = archived_predictive_values(data.reports.observed("catalogue"), prevalence)
    return render(request, "result", active="home", title="Bedaquiline result", record=record,
                  prevalence=prevalence, ppv=ppv, npv=npv)


@app.get("/queue", response_class=HTMLResponse)
def queue(request: Request, scope: str = "all", download: bool = False, batch: str | None = None) -> Response:
    if batch is not None:
        upload = uploaded_batch(request, batch)
        if upload is None:
            return render(request, "error", status=404, title="Upload unavailable", active="queue",
                          detail="Open the queue from your upload or submit the file again.")
        return render_uploaded_queue(request, upload, batch)
    data = archive(request)
    ids = [sid for sid, row in data.roster.items() if scope != "unresolved" or row["call"] == "U"]
    rows = data.queue(ids)
    if download:
        return queue_csv(rows)
    return render(request, "queue", title="Review queue", active="queue", rows=rows, scope=scope)


def render_uploaded_queue(request: Request, upload: UploadedBatch, batch_id: str) -> Response:
    return render(request, "queue", title="Review queue", active="queue", rows=upload.rows, scope="import",
                  upload_quality=upload.quality, batch_prevalence=upload.prevalence,
                  recurrence_alerts=upload.alerts, batch_id=batch_id, illustrative=upload.illustrative)


@app.get("/patient/{patient_id}", response_class=HTMLResponse)
def patient_timeline_view(request: Request, patient_id: str, batch: str | None = None) -> Response:
    upload = uploaded_batch(request, batch)
    if upload is None or patient_id not in upload.timelines:
        return render(request, "error", status=404, title="Patient timeline unavailable", active="queue",
                      detail="Choose a patient with multiple visits from your uploaded queue.")
    return render(request, "patient_timeline", title=f"Patient timeline for {patient_id}", active="queue",
                  timeline=upload.timelines[patient_id], patient_id=patient_id, batch_id=upload.id,
                  illustrative=upload.illustrative)


@app.get("/api/patient/{patient_id}")
def api_patient_timeline(request: Request, patient_id: str, batch: str | None = None) -> dict:
    upload = uploaded_batch(request, batch)
    if upload is None or patient_id not in upload.timelines:
        raise HTTPException(404, "Patient timeline unavailable.")
    return {**upload.timelines[patient_id], "banner": BANNER}


def longitudinal_example_path() -> Path:
    demo_csv = workspace() / "demo/longitudinal_example.csv"
    if not demo_csv.is_file():
        raise HTTPException(404, "Longitudinal example CSV not found.")
    return demo_csv


@app.get("/longitudinal-example.csv")
def longitudinal_example_csv() -> Response:
    return csv_response(longitudinal_example_path().read_text(encoding="utf-8"),
                        "discern-longitudinal-example.csv")


@app.get("/priority", response_class=HTMLResponse)
def priority(request: Request, page: int = Query(1, ge=1, le=1000)) -> Response:
    reports = archive(request).reports
    columns, rows, more = reports.table(REPORTS["priority"], "top_sets.csv", page)
    return render(request, "priority", title="Variant evidence priorities", active="priority",
                  columns=columns, rows=rows, page=page, more=more, overview=reports.observed("priority"))


@app.get("/performance", response_class=HTMLResponse)
def performance(
    request: Request, report: str = REPORTS["catalogue"], table: str = "",
    page: int = Query(1, ge=1, le=1000),
) -> Response:
    reports = archive(request).reports
    if report not in reports.paths:
        return render(request, "error", status=404, title="Report unavailable", active="performance",
                      detail="Choose a registered report from Performance.")
    record = reports.get(report)
    names = reports.tables(report, record)
    table = table if table in names else next(iter(names), "")
    columns, rows, more = reports.table(report, table, page) if table else ([], [], False)
    return render(request, "performance", title="Performance archive", active="performance",
                  approaches=reports.approaches, report=report, record=record, table=table,
                  names=names, columns=columns, rows=rows, more=more, page=page)


@app.get("/model-card", response_class=HTMLResponse)
def card(request: Request) -> Response:
    return render(request, "model_card", title="Model card", active="model-card",
                  sections=model_card(archive(request).reports))


@app.get("/panels", response_class=HTMLResponse)
def panels_page(request: Request) -> Response:
    panels_json_path = workspace() / "results/approaches/panel_audit/e24_panel_coverage/results.json"
    if not panels_json_path.is_file():
        return render(request, "error", status=404, title="Panel coverage unavailable", active="panels",
                      detail="The saved panel coverage report is missing.")
    panels_data = read_json(panels_json_path)
    return render(request, "panels", title="Panel coverage", active="panels", panels_data=panels_data)


@app.get("/about", response_class=HTMLResponse)
def about(request: Request) -> Response:
    return render(request, "about", title="About Discern", active="about")
