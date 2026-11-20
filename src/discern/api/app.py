"""Serve the fixed research archive and the local Jinja interface."""

from collections.abc import AsyncIterator, Awaitable, Callable
from contextlib import asynccontextmanager
from decimal import Decimal, InvalidOperation
from pathlib import Path
from urllib.parse import parse_qs, urlencode

from fastapi import FastAPI, HTTPException, Query, Request
from fastapi.exceptions import RequestValidationError
from fastapi.responses import HTMLResponse, JSONResponse, RedirectResponse, Response
from fastapi.staticfiles import StaticFiles
from fastapi.templating import Jinja2Templates
from starlette.exceptions import HTTPException as StarletteHTTPException

from ..decision import Predictor
from ..documentation import model_card
from ..importers import detect_and_parse_upload
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
    workspace,
)

UI = Path(__file__).resolve().parents[1] / "ui"
templates = Jinja2Templates(directory=UI / "templates")
templates.env.globals.update(banner=BANNER, mode=MODE, max_csv_bytes=MAX_CSV_BYTES)


@asynccontextmanager
async def lifespan(api: FastAPI) -> AsyncIterator[None]:
    api.state.archive = Archive(workspace())
    api.state.predictor = Predictor()
    api.state.live_results = {}
    api.state.live_quality = {}
    yield


app = FastAPI(title="Discern Research Explorer", lifespan=lifespan, docs_url=None, redoc_url=None, openapi_url=None)
app.mount("/static", StaticFiles(directory=UI / "static"), name="static")


@app.middleware("http")
async def research_header(request: Request, call_next: Callable[[Request], Awaitable[Response]]) -> Response:
    response = await call_next(request)
    response.headers["X-Research-Use"] = BANNER
    response.headers["X-Content-Type-Options"] = "nosniff"
    response.headers["Content-Security-Policy"] = "default-src 'self'; object-src 'none'; frame-ancestors 'none'; form-action 'self'"
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
    if not (set(form) == {"csv"} or set(form) == {"csv", "assay"}):
        raise UploadError("invalid_form", "Open the CSV upload form and submit the file again.")
    raw_csv = form["csv"][0].encode("utf-8")
    assay = form.get("assay", ["wgs"])[0].strip() or "wgs"
    return raw_csv, assay


@app.post("/review", response_class=HTMLResponse)
async def review(request: Request) -> Response:
    data = archive(request)
    pred_service = request.app.state.predictor
    try:
        raw_csv, assay = await form_csv(request)
    except UploadError as error:
        return sample_home(request, error)

    parsed = detect_and_parse_upload(raw_csv)
    sample_meta: dict[str, dict] = {}
    if parsed["format"] not in ("discern_csv", "empty"):
        if assay == "wgs" and parsed["inferred_assay"] != "wgs":
            assay = parsed["inferred_assay"]
        rows_by_sample: dict[str, list[dict]] = {}
        for r in parsed["rows"]:
            sid = r["sample_id"]
            sample_meta.setdefault(sid, {})
            if "patient_id" in r and r["patient_id"]:
                sample_meta[sid]["patient_id"] = r["patient_id"]
            if "collection_date" in r and r["collection_date"]:
                sample_meta[sid]["collection_date"] = r["collection_date"]
            rows_by_sample.setdefault(sid, []).append({
                "gene": r["gene"],
                "mutation": r["mutation"],
                "frs": float(r["frs"]),
            })
        for sid in parsed.get("sample_ids", []):
            rows_by_sample.setdefault(sid, [])
        quality = {
            "completeness_score": 100.0,
            "total_rows": len(parsed["rows"]),
            "valid_rows": len(parsed["rows"]),
            "missing_count": 0,
            "duplicate_count": 0,
            "formatting_count": 0,
            "total_issues": 0,
            "issues": [],
            "is_valid": True,
            "summary": f"Imported {len(parsed['rows'])} variants from {parsed['format'].replace('_', ' ').title()}.",
            "format": parsed["format"],
            "headers": ["sample_id", "gene", "mutation", "frs"],
        }
    else:
        try:
            quality = check_csv_quality(raw_csv)
            headers, rows = read_batch(raw_csv)
        except UploadError as error:
            return sample_home(request, error)

        first_lines: dict[str, int] = {}
        for line, row in rows:
            first_lines.setdefault(row["sample_id"], line)

        all_known = all(sid in data.roster for sid in first_lines)
        if all_known:
            try:
                ids = data.imported_ids(raw_csv)
            except UploadError as error:
                return sample_home(request, error)
            for sid in ids:
                request.app.state.live_quality[sid] = quality
            if len(ids) == 1:
                destination = "/result?" + urlencode({"sample_id": ids[0]})
                return RedirectResponse(destination, status_code=303)
            return render(request, "queue", title="Review queue", active="queue", rows=data.queue(ids), scope="import", upload_quality=quality)

        # Live scoring for new isolate uploads
        rows_by_sample = {}
        if headers == PROFILE_COLUMNS or (len(headers) >= 4 and headers[:4] == PROFILE_COLUMNS):
            for line, row in rows:
                sid = row["sample_id"]
                sample_meta.setdefault(sid, {})
                if "patient_id" in row and row["patient_id"]:
                    sample_meta[sid]["patient_id"] = row["patient_id"]
                if "collection_date" in row and row["collection_date"]:
                    sample_meta[sid]["collection_date"] = row["collection_date"]
                try:
                    fraction = Decimal(row["frs"])
                except InvalidOperation:
                    fraction = Decimal("NaN")
                if not fraction.is_finite() or not 0 <= fraction <= 1:
                    return sample_home(
                        request,
                        UploadError("invalid_fraction", "Enter a number between zero and one.", line=line, column="frs"),
                    )
                rows_by_sample.setdefault(sid, []).append({
                    "gene": row["gene"],
                    "mutation": row["mutation"],
                    "frs": float(fraction),
                })
        else:
            for line, row in rows:
                sid = row["sample_id"]
                rows_by_sample.setdefault(sid, [])

    new_rows = []
    sample_records_for_longitudinal = []
    for sid, variants in rows_by_sample.items():
        meta = sample_meta.get(sid, {})
        pred = pred_service.predict(variants, prevalence=0.05, sample_id=sid, assay=assay)
        evidence = []
        for card, var in zip(pred["evidence"], variants):
            evidence.append({**card, "frs": float(var.get("frs", 1.0))})
        if not evidence and pred["evidence"]:
            evidence = pred["evidence"]

        record = {
            "sample_id": sid,
            "call": pred["call"],
            "call_label": pred["call_label"],
            "probability": pred["model"]["probability_at_prevalence"],
            "probability_train": pred["model"]["probability_train_prevalence"],
            "score_raw": pred["model"]["score_raw"],
            "fold": "Live",
            "source": "Uploaded CSV",
            "phenotype": "Live prediction",
            "evidence": evidence,
            "tier_a": pred["tier_a"],
            "clinical_reflex": pred.get("clinical_reflex", []),
            "panel_audit": pred.get("panel_audit"),
            "assay": assay,
            "mode": "Live model evaluation (Discern Tier A + Tier B)",
            "banner": BANNER,
            "score_context": {
                "reference_prevalence": pred_service.pi_tr,
                "uncertainty": None,
                "uncertainty_status": "live",
            },
            "warnings": pred["warnings"],
            "is_live": True,
            "raw_variants": variants,
        }
        request.app.state.live_results[sid] = record
        request.app.state.live_quality[sid] = quality
        new_rows.append({
            "sample_id": sid,
            "call": pred["call"],
            "call_label": pred["call_label"],
            "probability": pred["model"]["probability_at_prevalence"],
            "fold": "Live",
            "source": "Uploaded CSV",
            "phenotype": "Live prediction",
            "is_live": True,
        })
        sample_records_for_longitudinal.append({
            "sample_id": sid,
            "patient_id": meta.get("patient_id"),
            "collection_date": meta.get("collection_date"),
            "call": pred["call"],
            "call_label": pred["call_label"],
            "probability": pred["model"]["probability_at_prevalence"],
            "variants": variants,
        })

    from ..longitudinal import analyze_patient_timelines
    timelines = analyze_patient_timelines(sample_records_for_longitudinal)
    if not hasattr(request.app.state, "live_timelines"):
        request.app.state.live_timelines = {}
    request.app.state.live_timelines.update(timelines)

    for sid in rows_by_sample:
        meta = sample_meta.get(sid, {})
        pid = meta.get("patient_id")
        if pid:
            request.app.state.live_results[sid]["patient_id"] = pid
            if pid in request.app.state.live_timelines:
                request.app.state.live_results[sid]["longitudinal_timeline"] = request.app.state.live_timelines[pid]

    for r in new_rows:
        meta = sample_meta.get(r["sample_id"], {})
        pid = meta.get("patient_id")
        if pid:
            r["patient_id"] = pid
            if pid in request.app.state.live_timelines:
                r["n_visits"] = request.app.state.live_timelines[pid]["n_visits"]
                r["has_longitudinal_alerts"] = request.app.state.live_timelines[pid]["has_alerts"]

    from ..label_shift import deduplicate_patient_probabilities, estimate_prevalence_mlls
    from ..surveillance import batch_recurrence_audit

    patient_map = {sid: sample_meta.get(sid, {}).get("patient_id") or sid for sid in rows_by_sample}
    recurrence_alerts = batch_recurrence_audit(rows_by_sample, pred_service.ev, patient_map=patient_map)
    request.app.state.recurrence_alerts = recurrence_alerts

    deduped_probs = deduplicate_patient_probabilities(sample_records_for_longitudinal)
    batch_prev = estimate_prevalence_mlls(deduped_probs, pred_service.pi_tr)
    if batch_prev["status"] == "estimated":
        p_hat = batch_prev["pi_hat"]
        ppv_hat, npv_hat = archived_predictive_values(pred_service.cat_metrics, p_hat)
        batch_prev["ppv_hat"] = ppv_hat
        batch_prev["npv_hat"] = npv_hat
        ppv_low, _ = archived_predictive_values(pred_service.cat_metrics, batch_prev["ci95_low"])
        ppv_high, _ = archived_predictive_values(pred_service.cat_metrics, batch_prev["ci95_high"])
        batch_prev["ppv_low"] = ppv_low
        batch_prev["ppv_high"] = ppv_high
        curve_pts = []
        for x_val in np.linspace(0.005, 0.25, 40):
            px, _ = archived_predictive_values(pred_service.cat_metrics, float(x_val))
            curve_pts.append({"pi": float(x_val), "ppv": float(px)})
        batch_prev["curve_points"] = curve_pts
    request.app.state.batch_prevalence = batch_prev

    if len(rows_by_sample) == 1:
        sid = next(iter(rows_by_sample))
        destination = "/result?" + urlencode({"sample_id": sid})
        return RedirectResponse(destination, status_code=303)

    new_rows.sort(key=lambda r: (r["call"] != "U", -r["probability"], r["sample_id"]))
    for rank, r in enumerate(new_rows, 1):
        r["rank"] = rank if r["call"] == "U" else None

    return render(request, "queue", title="Review queue", active="queue", rows=new_rows, scope="import",
                  upload_quality=quality, batch_prevalence=batch_prev, recurrence_alerts=recurrence_alerts)


@app.get("/result", response_class=HTMLResponse)
def result(request: Request, sample_id: str, prevalence: float = Query(0.05, ge=0.001, le=0.5)) -> Response:
    data = archive(request)
    pred_service = request.app.state.predictor
    upload_quality = getattr(request.app.state, "live_quality", {}).get(sample_id)
    batch_prev = getattr(request.app.state, "batch_prevalence", None)
    prefilled_from_batch = False
    if prevalence == 0.05 and batch_prev and batch_prev.get("status") == "estimated":
        prevalence = batch_prev["pi_hat"]
        prefilled_from_batch = True

    if sample_id in getattr(request.app.state, "live_results", {}):
        record = request.app.state.live_results[sample_id].copy()
        p_train = record["probability_train"]
        p_adj = float(adjust(p_train, pred_service.pi_tr, prevalence))
        record["probability"] = p_adj
        ppv, npv = archived_predictive_values(pred_service.cat_metrics, prevalence)
        return render(request, "result", active="home", title="Bedaquiline result", record=record,
                      prevalence=prevalence, ppv=ppv, npv=npv, upload_quality=upload_quality,
                      prefilled_from_batch=prefilled_from_batch)

    try:
        record = data.replay(sample_id)
    except KeyError:
        return render(request, "error", status=404, title="Recorded isolate unavailable", active="home",
                      detail="Choose an isolate from the replay collection.")
    if "panel_audit" not in record:
        from ..panels import audit_sample_panel
        v_list = [{"gene": item["variant"].split("@")[0], "mutation": item["variant"].split("@")[-1], "frs": float(item.get("frs", 1.0))} for item in record.get("evidence", []) if "@" in item.get("variant", "")]
        record["panel_audit"] = audit_sample_panel(v_list, "wgs")
    cat = data.reports.observed("catalogue")
    ppv, npv = archived_predictive_values(cat, prevalence)
    return render(request, "result", active="home", title="Bedaquiline result", record=record,
                  prevalence=prevalence, ppv=ppv, npv=npv, upload_quality=upload_quality,
                  prefilled_from_batch=prefilled_from_batch)


@app.get("/queue", response_class=HTMLResponse)
def queue(request: Request, scope: str = "all", download: bool = False) -> Response:
    data = archive(request)
    ids = [sid for sid, row in data.roster.items() if scope != "unresolved" or row["call"] == "U"]
    rows = data.queue(ids)
    if download:
        return queue_csv(rows)
    batch_prevalence = getattr(request.app.state, "batch_prevalence", None)
    recurrence_alerts = getattr(request.app.state, "recurrence_alerts", [])
    return render(request, "queue", title="Review queue", active="queue", rows=rows, scope=scope,
                  batch_prevalence=batch_prevalence, recurrence_alerts=recurrence_alerts)


@app.post("/queue", response_class=HTMLResponse)
async def queue_import(request: Request) -> Response:
    data = archive(request)
    pred_service = request.app.state.predictor
    try:
        raw_csv, assay = await form_csv(request)
    except UploadError as error:
        return render(request, "queue", status=error.status, title="Review queue", active="queue", scope="all",
                      rows=data.queue(list(data.roster)), error=str(error))

    parsed = detect_and_parse_upload(raw_csv)
    sample_meta: dict[str, dict] = {}
    if parsed["format"] not in ("discern_csv", "empty"):
        if assay == "wgs" and parsed["inferred_assay"] != "wgs":
            assay = parsed["inferred_assay"]
        rows_by_sample: dict[str, list[dict]] = {}
        for r in parsed["rows"]:
            sid = r["sample_id"]
            sample_meta.setdefault(sid, {})
            if "patient_id" in r and r["patient_id"]:
                sample_meta[sid]["patient_id"] = r["patient_id"]
            if "collection_date" in r and r["collection_date"]:
                sample_meta[sid]["collection_date"] = r["collection_date"]
            rows_by_sample.setdefault(sid, []).append({
                "gene": r["gene"],
                "mutation": r["mutation"],
                "frs": float(r["frs"]),
            })
        for sid in parsed.get("sample_ids", []):
            rows_by_sample.setdefault(sid, [])
        quality = {
            "completeness_score": 100.0,
            "total_rows": len(parsed["rows"]),
            "valid_rows": len(parsed["rows"]),
            "missing_count": 0,
            "duplicate_count": 0,
            "formatting_count": 0,
            "total_issues": 0,
            "issues": [],
            "is_valid": True,
            "summary": f"Imported {len(parsed['rows'])} variants from {parsed['format'].replace('_', ' ').title()}.",
            "format": parsed["format"],
            "headers": ["sample_id", "gene", "mutation", "frs"],
        }
    else:
        try:
            quality = check_csv_quality(raw_csv)
            headers, batch_rows = read_batch(raw_csv)
        except UploadError as error:
            return render(request, "queue", status=error.status, title="Review queue", active="queue", scope="all",
                          rows=data.queue(list(data.roster)), error=str(error))

        first_lines: dict[str, int] = {}
        for line, row in batch_rows:
            first_lines.setdefault(row["sample_id"], line)

        all_known = all(sid in data.roster for sid in first_lines)
        if all_known:
            try:
                ids = data.imported_ids(raw_csv)
                rows = data.queue(ids)
            except UploadError as error:
                return render(request, "queue", status=error.status, title="Review queue", active="queue", scope="all",
                              rows=data.queue(list(data.roster)), error=str(error))
            return render(request, "queue", title="Review queue", active="queue", rows=rows, scope="import", upload_quality=quality)

        # Live scoring for batch
        rows_by_sample = {}
        if headers == PROFILE_COLUMNS or (len(headers) >= 4 and headers[:4] == PROFILE_COLUMNS):
            for line, row in batch_rows:
                sid = row["sample_id"]
                sample_meta.setdefault(sid, {})
                if "patient_id" in row and row["patient_id"]:
                    sample_meta[sid]["patient_id"] = row["patient_id"]
                if "collection_date" in row and row["collection_date"]:
                    sample_meta[sid]["collection_date"] = row["collection_date"]
                try:
                    fraction = Decimal(row["frs"])
                except InvalidOperation:
                    fraction = Decimal("NaN")
                if not fraction.is_finite() or not 0 <= fraction <= 1:
                    return render(
                        request, "queue", status=422, title="Review queue", active="queue", scope="all",
                        rows=data.queue(list(data.roster)),
                        error=f"CSV line {line}, column frs. Enter a number between zero and one.",
                    )
                rows_by_sample.setdefault(sid, []).append({
                    "gene": row["gene"],
                    "mutation": row["mutation"],
                    "frs": float(fraction),
                })
        else:
            for line, row in batch_rows:
                sid = row["sample_id"]
                rows_by_sample.setdefault(row["sample_id"], [])

    new_rows = []
    sample_records_for_longitudinal = []
    for sid, variants in rows_by_sample.items():
        meta = sample_meta.get(sid, {})
        pred = pred_service.predict(variants, prevalence=0.05, sample_id=sid, assay=assay)
        evidence = [{**card, "frs": float(var.get("frs", 1.0))} for card, var in zip(pred["evidence"], variants)]
        if not evidence and pred["evidence"]:
            evidence = pred["evidence"]
        record = {
            "sample_id": sid,
            "call": pred["call"],
            "call_label": pred["call_label"],
            "probability": pred["model"]["probability_at_prevalence"],
            "probability_train": pred["model"]["probability_train_prevalence"],
            "score_raw": pred["model"]["score_raw"],
            "fold": "Live",
            "source": "Uploaded CSV",
            "phenotype": "Live prediction",
            "evidence": evidence,
            "tier_a": pred["tier_a"],
            "clinical_reflex": pred.get("clinical_reflex", []),
            "panel_audit": pred.get("panel_audit"),
            "assay": assay,
            "mode": "Live model evaluation (Discern Tier A + Tier B)",
            "banner": BANNER,
            "score_context": {
                "reference_prevalence": pred_service.pi_tr,
                "uncertainty": None,
                "uncertainty_status": "live",
            },
            "warnings": pred["warnings"],
            "is_live": True,
            "raw_variants": variants,
        }
        request.app.state.live_results[sid] = record
        request.app.state.live_quality[sid] = quality
        new_rows.append({
            "sample_id": sid,
            "call": pred["call"],
            "call_label": pred["call_label"],
            "probability": pred["model"]["probability_at_prevalence"],
            "fold": "Live",
            "source": "Uploaded CSV",
            "phenotype": "Live prediction",
            "is_live": True,
        })
        sample_records_for_longitudinal.append({
            "sample_id": sid,
            "patient_id": meta.get("patient_id"),
            "collection_date": meta.get("collection_date"),
            "call": pred["call"],
            "call_label": pred["call_label"],
            "probability": pred["model"]["probability_at_prevalence"],
            "variants": variants,
        })

    from ..longitudinal import analyze_patient_timelines
    timelines = analyze_patient_timelines(sample_records_for_longitudinal)
    if not hasattr(request.app.state, "live_timelines"):
        request.app.state.live_timelines = {}
    request.app.state.live_timelines.update(timelines)

    for sid in rows_by_sample:
        meta = sample_meta.get(sid, {})
        pid = meta.get("patient_id")
        if pid:
            request.app.state.live_results[sid]["patient_id"] = pid
            if pid in request.app.state.live_timelines:
                request.app.state.live_results[sid]["longitudinal_timeline"] = request.app.state.live_timelines[pid]

    for r in new_rows:
        meta = sample_meta.get(r["sample_id"], {})
        pid = meta.get("patient_id")
        if pid:
            r["patient_id"] = pid
            if pid in request.app.state.live_timelines:
                r["n_visits"] = request.app.state.live_timelines[pid]["n_visits"]
                r["has_longitudinal_alerts"] = request.app.state.live_timelines[pid]["has_alerts"]

    from ..label_shift import deduplicate_patient_probabilities, estimate_prevalence_mlls
    from ..surveillance import batch_recurrence_audit

    patient_map = {sid: sample_meta.get(sid, {}).get("patient_id") or sid for sid in rows_by_sample}
    recurrence_alerts = batch_recurrence_audit(rows_by_sample, pred_service.ev, patient_map=patient_map)
    request.app.state.recurrence_alerts = recurrence_alerts

    deduped_probs = deduplicate_patient_probabilities(sample_records_for_longitudinal)
    batch_prev = estimate_prevalence_mlls(deduped_probs, pred_service.pi_tr)
    if batch_prev["status"] == "estimated":
        p_hat = batch_prev["pi_hat"]
        ppv_hat, npv_hat = archived_predictive_values(pred_service.cat_metrics, p_hat)
        batch_prev["ppv_hat"] = ppv_hat
        batch_prev["npv_hat"] = npv_hat
        ppv_low, _ = archived_predictive_values(pred_service.cat_metrics, batch_prev["ci95_low"])
        ppv_high, _ = archived_predictive_values(pred_service.cat_metrics, batch_prev["ci95_high"])
        batch_prev["ppv_low"] = ppv_low
        batch_prev["ppv_high"] = ppv_high
        curve_pts = []
        for x_val in np.linspace(0.005, 0.25, 40):
            px, _ = archived_predictive_values(pred_service.cat_metrics, float(x_val))
            curve_pts.append({"pi": float(x_val), "ppv": float(px)})
        batch_prev["curve_points"] = curve_pts
    request.app.state.batch_prevalence = batch_prev

    new_rows.sort(key=lambda r: (r["call"] != "U", -r["probability"], r["sample_id"]))
    for rank, r in enumerate(new_rows, 1):
        r["rank"] = rank if r["call"] == "U" else None

    return render(request, "queue", title="Review queue", active="queue", rows=new_rows, scope="import",
                  upload_quality=quality, batch_prevalence=batch_prev, recurrence_alerts=recurrence_alerts)


@app.get("/patient/{patient_id}", response_class=HTMLResponse)
def patient_timeline_view(request: Request, patient_id: str) -> Response:
    timelines = getattr(request.app.state, "live_timelines", {})
    if patient_id not in timelines:
        from ..longitudinal import analyze_patient_timelines
        demo_csv = Path(__file__).resolve().parents[3] / "demo/longitudinal_example.csv"
        if demo_csv.is_file():
            import csv
            by_sid = {}
            meta = {}
            with demo_csv.open(encoding="utf-8") as f:
                reader = csv.DictReader(f)
                for row in reader:
                    sid = row["sample_id"]
                    by_sid.setdefault(sid, []).append({"gene": row["gene"], "mutation": row["mutation"], "frs": float(row["frs"])})
                    meta[sid] = {"patient_id": row["patient_id"], "collection_date": row["collection_date"]}
            demo_rows = []
            for sid, vars in by_sid.items():
                pred = request.app.state.predictor.predict(vars, sample_id=sid)
                demo_rows.append({
                    "sample_id": sid,
                    "patient_id": meta[sid]["patient_id"],
                    "collection_date": meta[sid]["collection_date"],
                    "call": pred["call"],
                    "call_label": pred["call_label"],
                    "probability": pred["model"]["probability_at_prevalence"],
                    "variants": vars,
                })
            timelines.update(analyze_patient_timelines(demo_rows))
            if not hasattr(request.app.state, "live_timelines"):
                request.app.state.live_timelines = {}
            request.app.state.live_timelines.update(timelines)

    if patient_id not in timelines:
        return render(request, "error", status=404, title="Patient timeline unavailable", active="queue",
                      detail=f"No multi-visit longitudinal data recorded for patient {patient_id}.")
    return render(request, "patient_timeline", title=f"Patient Timeline — {patient_id}", active="queue",
                  timeline=timelines[patient_id], patient_id=patient_id)


@app.get("/api/patient/{patient_id}")
def api_patient_timeline(request: Request, patient_id: str) -> dict:
    timelines = getattr(request.app.state, "live_timelines", {})
    if patient_id not in timelines:
        raise HTTPException(404, f"Patient {patient_id} not found in longitudinal registry.")
    return timelines[patient_id]


@app.get("/longitudinal-example.csv")
def longitudinal_example_csv(request: Request) -> Response:
    demo_csv = Path(__file__).resolve().parents[3] / "demo/longitudinal_example.csv"
    if not demo_csv.is_file():
        raise HTTPException(404, "Longitudinal example CSV not found.")
    return csv_response(demo_csv.read_text(encoding="utf-8"), "discern-longitudinal-example.csv")


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
    panels_json_path = Path("results/approaches/panel_audit/e24_panel_coverage/results.json")
    if panels_json_path.is_file():
        with open(panels_json_path, encoding="utf-8") as f:
            panels_data = json.load(f)
    else:
        from ..panels import PANELS, compute_cohort_panel_coverage
        panels_data = {
            "panels": {
                k: {
                    "panel_label": p.label,
                    "platform": p.platform,
                    "vaf_floor": p.vaf_floor,
                    "covered_genes": ",".join(sorted(list(p.covered_genes))),
                    "notes": p.notes,
                    **compute_cohort_panel_coverage(p.covered_genes),
                }
                for k, p in PANELS.items()
            }
        }
    return render(request, "panels", title="Panel coverage", active="panels", panels_data=panels_data)


@app.get("/about", response_class=HTMLResponse)
def about(request: Request) -> Response:
    return render(request, "about", title="About Discern", active="about")
