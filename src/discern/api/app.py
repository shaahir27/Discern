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


async def form_csv(request: Request) -> bytes:
    if request.headers.get("content-type", "").split(";", 1)[0].strip().lower() != "application/x-www-form-urlencoded":
        raise UploadError("unsupported_media_type", "Use the CSV upload form.", status=415)
    # Each CSV byte can expand to three bytes when the browser percent-encodes the form.
    body = await small_body(request, limit=3 * MAX_CSV_BYTES + len("csv="))
    try:
        form = parse_qs(body.decode("utf-8"), max_num_fields=1, strict_parsing=True,
                        keep_blank_values=True, errors="strict")
    except ValueError:
        raise UploadError("invalid_form", "Open the CSV upload form and submit the file again.") from None
    if set(form) != {"csv"}:
        raise UploadError("invalid_form", "Open the CSV upload form and submit the file again.")
    return form["csv"][0].encode("utf-8")


@app.post("/review", response_class=HTMLResponse)
async def review(request: Request) -> Response:
    data = archive(request)
    pred_service = request.app.state.predictor
    try:
        raw_csv = await form_csv(request)
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
    rows_by_sample: dict[str, list[dict]] = {}
    if headers == PROFILE_COLUMNS:
        for line, row in rows:
            sid = row["sample_id"]
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
    for sid, variants in rows_by_sample.items():
        pred = pred_service.predict(variants, prevalence=0.05, sample_id=sid)
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

    if len(rows_by_sample) == 1:
        sid = next(iter(rows_by_sample))
        destination = "/result?" + urlencode({"sample_id": sid})
        return RedirectResponse(destination, status_code=303)

    new_rows.sort(key=lambda r: (r["call"] != "U", -r["probability"], r["sample_id"]))
    for rank, r in enumerate(new_rows, 1):
        r["rank"] = rank if r["call"] == "U" else None

    return render(request, "queue", title="Review queue", active="queue", rows=new_rows, scope="import", upload_quality=quality)


@app.get("/result", response_class=HTMLResponse)
def result(request: Request, sample_id: str, prevalence: float = Query(0.05, ge=0.001, le=0.5)) -> Response:
    data = archive(request)
    pred_service = request.app.state.predictor
    upload_quality = getattr(request.app.state, "live_quality", {}).get(sample_id)

    if sample_id in getattr(request.app.state, "live_results", {}):
        record = request.app.state.live_results[sample_id].copy()
        p_train = record["probability_train"]
        p_adj = float(adjust(p_train, pred_service.pi_tr, prevalence))
        record["probability"] = p_adj
        ppv, npv = archived_predictive_values(pred_service.cat_metrics, prevalence)
        return render(request, "result", active="home", title="Bedaquiline result", record=record,
                      prevalence=prevalence, ppv=ppv, npv=npv, upload_quality=upload_quality)

    try:
        record = data.replay(sample_id)
    except KeyError:
        return render(request, "error", status=404, title="Recorded isolate unavailable", active="home",
                      detail="Choose an isolate from the replay collection.")
    cat = data.reports.observed("catalogue")
    ppv, npv = archived_predictive_values(cat, prevalence)
    return render(request, "result", active="home", title="Bedaquiline result", record=record,
                  prevalence=prevalence, ppv=ppv, npv=npv, upload_quality=upload_quality)


@app.get("/queue", response_class=HTMLResponse)
def queue(request: Request, scope: str = "all", download: bool = False) -> Response:
    data = archive(request)
    ids = [sid for sid, row in data.roster.items() if scope != "unresolved" or row["call"] == "U"]
    rows = data.queue(ids)
    if download:
        return queue_csv(rows)
    return render(request, "queue", title="Review queue", active="queue", rows=rows, scope=scope)


@app.post("/queue", response_class=HTMLResponse)
async def queue_import(request: Request) -> Response:
    data = archive(request)
    pred_service = request.app.state.predictor
    try:
        raw_csv = await form_csv(request)
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
    rows_by_sample: dict[str, list[dict]] = {}
    if headers == PROFILE_COLUMNS:
        for line, row in batch_rows:
            sid = row["sample_id"]
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
            rows_by_sample.setdefault(row["sample_id"], [])

    new_rows = []
    for sid, variants in rows_by_sample.items():
        pred = pred_service.predict(variants, prevalence=0.05, sample_id=sid)
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

    new_rows.sort(key=lambda r: (r["call"] != "U", -r["probability"], r["sample_id"]))
    for rank, r in enumerate(new_rows, 1):
        r["rank"] = rank if r["call"] == "U" else None

    return render(request, "queue", title="Review queue", active="queue", rows=new_rows, scope="import", upload_quality=quality)


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


@app.get("/about", response_class=HTMLResponse)
def about(request: Request) -> Response:
    return render(request, "about", title="About Discern", active="about")
