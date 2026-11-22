# Discern

Discern is an offline research explorer for archived evaluations of bedaquiline resistance in tuberculosis. Open a recorded isolate, compare its saved catalogue call with the laboratory phenotype, and inspect variant evidence and its uncertainty.

**Research use only. Clinical decisions require independent laboratory confirmation.**

## Run locally

From the repository root, using the Discern Python environment:

```bash
# 1. Install dependencies
python -m pip install -r requirements.txt
python -m pip install -e . --no-deps

# 2. Launch using the shortcut:
make app

# Or run directly:
python -m uvicorn discern.api.app:app --app-dir src --host 127.0.0.1 --port 8000 --no-access-log
```

Open [the local explorer](http://127.0.0.1:8000/). On **Check a sample**, upload the example CSV or open a saved example directly. The result shows the catalogue call, saved model score, and lab finding. **Review queue** organizes recorded cases and **Study results** contains the evaluation tables.

On Windows PowerShell, use the included shortcut script:

```powershell
conda activate discern
.\serve.ps1

# Or run directly:
python -m uvicorn discern.api.app:app --app-dir src --host 127.0.0.1 --port 8000 --no-access-log
```

Use `sample_id,gene,mutation,frs` for a Discern CSV. The upload forms also accept supported TB-Profiler JSON, Deeplex CSV, and AmPORE-TB CSV or TSV exports. A file containing only recorded sample IDs opens archived results. A new sample ID opens a live research result. Add `patient_id` and `collection_date` to connect dated visits in an uploaded batch. The Review queue includes an illustrative timeline.

## CLI and API

Use a recorded sample ID from the roster available in Review queue.

```bash
discern replay --sample-id YOUR_RECORDED_ID
discern batch discern-example.csv ranked.csv
discern verify
```

The explorer API offers `GET /health`, `POST /predict` (single-sample JSON prediction with Pydantic validation), `GET /replay/{sample_id}`, `GET /example.csv`, `GET /replay-samples.csv`, `GET /evidence/{gene}/{mutation}`, and `POST /batch` (prioritized triage queue CSV).

## Docker

```bash
docker build -t discern-research .
docker run --rm -p 127.0.0.1:8000:8000 discern-research
```

The built container reads local evidence and saved reports at runtime. It needs no network connection for the archive, interface assets, or documentation.

## Reproduce and verify

Run `make verify` from the repository root. The development evaluation command is `discern eval --all` and requires the local evaluation data. Regenerate the public documents with `python eval/generate_docs.py`.

[Codebase guide](CODEBASE_GUIDE.md) · [Model card](docs/MODEL_CARD.md) · [Results](docs/RESULTS.md) · [Demo](docs/DEMO_SCRIPT.md) · [Glossary](docs/GLOSSARY.md) · [Backend review](docs/BACKEND_REVIEW.md)

[Final novelty approach](docs/FINAL_NOVELTY_APPROACH.md) · [Novelty upgrades](docs/NOVELTY_UPGRADES.md)
