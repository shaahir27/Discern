# Discern

Discern is an offline research explorer for archived evaluations of bedaquiline resistance in tuberculosis. Open a recorded isolate, compare its saved catalogue call with the laboratory phenotype, and inspect variant evidence and its uncertainty.

**Research use only. Clinical decisions require independent laboratory confirmation.**

## Run locally

From the repository root, using the Discern Python environment:

```bash
python -m pip install -r requirements.txt
python -m pip install -e . --no-deps
python -m uvicorn discern.api.app:app --host 127.0.0.1 --port 8000 --no-access-log
```

Open [the local explorer](http://127.0.0.1:8000/). On **Check a sample**, upload the example CSV or open a saved example directly. The result shows the catalogue call, saved model score, and lab finding. **Review queue** organizes recorded cases and **Study results** contains the evaluation tables.

The sample CSV format is `sample_id,gene,mutation,frs`. A CSV with only the `sample_id` column also works. Uploaded variant rows must match the full profile of a recorded sample. New sequences and changed profiles are outside the explorer’s scope.

## CLI and API

Use a recorded sample ID from the roster available in Review queue.

```bash
discern replay --sample-id YOUR_RECORDED_ID
discern batch discern-example.csv ranked.csv
discern verify
```

The explorer API offers `GET /health`, `GET /replay/{sample_id}`, `GET /example.csv`, `GET /replay-samples.csv`, `GET /evidence/{gene}/{mutation}`, and `POST /batch` with a `text/csv` body. Batches retrieve saved results for the bundled cohort.

## Docker

```bash
docker build -t discern-research .
docker run --rm -p 127.0.0.1:8000:8000 discern-research
```

The built container reads local evidence and saved reports at runtime. It needs no network connection for the archive, interface assets, or documentation.

## Reproduce and verify

Run `make verify` from the repository root. The development evaluation command is `discern eval --all` and requires the local evaluation data. Regenerate the public documents with `python eval/generate_docs.py`.

[Model card](docs/MODEL_CARD.md) · [Results](docs/RESULTS.md) · [Demo](docs/DEMO_SCRIPT.md) · [Glossary](docs/GLOSSARY.md) · [Backend review](docs/BACKEND_REVIEW.md)
