PYTHON ?= python3.12
PYTEST_ARGS ?= -q

.PHONY: verify test docs app

verify:
	"$(PYTHON)" eval/verify.py
	"$(PYTHON)" -m pytest tests/ $(PYTEST_ARGS)

test:
	"$(PYTHON)" -m pytest tests/ $(PYTEST_ARGS)

docs:
	"$(PYTHON)" eval/generate_docs.py

app:
	"$(PYTHON)" -m uvicorn discern.api.app:app --app-dir src --host 127.0.0.1 --port 8000 --no-access-log
