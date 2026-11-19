PYTHON ?= python3.12
PYTEST_ARGS ?= -q

.PHONY: verify test

verify:
	"$(PYTHON)" eval/verify.py
	"$(PYTHON)" -m pytest tests/ $(PYTEST_ARGS)

test:
	"$(PYTHON)" -m pytest tests/ $(PYTEST_ARGS)
