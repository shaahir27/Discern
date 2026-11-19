FROM python:3.12-slim
WORKDIR /app
ENV PYTHONUNBUFFERED=1 PYTHONPATH=/app/src DISCERN_ROOT=/app
COPY requirements.txt .
RUN grep -E '^(pandas|numpy|pyarrow|fastapi|uvicorn|jinja2)>=' requirements.txt > /tmp/runtime-requirements.txt \
    && pip install --no-cache-dir -r /tmp/runtime-requirements.txt
COPY src ./src
COPY config.yaml .
COPY demo/held_out_samples.csv ./demo/held_out_samples.csv
COPY models/lr_model/evidence_store.csv models/lr_model/metadata.json ./models/lr_model/
COPY results ./results
COPY docs ./docs
RUN adduser --disabled-password --gecos "" discern
USER discern
EXPOSE 8000
CMD ["python", "-m", "uvicorn", "discern.api.app:app", "--host", "0.0.0.0", "--port", "8000", "--no-access-log"]
