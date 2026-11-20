FROM python:3.12-slim
WORKDIR /app
ENV PYTHONUNBUFFERED=1 PYTHONPATH=/app/src DISCERN_ROOT=/app
COPY requirements.txt .
RUN grep -E '^(pandas|numpy|pyarrow|fastapi|uvicorn|jinja2)>=' requirements.txt > /tmp/runtime-requirements.txt \
    && pip install --no-cache-dir -r /tmp/runtime-requirements.txt
RUN grep -E '^(scikit-learn|scipy|joblib|piezo|pyyaml)[=<>]' requirements.txt > /tmp/scoring-requirements.txt \
    && pip install --no-cache-dir --timeout 60 -r /tmp/scoring-requirements.txt
COPY src ./src
COPY config.yaml .
COPY demo/held_out_samples.csv demo/longitudinal_example.csv ./demo/
COPY models/lr_model ./models/lr_model
COPY models/lr_mechanism/lof_effects.csv ./models/lr_mechanism/
COPY results ./results
COPY docs ./docs
RUN adduser --disabled-password --gecos "" discern
USER discern
EXPOSE 8000
CMD ["python", "-m", "uvicorn", "discern.api.app:app", "--host", "0.0.0.0", "--port", "8000", "--no-access-log"]
