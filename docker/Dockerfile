FROM python:3.12-slim

LABEL org.opencontainers.image.title="CellFlowX"
LABEL org.opencontainers.image.description="Reproducible single-cell RNA-seq analysis environment for CellFlowX"

ENV PYTHONDONTWRITEBYTECODE=1
ENV PYTHONUNBUFFERED=1
ENV MPLBACKEND=Agg

WORKDIR /app

RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    git \
    curl \
    procps \
    && rm -rf /var/lib/apt/lists/*

COPY requirements.txt /app/requirements.txt

RUN pip install --no-cache-dir --upgrade pip && \
    pip install --no-cache-dir -r /app/requirements.txt

COPY scripts /app/scripts
COPY config /app/config
COPY manifests /app/manifests
COPY streamlit_app /app/streamlit_app
COPY database/schema.sql /app/database/schema.sql

CMD ["python", "--version"]
