FROM python:3.12-slim-bookworm@sha256:d50fb7611f86d04a3b0471b46d7557818d88983fc3136726336b2a4c657aa30b

ENV PYTHONDONTWRITEBYTECODE=1 PYTHONUNBUFFERED=1

RUN python -m pip install --no-cache-dir \
    "spatialdata==0.8.0" \
    "spatialdata-io==0.7.1" \
    "squidpy[leiden]==1.8.3"

COPY spatial.py /opt/txsuite/spatial.py

ENTRYPOINT ["python", "/opt/txsuite/spatial.py"]
