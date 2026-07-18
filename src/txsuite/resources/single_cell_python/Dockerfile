FROM python:3.12-slim-bookworm@sha256:d50fb7611f86d04a3b0471b46d7557818d88983fc3136726336b2a4c657aa30b

ENV PYTHONDONTWRITEBYTECODE=1 PYTHONUNBUFFERED=1

RUN python -m pip install --no-cache-dir "scanpy[leiden]==1.12.2"

COPY single_cell.py /opt/txsuite/single_cell.py

CMD ["python", "/opt/txsuite/single_cell.py", "--help"]
