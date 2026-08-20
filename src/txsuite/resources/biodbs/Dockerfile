FROM python:3.12-slim-bookworm@sha256:d50fb7611f86d04a3b0471b46d7557818d88983fc3136726336b2a4c657aa30b

ENV PYTHONDONTWRITEBYTECODE=1 PYTHONUNBUFFERED=1

# biodbs is a REST client for KEGG, QuickGO, Reactome and BioMart. Pinned
# exactly because it is pre-1.0 and its API may still move; the stage that uses
# it records the resolved version in its provenance file.
RUN python -m pip install --no-cache-dir "biodbs==0.4.1"

COPY fetch_genesets.py /opt/txsuite/fetch_genesets.py

WORKDIR /work
CMD ["python", "/opt/txsuite/fetch_genesets.py", "--help"]
