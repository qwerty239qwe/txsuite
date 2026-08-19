FROM condaforge/miniforge3:26.3.2-3@sha256:532f6ee7a858b009dc895f8313eb6ed875f05a455ec57d832aa6b4c66e2799b9

# Selective-alignment quantification for bulk (salmon) and single-cell
# (simpleaf, which drives piscem/salmon plus alevin-fry). gffread derives the
# transcriptome FASTA and the transcript-to-gene table from the annotation, so
# a gene-level matrix needs no input beyond --fasta and --gtf.
#
# miniforge3 keeps the base environment on PATH and starts no wrapper process of
# its own, so every process in the native DAGs can pass an explicit command,
# matching the other TxSuite-owned images.
RUN conda install --yes --name base --channel conda-forge --channel bioconda \
        "salmon=2.5.1" \
        "simpleaf=0.28.0" \
        "alevin-fry=0.18.0" \
        "piscem=0.22.1" \
        "gffread=0.12.9" \
    && conda clean --all --yes \
    && salmon --version \
    && simpleaf --version \
    && alevin-fry --version \
    && gffread --version

COPY merge_quants.py /opt/txsuite/merge_quants.py

# simpleaf writes its configuration under this directory; keeping it writable
# and outside the mounted work directory avoids polluting published results.
ENV ALEVIN_FRY_HOME=/tmp/alevin-fry-home

WORKDIR /work
CMD ["salmon", "--help"]
