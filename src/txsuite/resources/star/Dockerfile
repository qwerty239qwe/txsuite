FROM condaforge/miniforge3:26.3.2-3@sha256:532f6ee7a858b009dc895f8313eb6ed875f05a455ec57d832aa6b4c66e2799b9

# Reference preparation for splice-aware alignment and GATK short variant
# calling. The three tools ship together because they build one artifact set
# from the same FASTA and GTF: STAR needs the index, GATK needs the sequence
# dictionary, and both need the FASTA index.
RUN conda install --yes --name base --channel conda-forge --channel bioconda \
        "star=2.7.11b" \
        "samtools=1.24" \
        "gatk4=4.6.2.0" \
    && conda clean --all --yes \
    && STAR --version \
    && samtools --version | head -1 \
    && gatk --version

COPY merge_star_counts.py /opt/txsuite/merge_star_counts.py

WORKDIR /work
CMD ["STAR", "--help"]
