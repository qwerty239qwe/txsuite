FROM bioconductor/bioconductor_docker:RELEASE_3_23@sha256:1d871e1ca9cca76b220eb16e22677e728f4352f81a9ee91aaf29e24aea43e624

RUN R -e 'BiocManager::install("DESeq2", ask = FALSE, update = FALSE)' \
    && R -e 'stopifnot(requireNamespace("DESeq2", quietly = TRUE))'

COPY deseq2.R /opt/txsuite/deseq2.R
ENTRYPOINT ["Rscript", "/opt/txsuite/deseq2.R"]
