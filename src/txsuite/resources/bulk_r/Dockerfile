FROM bioconductor/bioconductor_docker:RELEASE_3_23@sha256:1d871e1ca9cca76b220eb16e22677e728f4352f81a9ee91aaf29e24aea43e624

RUN R -e 'BiocManager::install(c("DESeq2", "edgeR", "limma", "clusterProfiler"), ask = FALSE, update = FALSE)' \
    && R -e 'stopifnot(all(vapply(c("DESeq2", "edgeR", "limma", "clusterProfiler"), requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))))'

COPY deseq2.R alternative_de.R enrichment.R /opt/txsuite/
CMD ["Rscript", "/opt/txsuite/deseq2.R"]
