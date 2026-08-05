# CELLRANGER_TARBALL must already be present in the build context: TxSuite
# never downloads or bundles licensed 10x Genomics software, so the caller
# supplies their own EULA-accepted cellranger-*.tar.gz (see `env build cellranger`).
FROM ubuntu:24.04@sha256:561618e2c15bf2397621dd04f96926663a3b5616c189cf7e38db7e82f5c538ea

ARG CELLRANGER_TARBALL

RUN apt-get update && apt-get install -y --no-install-recommends tzdata \
    && rm -rf /var/lib/apt/lists/*

COPY ${CELLRANGER_TARBALL} /tmp/cellranger.tar.gz

RUN tar -xzf /tmp/cellranger.tar.gz -C /opt \
    && mv /opt/cellranger-*/ /opt/cellranger \
    && rm /tmp/cellranger.tar.gz

ENV PATH="/opt/cellranger:${PATH}"

CMD ["cellranger", "--help"]
