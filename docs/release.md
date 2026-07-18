# Release checklist

1. Run `python -m unittest discover -s tests`, Ruff, and `uv build`.
2. Build all three owned images and repeat their documented smoke analyses.
3. Run the bulk and single-cell nf-core fixtures under Docker and Apptainer.
4. Run a public Visium fixture through SpatialData and a licensed Space Ranger
   fixture on a supported Linux host.
5. Confirm the compatibility table and roadmap reflect the observed results.
6. Tag only after the working tree is clean and CI passes. A `v*` tag publishes
   the owned images; the `containers` workflow can also be run manually.
7. Copy the three `ghcr.io/<owner>/txsuite-<image>@sha256:<digest>` references
   from the job summaries into `txsuite.toml`, then require
   `txsuite env verify-images` to pass.

Image tags are convenient during development but are not release locks. Docker
digest references are immutable; TxSuite's verification command intentionally
returns nonzero while any configured image is tag-only.
