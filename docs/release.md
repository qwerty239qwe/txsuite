# Release checklist

1. Run `python -m unittest discover -s tests`, Ruff, and `uv build`.
2. Build all three owned images and repeat their documented smoke analyses.
3. Push owned images to a registry, replace mutable tags in `txsuite.toml` with
   `name@sha256:digest`, and require `txsuite env verify-images` to pass.
4. Run the bulk and single-cell nf-core fixtures under Docker and Apptainer.
5. Run a public Visium fixture through SpatialData and a licensed Space Ranger
   fixture on a supported Linux host.
6. Confirm the compatibility table and roadmap reflect the observed results.
7. Tag only after the working tree is clean and CI passes.

Image tags are convenient during development but are not release locks. Docker
digest references are immutable; TxSuite's verification command intentionally
returns nonzero while any configured image is tag-only.
