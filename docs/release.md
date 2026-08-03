# Release checklist

## Local and CI gates

1. Run `python -m unittest discover -s tests -v` on every supported Python
   version and run the configured formatting/static checks.
2. Run `uv build`, inspect both source and wheel archives, and confirm the wheel
   contains `txsuite/resources/project/workflow.schema.json` plus every file in
   all three preset directories.
3. Install the wheel into a clean virtual environment without project
   dependencies, change outside the source checkout, and load the schema and a
   preset through `importlib.resources`.
4. Scaffold each preset, then run the exact user-facing commands:

   ```bash
   txsuite project validate workflow.toml
   txsuite project plan workflow.toml --json
   txsuite project run workflow.toml --dry-run
   ```

5. Exercise a synthetic project run, failed attempt, immutable retry record,
   resume skip, deferred nf-core artifact handoff, and missing/ambiguous
   postflight failure. Confirm persisted snapshots and argv redact secrets.
6. Confirm CI's unit-test, wheel, Nextflow stub, and container workflow results
   are green. Stub and planner tests do not replace the external gates below.

## External release blockers

These checks remain blockers even when the repository CI is green:

1. Build all three owned images and repeat their documented real-container smoke
   analyses.
2. Run the bulk and single-cell nf-core fixtures end to end under both Docker and
   Apptainer, including postflight artifact discovery and one resumed run.
3. Run a public Visium fixture through SpatialData and a licensed Space Ranger
   fixture on a supported Linux host. Record the Space Ranger version and host
   assumptions without publishing licensed binaries or data.
4. Exercise the documented SLURM/Apptainer path on an actual cluster. Site
   account, partition, bind, cache, and registry policy remain external.
5. Treat any failed, unavailable, or skipped external smoke as a blocker, not as
   an implicit pass. Record who ran it, when, the immutable inputs, and where the
   non-sensitive evidence is archived.

## Digest and publication gate

Image tags are convenient during development but are not release locks. Before
tagging, copy the three
`ghcr.io/<owner>/txsuite-<image>@sha256:<digest>` references from the image job
summaries into the release configuration and require
`txsuite env verify-images` to pass. A mutable tag or missing digest is a release
blocker.

Confirm the compatibility table and roadmap match the observed results, update
the package version and release notes, and audit the source archive for secrets
or generated run data. Tag only from a clean working tree after every applicable
gate passes. A `v*` tag publishes owned images; the `containers` workflow can
also be run manually, but manual publication does not waive any gate.
