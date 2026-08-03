# Project run bundles and results

Project runs live below:

```text
<output_root>/<project.id>/runs/<run-id>/
├── workflow.json
├── resolved-config.json
├── command-plan.json
├── run-manifest.json
├── logs/
│   ├── <stage>.attempt-1.stdout.log
│   └── <stage>.attempt-1.stderr.log
└── stage-results/
    ├── <stage>.json
    └── <stage>.attempt-1.<run-id>.json
```

The run ID combines a UTC timestamp with random entropy unless supplied with
`--run-id`. Bundle creation is exclusive: an existing run directory is never
silently replaced.

## Immutable and mutable records

`workflow.json`, `resolved-config.json`, and `command-plan.json` are write-once
snapshots. `run-manifest.json` is atomically replaced as run and stage states
change. `<stage>.json` is the current stage summary. Every finalized attempt has
a distinct, write-once record containing its timing, redacted argv, exit code,
output checks, artifact evidence, and log paths. Retrying a stage creates attempt
2 and leaves attempt 1 byte-for-byte intact.

Stage states are `pending`, `running`, `completed`, `failed`, and `skipped`.
Subprocess and postflight failures are persisted before an optional CLI error is
raised. Standard output and standard error go directly to the attempt's
`stdout.log` and `stderr.log`; commands are executed as argument lists with
`shell=False`.

The manifest indexes current stage states, attempt-record paths, canonical input
hashes, and resolved artifacts. nf-core artifact entries contain the absolute
path plus adapter, pipeline, release, match pattern, and whether an explicit
override or adapter rule selected the path. Concrete planned outputs enter the
same artifact map after existence checks.

## Resume

Resume is explicit and requires the original ID:

```bash
txsuite project run workflow.toml --resume --run-id RUN_ID
```

A stage is skipped only when its prior state is reusable, its input and resolved
configuration hashes match, its fully materialized command and parameters match,
and all required outputs and recorded artifacts still exist. Changed commands,
configuration, declared inputs, or missing artifacts create a new attempt. Stage
selection does not broaden resume eligibility.

## Redaction and security

Persisted workflow, resolved configuration, command plans, and attempt argv are
redacted recursively for common credential, password, secret, authorization,
private-key, API-key, and token names. Split arguments such as `--token VALUE`,
inline arguments such as `--password=VALUE`, environment-style assignments, and
URL user information are also redacted.

Redaction is defense in depth, not a credential store. Do not put secrets in
`workflow.toml`, command-line arguments, filenames, output data, or logs. A
canonical SHA-256 of the original snapshot is retained for change detection;
hashes can reveal equality and do not protect low-entropy secrets from guessing.
Use the platform's environment or secret manager, restrict filesystem access to
run bundles, and review logs before sharing or archiving them.

The older single-command CLI writes `run.json` and `txsuite-results.json` in its
task `.txsuite` directory. That legacy layout is intentionally separate from the
multi-stage project run bundle described here.
