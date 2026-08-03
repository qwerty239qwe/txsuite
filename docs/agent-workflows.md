# Agent and automation workflows

TxSuite's project commands are designed so people, CI jobs, and coding agents can
inspect a complete plan before authorizing expensive or external execution.

## Safe sequence

Start from a packaged scaffold or a reviewed workflow:

```bash
txsuite project init --preset bulk-rnaseq analysis
cd analysis
txsuite project validate workflow.toml
txsuite project plan workflow.toml --json
txsuite project run workflow.toml --dry-run
txsuite project run workflow.toml --run-id reviewed-run
txsuite project status results/PROJECT_ID/runs/reviewed-run
```

Automation should parse JSON plan/status output rather than terminal-oriented
human plan text. Validation and planning are side-effect free. Dry-run validates
stage selection and prints the plan but does not create a run bundle or launch a
subprocess.

Before `run`, an agent should surface the exact workflow, configuration source,
execution profile, image references, external executables, output root, selected
stages, and unresolved nf-core artifacts for review. A plan is evidence of what
would run; it is not evidence that executables, licensed tools, input datasets,
containers, cluster permissions, or remote registries are available.

## Selection and resume

Use only one selection form at a time:

```bash
txsuite project run workflow.toml --from rnaseq --to differential
txsuite project run workflow.toml --stages rnaseq,differential
```

Selection follows deterministic plan order. When one selected stage fails,
later selected stages are skipped by the default fail-fast executor. Do not infer
success from the existence of a directory; inspect `run-manifest.json` or use
`txsuite project status RUN_DIR`.

Resume always identifies an existing bundle explicitly:

```bash
txsuite project run workflow.toml --resume --run-id reviewed-run
```

Agents must not invent a run ID, delete a bundle, edit stage summaries, or mark a
stage complete. TxSuite decides reuse from hashes and required artifacts. If a
stage must rerun, it appends an immutable attempt rather than rewriting history.

## Deferred artifacts

Plans may display `${producer.artifact}` in downstream argv because nf-core
chooses some result paths at runtime. This is expected. After the producer
completes, TxSuite resolves every declared artifact against the adapter for the
pinned release, records its evidence, and materializes the consumer argv as a
list. Missing or ambiguous matches fail with provenance before the consumer is
started. Prefer a reviewed explicit output override when the result tree contains
multiple valid candidates.

## Security boundaries

- Never place credentials in workflow files or prompts sent to an agent.
- Treat plan and status commands as read-only; `run` launches external tools and
  can consume substantial compute and storage.
- Require separate authorization for licensed tools, remote downloads, registry
  pushes, cluster submissions, or destructive cleanup.
- Use immutable `image@sha256:` references for releases.
- Restrict and review run logs. Redaction covers common argument/config forms but
  cannot sanitize arbitrary tool output or biological data.
- Preserve the original bundle when diagnosing a failure; copy it before any
  manual experimentation.
