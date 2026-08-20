# Bulk RNA-seq alignment project

Splice-aware genome alignment with STAR, then differential expression. The
`ref` stage builds the STAR index once and `align` reuses it through
`${ref.star_index}`; the index is required, so the expensive artifact is never
rebuilt silently.

STAR emits gene counts as a by-product of alignment, so `align` produces both
sorted BAMs and a count matrix that feeds `bulk.de` directly — no separate
counting stage.

Replace every `data/...` path with real inputs. Keep `samplesheet.csv` sample
IDs identical to the `sample` column in `metadata.tsv`.

```console
txsuite project validate workflow.toml
txsuite project plan workflow.toml
txsuite project run workflow.toml
```

Run these commands from this directory.

`strandedness = "auto"` infers the library orientation from the balance between
STAR's forward and reverse count columns and records the evidence in
`counts/strandedness.tsv`. Check that file after the first run: an inference
that disagrees with your library prep means something is wrong upstream. Pin it
to `unstranded`, `forward`, or `reverse` when you already know. The run fails
rather than guessing if the split is unclear, or if samples disagree with each
other.

Set `read_length` in the `ref` stage to your actual read length; it drives the
splice-junction overhang.
