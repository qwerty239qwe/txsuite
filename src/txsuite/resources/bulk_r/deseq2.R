suppressPackageStartupMessages(library(DESeq2))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 6) {
    stop("usage: deseq2.R COUNTS METADATA DESIGN REFERENCE TEST OUTDIR")
}

counts_path <- args[[1]]
metadata_path <- args[[2]]
design <- args[[3]]
reference <- args[[4]]
test <- args[[5]]
outdir <- args[[6]]

counts <- read.delim(counts_path, check.names = FALSE, stringsAsFactors = FALSE)
metadata <- read.delim(metadata_path, check.names = FALSE, stringsAsFactors = FALSE)
if (ncol(counts) < 2 || ncol(metadata) < 2) {
    stop("counts and metadata must each contain an ID column plus data columns")
}

gene_ids <- counts[[1]]
if (anyNA(gene_ids) || anyDuplicated(gene_ids)) {
    stop("the first counts column must contain unique gene IDs")
}
count_matrix <- as.matrix(counts[-1])
suppressWarnings(storage.mode(count_matrix) <- "numeric")
if (anyNA(count_matrix) || any(count_matrix < 0) || any(count_matrix != round(count_matrix))) {
    stop("counts must be non-negative integers")
}
storage.mode(count_matrix) <- "integer"
rownames(count_matrix) <- gene_ids

sample_ids <- metadata[[1]]
if (anyNA(sample_ids) || anyDuplicated(sample_ids)) {
    stop("the first metadata column must contain unique sample IDs")
}
rownames(metadata) <- sample_ids
missing_samples <- setdiff(colnames(count_matrix), rownames(metadata))
if (length(missing_samples)) {
    stop(paste("metadata is missing samples:", paste(missing_samples, collapse = ", ")))
}
metadata <- metadata[colnames(count_matrix), , drop = FALSE]
if (!(design %in% colnames(metadata))) {
    stop(paste("metadata has no design column:", design))
}

metadata[[design]] <- factor(metadata[[design]])
levels_present <- levels(metadata[[design]])
if (!(reference %in% levels_present) || !(test %in% levels_present)) {
    stop("reference and test must both occur in the design column")
}
metadata[[design]] <- relevel(metadata[[design]], ref = reference)

dds <- DESeqDataSetFromMatrix(
    countData = count_matrix,
    colData = metadata,
    design = reformulate(design)
)
dds <- DESeq(dds)
result <- as.data.frame(results(dds, contrast = c(design, test, reference)))
result <- data.frame(gene_id = rownames(result), result, check.names = FALSE)

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
write.table(result, file.path(outdir, "deseq2-results.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
normalized <- data.frame(gene_id = rownames(dds), counts(dds, normalized = TRUE), check.names = FALSE)
write.table(
    normalized,
    file.path(outdir, "normalized-counts.tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)
writeLines(capture.output(sessionInfo()), file.path(outdir, "session-info.txt"))

