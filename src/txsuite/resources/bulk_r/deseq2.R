suppressPackageStartupMessages(library(DESeq2))

args <- commandArgs(trailingOnly = TRUE)
if (!(length(args) %in% c(6, 10))) {
    stop(paste(
        "usage: deseq2.R COUNTS METADATA DESIGN REFERENCE TEST OUTDIR",
        "[PADJ LFC TOP_GENES COVARIATES]"
    ))
}

counts_path <- args[[1]]
metadata_path <- args[[2]]
design <- args[[3]]
reference <- args[[4]]
test <- args[[5]]
outdir <- args[[6]]
padj_threshold <- if (length(args) == 10) as.numeric(args[[7]]) else 0.05
lfc_threshold <- if (length(args) == 10) as.numeric(args[[8]]) else 1
top_genes <- if (length(args) == 10) as.integer(args[[9]]) else 50L
covariates <- if (length(args) == 10 && nzchar(args[[10]])) {
    strsplit(args[[10]], ",", fixed = TRUE)[[1]]
} else {
    character()
}

if (!grepl("^[A-Za-z][A-Za-z0-9_.]*$", design)) {
    stop("design must be a simple metadata column name")
}
if (any(!grepl("^[A-Za-z][A-Za-z0-9_.]*$", covariates))) {
    stop("covariates must be simple metadata column names")
}
if (design %in% covariates || anyDuplicated(covariates)) {
    stop("covariates must be unique and different from the design column")
}
if (!is.finite(padj_threshold) || padj_threshold <= 0 || padj_threshold > 1) {
    stop("adjusted p-value threshold must be in (0, 1]")
}
if (!is.finite(lfc_threshold) || lfc_threshold < 0) {
    stop("absolute log2 fold-change threshold must be non-negative")
}
if (is.na(top_genes) || top_genes < 1) {
    stop("top genes must be positive")
}

counts <- read.delim(counts_path, check.names = FALSE, stringsAsFactors = FALSE)
metadata <- read.delim(metadata_path, check.names = FALSE, stringsAsFactors = FALSE)
if (ncol(counts) < 2 || ncol(metadata) < 2) {
    stop("counts and metadata must each contain an ID column plus data columns")
}

gene_ids <- as.character(counts[[1]])
if (anyNA(gene_ids) || any(!nzchar(gene_ids)) || anyDuplicated(gene_ids)) {
    stop("the first counts column must contain unique, non-empty gene IDs")
}
count_matrix <- as.matrix(counts[-1])
if (any(!nzchar(colnames(count_matrix))) || anyDuplicated(colnames(count_matrix))) {
    stop("count sample columns must have unique, non-empty names")
}
suppressWarnings(storage.mode(count_matrix) <- "numeric")
if (anyNA(count_matrix) || any(count_matrix < 0) || any(count_matrix != round(count_matrix))) {
    stop("counts must be non-negative integers")
}
storage.mode(count_matrix) <- "integer"
rownames(count_matrix) <- gene_ids
count_matrix <- count_matrix[rowSums(count_matrix) > 0, , drop = FALSE]
if (!nrow(count_matrix)) {
    stop("counts contain no expressed genes")
}

sample_ids <- as.character(metadata[[1]])
if (anyNA(sample_ids) || any(!nzchar(sample_ids)) || anyDuplicated(sample_ids)) {
    stop("the first metadata column must contain unique, non-empty sample IDs")
}
rownames(metadata) <- sample_ids
missing_samples <- setdiff(colnames(count_matrix), rownames(metadata))
if (length(missing_samples)) {
    stop(paste("metadata is missing samples:", paste(missing_samples, collapse = ", ")))
}
metadata <- metadata[colnames(count_matrix), , drop = FALSE]
model_columns <- c(covariates, design)
missing_columns <- setdiff(model_columns, colnames(metadata))
if (length(missing_columns)) {
    stop(paste("metadata is missing model columns:", paste(missing_columns, collapse = ", ")))
}
for (column in model_columns) {
    if (is.character(metadata[[column]])) {
        metadata[[column]] <- factor(metadata[[column]])
    }
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
    design = reformulate(model_columns)
)
dds <- DESeq(dds)
de_result <- results(
    dds,
    contrast = c(design, test, reference),
    alpha = padj_threshold
)
result <- data.frame(
    gene_id = rownames(de_result),
    as.data.frame(de_result),
    check.names = FALSE
)
significant <- result[
    !is.na(result$padj) &
        result$padj <= padj_threshold &
        abs(result$log2FoldChange) >= lfc_threshold,
    ,
    drop = FALSE
]

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
write.table(
    result,
    file.path(outdir, "deseq2-results.tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)
write.table(
    significant,
    file.path(outdir, "significant-genes.tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)

normalized_matrix <- counts(dds, normalized = TRUE)
normalized <- data.frame(
    gene_id = rownames(normalized_matrix),
    normalized_matrix,
    check.names = FALSE
)
write.table(
    normalized,
    file.path(outdir, "normalized-counts.tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)

vst_matrix <- assay(varianceStabilizingTransformation(dds, blind = TRUE))
write.table(
    data.frame(gene_id = rownames(vst_matrix), vst_matrix, check.names = FALSE),
    file.path(outdir, "vst-counts.tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)

sample_qc <- data.frame(
    sample = colnames(count_matrix),
    library_size = colSums(count_matrix),
    detected_genes = colSums(count_matrix > 0),
    metadata[colnames(count_matrix), model_columns, drop = FALSE],
    check.names = FALSE,
    row.names = NULL
)
write.table(
    sample_qc,
    file.path(outdir, "sample-qc.tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)

pca <- prcomp(t(vst_matrix))
pca_scores <- pca$x
if (ncol(pca_scores) < 2) {
    pca_scores <- cbind(pca_scores, PC2 = 0)
}
pca_table <- data.frame(
    sample = rownames(pca_scores),
    pca_scores[, 1:2, drop = FALSE],
    metadata[rownames(pca_scores), model_columns, drop = FALSE],
    check.names = FALSE,
    row.names = NULL
)
write.table(
    pca_table,
    file.path(outdir, "pca.tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)

sample_correlation <- cor(vst_matrix, method = "pearson")
write.table(
    data.frame(sample = rownames(sample_correlation), sample_correlation, check.names = FALSE),
    file.path(outdir, "sample-correlation.tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)

pdf(file.path(outdir, "pca.pdf"), width = 7, height = 6)
groups <- metadata[rownames(pca_scores), design]
group_factor <- factor(groups)
variance <- 100 * pca$sdev^2 / sum(pca$sdev^2)
if (length(variance) < 2) {
    variance <- c(variance, 0)
}
plot(
    pca_scores[, 1],
    pca_scores[, 2],
    col = as.integer(group_factor),
    pch = 19,
    xlab = sprintf("PC1 (%.1f%%)", variance[[1]]),
    ylab = sprintf("PC2 (%.1f%%)", variance[[2]]),
    main = paste(test, "vs", reference)
)
text(pca_scores[, 1], pca_scores[, 2], labels = rownames(pca_scores), pos = 3, cex = 0.7)
legend("topright", legend = levels(group_factor), col = seq_along(levels(group_factor)), pch = 19)
dev.off()

pdf(file.path(outdir, "sample-correlation.pdf"), width = 8, height = 8)
heatmap(sample_correlation, scale = "none", symm = TRUE, margins = c(9, 9))
dev.off()

pdf(file.path(outdir, "ma.pdf"), width = 7, height = 6)
plotMA(de_result, alpha = padj_threshold, main = paste(test, "vs", reference))
abline(h = c(-lfc_threshold, lfc_threshold), col = "steelblue", lty = 2)
dev.off()

pdf(file.path(outdir, "volcano.pdf"), width = 7, height = 6)
volcano_y <- -log10(pmax(result$padj, .Machine$double.xmin))
is_significant <- !is.na(result$padj) &
    result$padj <= padj_threshold &
    abs(result$log2FoldChange) >= lfc_threshold
plot(
    result$log2FoldChange,
    volcano_y,
    col = ifelse(is_significant, "firebrick", "grey70"),
    pch = 19,
    cex = 0.6,
    xlab = "log2 fold change",
    ylab = "-log10 adjusted p-value",
    main = paste(test, "vs", reference)
)
abline(v = c(-lfc_threshold, lfc_threshold), h = -log10(padj_threshold), lty = 2)
label_rows <- head(order(result$padj, na.last = NA), 10L)
text(
    result$log2FoldChange[label_rows],
    volcano_y[label_rows],
    labels = result$gene_id[label_rows],
    pos = 3,
    cex = 0.6
)
dev.off()

pdf(file.path(outdir, "top-genes-heatmap.pdf"), width = 9, height = 9)
if (nrow(vst_matrix) >= 2 && ncol(vst_matrix) >= 2) {
    gene_variance <- apply(vst_matrix, 1, var)
    selected <- head(order(gene_variance, decreasing = TRUE), min(top_genes, nrow(vst_matrix)))
    heatmap(
        vst_matrix[selected, , drop = FALSE],
        scale = "row",
        margins = c(9, 9),
        main = "Top variable genes"
    )
} else {
    plot.new()
    text(0.5, 0.5, "At least two genes and samples are required")
}
dev.off()

write.table(
    data.frame(
        metric = c("genes_tested", "significant_genes", "samples", "padj", "abs_log2fc"),
        value = c(nrow(result), nrow(significant), ncol(count_matrix), padj_threshold, lfc_threshold)
    ),
    file.path(outdir, "analysis-summary.tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)
writeLines(capture.output(sessionInfo()), file.path(outdir, "session-info.txt"))
