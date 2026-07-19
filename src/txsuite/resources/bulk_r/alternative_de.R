suppressPackageStartupMessages({
    library(edgeR)
    library(limma)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 11) {
    stop(paste(
        "usage: alternative_de.R METHOD COUNTS METADATA DESIGN REFERENCE TEST",
        "OUTDIR PADJ LFC TOP_GENES COVARIATES"
    ))
}

method <- args[[1]]
counts_path <- args[[2]]
metadata_path <- args[[3]]
design <- args[[4]]
reference <- args[[5]]
test <- args[[6]]
outdir <- args[[7]]
padj_threshold <- as.numeric(args[[8]])
lfc_threshold <- as.numeric(args[[9]])
top_genes <- as.integer(args[[10]])
covariates <- if (nzchar(args[[11]])) {
    strsplit(args[[11]], ",", fixed = TRUE)[[1]]
} else {
    character()
}

if (!(method %in% c("edger", "limma"))) {
    stop("method must be 'edger' or 'limma'")
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
suppressWarnings(storage.mode(count_matrix) <- "numeric")
if (anyNA(count_matrix) || any(count_matrix < 0) || any(count_matrix != round(count_matrix))) {
    stop("counts must be non-negative integers")
}
storage.mode(count_matrix) <- "integer"
rownames(count_matrix) <- gene_ids

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
if (anyNA(metadata[, model_columns, drop = FALSE])) {
    stop("model columns cannot contain missing values")
}

selected <- metadata[[design]] %in% c(reference, test)
metadata <- droplevels(metadata[selected, , drop = FALSE])
count_matrix <- count_matrix[, rownames(metadata), drop = FALSE]
group <- factor(metadata[[design]], levels = c(reference, test))
if (anyNA(group) || any(table(group) < 2)) {
    stop("reference and test must each contain at least two samples")
}
for (column in covariates) {
    if (is.character(metadata[[column]])) {
        metadata[[column]] <- factor(metadata[[column]])
    }
}
metadata$.txsuite_test <- as.integer(group == test)
model <- model.matrix(reformulate(c(covariates, ".txsuite_test")), metadata)
if (qr(model)$rank < ncol(model)) {
    stop("the design matrix is not full rank; check covariates and groups")
}
coefficient <- match(".txsuite_test", colnames(model))

dge <- DGEList(counts = count_matrix)
keep <- filterByExpr(dge, design = model)
if (!any(keep)) {
    stop("no genes pass expression filtering")
}
dge <- calcNormFactors(dge[keep, , keep.lib.sizes = FALSE])
normalized_matrix <- cpm(dge, normalized.lib.sizes = TRUE, log = FALSE)
log_cpm <- cpm(dge, normalized.lib.sizes = TRUE, log = TRUE, prior.count = 2)
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

if (method == "edger") {
    dge <- estimateDisp(dge, model)
    fit <- glmQLFit(dge, model)
    test_result <- glmQLFTest(fit, coef = coefficient)
    table_result <- topTags(test_result, n = Inf, sort.by = "none")$table
    result <- data.frame(
        gene_id = rownames(table_result),
        baseMean = rowMeans(normalized_matrix[rownames(table_result), , drop = FALSE]),
        log2FoldChange = table_result$logFC,
        stat = sign(table_result$logFC) * sqrt(table_result$F),
        pvalue = table_result$PValue,
        padj = table_result$FDR,
        check.names = FALSE,
        row.names = NULL
    )
    pdf(file.path(outdir, "method-diagnostics.pdf"), width = 7, height = 6)
    plotBCV(dge)
    dev.off()
} else {
    pdf(file.path(outdir, "method-diagnostics.pdf"), width = 7, height = 6)
    voom_data <- voom(dge, model, plot = TRUE)
    dev.off()
    fit <- eBayes(lmFit(voom_data, model))
    table_result <- topTable(fit, coef = coefficient, number = Inf, sort.by = "none")
    result <- data.frame(
        gene_id = rownames(table_result),
        baseMean = rowMeans(normalized_matrix[rownames(table_result), , drop = FALSE]),
        log2FoldChange = table_result$logFC,
        stat = table_result$t,
        pvalue = table_result$P.Value,
        padj = table_result$adj.P.Val,
        check.names = FALSE,
        row.names = NULL
    )
}

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
    file.path(outdir, paste0(method, "-results.tsv")),
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
write.table(
    data.frame(gene_id = rownames(normalized_matrix), normalized_matrix, check.names = FALSE),
    file.path(outdir, "normalized-counts.tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)

sample_qc <- data.frame(
    sample = colnames(count_matrix),
    library_size = colSums(count_matrix),
    normalization_factor = dge$samples$norm.factors,
    group = group,
    metadata[, covariates, drop = FALSE],
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

colors <- as.integer(group)
pdf(file.path(outdir, "mds.pdf"), width = 7, height = 6)
plotMDS(dge, labels = colnames(dge), col = colors, main = paste(test, "vs", reference))
legend("topright", legend = levels(group), col = seq_along(levels(group)), pch = 1)
dev.off()

pdf(file.path(outdir, "ma.pdf"), width = 7, height = 6)
plot(
    log10(result$baseMean + 1),
    result$log2FoldChange,
    pch = 16,
    cex = 0.6,
    xlab = "log10 mean normalized CPM + 1",
    ylab = "log2 fold-change",
    main = paste(test, "vs", reference)
)
abline(h = 0, col = "grey50")
dev.off()

pdf(file.path(outdir, "volcano.pdf"), width = 7, height = 6)
plot(
    result$log2FoldChange,
    -log10(pmax(result$pvalue, .Machine$double.xmin)),
    col = ifelse(result$gene_id %in% significant$gene_id, "firebrick", "grey50"),
    pch = 16,
    cex = 0.6,
    xlab = "log2 fold-change",
    ylab = "-log10 p-value",
    main = paste(test, "vs", reference)
)
dev.off()

ranked <- order(result$padj, result$pvalue, na.last = NA)
ranked <- head(ranked, min(top_genes, length(ranked)))
pdf(file.path(outdir, "top-genes-heatmap.pdf"), width = 9, height = 9)
if (length(ranked) >= 2) {
    heatmap(
        log_cpm[result$gene_id[ranked], , drop = FALSE],
        scale = "row",
        ColSideColors = c("steelblue", "firebrick")[colors],
        margins = c(8, 8)
    )
} else {
    plot.new()
    text(0.5, 0.5, "Too few tested genes for a heatmap")
}
dev.off()

summary <- data.frame(
    method = if (method == "limma") "limma-voom" else "edgeR quasi-likelihood",
    contrast = paste(test, "vs", reference),
    samples = ncol(count_matrix),
    genes_tested = nrow(result),
    significant_genes = nrow(significant),
    padj_threshold = padj_threshold,
    abs_log2fc_threshold = lfc_threshold
)
write.table(
    summary,
    file.path(outdir, "analysis-summary.tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)
writeLines(capture.output(sessionInfo()), file.path(outdir, "session-info.txt"))
