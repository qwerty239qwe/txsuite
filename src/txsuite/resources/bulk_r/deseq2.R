suppressPackageStartupMessages(library(DESeq2))

args <- commandArgs(trailingOnly = TRUE)
if (!(length(args) %in% c(6, 10, 12, 13))) {
    stop(paste(
        "usage: deseq2.R COUNTS METADATA DESIGN REFERENCE TEST OUTDIR",
        "[PADJ LFC TOP_GENES COVARIATES [FORMULA COEFFICIENT [CONTRASTS]]]"
    ))
}

counts_path <- args[[1]]
metadata_path <- args[[2]]
design <- args[[3]]
reference <- args[[4]]
test <- args[[5]]
outdir <- args[[6]]
padj_threshold <- if (length(args) >= 10) as.numeric(args[[7]]) else 0.05
lfc_threshold <- if (length(args) >= 10) as.numeric(args[[8]]) else 1
top_genes <- if (length(args) >= 10) as.integer(args[[9]]) else 50L
covariates <- if (length(args) >= 10 && nzchar(args[[10]])) {
    strsplit(args[[10]], ",", fixed = TRUE)[[1]]
} else {
    character()
}
formula_text <- if (length(args) >= 12) args[[11]] else ""
coefficient_name <- if (length(args) >= 12) args[[12]] else ""
contrast_mode <- if (length(args) == 13 && nzchar(args[[13]])) args[[13]] else "single"
advanced <- nzchar(formula_text) || nzchar(coefficient_name)

if (advanced && (!nzchar(formula_text) || !nzchar(coefficient_name))) {
    stop("formula and coefficient must be provided together")
}
if (advanced && (nzchar(design) || nzchar(reference) || nzchar(test) || length(covariates))) {
    stop("formula mode cannot be combined with design, reference, test, or covariates")
}
if (advanced && (!grepl("^~[A-Za-z0-9_.+*: ]+$", formula_text) ||
    !grepl("[A-Za-z]", formula_text))) {
    stop("formula may contain only metadata names and +, *, or : operators")
}
if (advanced && !grepl("^[A-Za-z][A-Za-z0-9_.:]*$", coefficient_name)) {
    stop("coefficient must be a model coefficient name")
}
if (!(contrast_mode %in% c("single", "vs-reference", "all-pairs"))) {
    stop("contrasts must be 'single', 'vs-reference', or 'all-pairs'")
}
if (advanced && contrast_mode != "single") {
    stop("formula mode has no design levels to expand; contrasts must be 'single'")
}
if (!advanced && !grepl("^[A-Za-z][A-Za-z0-9_.]*$", design)) {
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
model_formula <- if (advanced) as.formula(formula_text) else reformulate(c(covariates, design))
model_columns <- all.vars(model_formula)
missing_columns <- setdiff(model_columns, colnames(metadata))
if (length(missing_columns)) {
    stop(paste("metadata is missing model columns:", paste(missing_columns, collapse = ", ")))
}
for (column in model_columns) {
    if (is.character(metadata[[column]])) {
        metadata[[column]] <- factor(metadata[[column]])
    }
}
if (anyNA(metadata[, model_columns, drop = FALSE])) {
    stop("model columns cannot contain missing values")
}
if (!advanced) {
    metadata[[design]] <- factor(metadata[[design]])
    levels_present <- levels(metadata[[design]])
    if (!(reference %in% levels_present) || !(test %in% levels_present)) {
        stop("reference and test must both occur in the design column")
    }
    metadata[[design]] <- relevel(metadata[[design]], ref = reference)
}
plot_group <- factor(metadata[[model_columns[[1]]]])
contrast_label <- if (advanced) coefficient_name else paste(test, "vs", reference)

dds <- DESeqDataSetFromMatrix(
    countData = count_matrix,
    colData = metadata,
    design = model_formula
)
dds <- DESeq(dds)
if (advanced) {
    if (!(coefficient_name %in% resultsNames(dds))) {
        stop(paste(
            "coefficient not found; available coefficients:",
            paste(resultsNames(dds), collapse = ", ")
        ))
    }
    de_result <- results(dds, name = coefficient_name, alpha = padj_threshold)
} else {
    de_result <- results(
        dds,
        contrast = c(design, test, reference),
        alpha = padj_threshold
    )
}

# Every contrast is read off the single fit above. Sharing one dispersion
# estimate is what makes expanded contrasts both cheaper and mutually
# consistent, unlike running one model per comparison.
as_result_table <- function(de) {
    data.frame(gene_id = rownames(de), as.data.frame(de), check.names = FALSE)
}
select_significant <- function(table) {
    table[
        !is.na(table$padj) &
            table$padj <= padj_threshold &
            abs(table$log2FoldChange) >= lfc_threshold,
        ,
        drop = FALSE
    ]
}
expanded_contrasts <- function(mode, levels_present, reference, test) {
    if (mode == "single") {
        return(list())
    }
    pairs <- list()
    if (mode == "vs-reference") {
        for (level in setdiff(levels_present, reference)) {
            pairs[[length(pairs) + 1]] <- c(level, reference)
        }
    } else {
        combinations <- combn(levels_present, 2, simplify = FALSE)
        for (pair in combinations) {
            pairs[[length(pairs) + 1]] <- c(pair[[2]], pair[[1]])
        }
    }
    # The primary contrast is written at the output root, so drop it here.
    Filter(function(pair) !(pair[[1]] == test && pair[[2]] == reference), pairs)
}

result <- as_result_table(de_result)
significant <- select_significant(result)

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

contrast_index <- data.frame(
    contrast_id = if (advanced) coefficient_name else paste0(test, "_vs_", reference),
    reference = if (advanced) "" else reference,
    test = if (advanced) "" else test,
    genes_tested = nrow(result),
    significant_genes = nrow(significant),
    path = "deseq2-results.tsv",
    stringsAsFactors = FALSE
)
if (!advanced) {
    extra_contrasts <- expanded_contrasts(
        contrast_mode, levels(metadata[[design]]), reference, test
    )
    if (length(extra_contrasts) + 1L > 50L) {
        stop(paste0(
            "contrasts='", contrast_mode, "' expands column '", design, "' with ",
            nlevels(metadata[[design]]), " levels into ", length(extra_contrasts) + 1L,
            " comparisons; select levels or use 'vs-reference'"
        ))
    }
    if (length(extra_contrasts)) {
        dir.create(file.path(outdir, "contrasts"), recursive = TRUE, showWarnings = FALSE)
    }
    for (pair in extra_contrasts) {
        pair_table <- as_result_table(
            results(dds, contrast = c(design, pair[[1]], pair[[2]]), alpha = padj_threshold)
        )
        relative <- file.path("contrasts", paste0("DE_", pair[[1]], "_vs_", pair[[2]], ".tsv"))
        write.table(
            pair_table,
            file.path(outdir, relative),
            sep = "\t",
            quote = FALSE,
            row.names = FALSE
        )
        contrast_index <- rbind(contrast_index, data.frame(
            contrast_id = paste0(pair[[1]], "_vs_", pair[[2]]),
            reference = pair[[2]],
            test = pair[[1]],
            genes_tested = nrow(pair_table),
            significant_genes = nrow(select_significant(pair_table)),
            path = relative,
            stringsAsFactors = FALSE
        ))
    }
}
write.table(
    contrast_index,
    file.path(outdir, "contrasts.tsv"),
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
group_factor <- plot_group[match(rownames(pca_scores), rownames(metadata))]
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
    main = contrast_label
)
text(pca_scores[, 1], pca_scores[, 2], labels = rownames(pca_scores), pos = 3, cex = 0.7)
legend("topright", legend = levels(group_factor), col = seq_along(levels(group_factor)), pch = 19)
dev.off()

pdf(file.path(outdir, "sample-correlation.pdf"), width = 8, height = 8)
heatmap(sample_correlation, scale = "none", symm = TRUE, margins = c(9, 9))
dev.off()

pdf(file.path(outdir, "ma.pdf"), width = 7, height = 6)
plotMA(de_result, alpha = padj_threshold, main = contrast_label)
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
    main = contrast_label
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
