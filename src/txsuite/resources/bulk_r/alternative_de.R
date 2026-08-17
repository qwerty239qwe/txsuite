suppressPackageStartupMessages({
    library(edgeR)
    library(limma)
})

args <- commandArgs(trailingOnly = TRUE)
if (!(length(args) %in% c(11, 13, 14))) {
    stop(paste(
        "usage: alternative_de.R METHOD COUNTS METADATA DESIGN REFERENCE TEST",
        "OUTDIR PADJ LFC TOP_GENES COVARIATES [FORMULA COEFFICIENT [CONTRASTS]]"
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
formula_text <- if (length(args) >= 13) args[[12]] else ""
coefficient_name <- if (length(args) >= 13) args[[13]] else ""
contrast_mode <- if (length(args) == 14 && nzchar(args[[14]])) args[[14]] else "single"
advanced <- nzchar(formula_text) || nzchar(coefficient_name)

if (!(method %in% c("edger", "limma"))) {
    stop("method must be 'edger' or 'limma'")
}
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
model_formula <- if (advanced) as.formula(formula_text) else NULL
model_columns <- if (advanced) all.vars(model_formula) else c(covariates, design)
missing_columns <- setdiff(model_columns, colnames(metadata))
if (length(missing_columns)) {
    stop(paste("metadata is missing model columns:", paste(missing_columns, collapse = ", ")))
}
if (anyNA(metadata[, model_columns, drop = FALSE])) {
    stop("model columns cannot contain missing values")
}

for (column in model_columns) {
    if (is.character(metadata[[column]])) {
        metadata[[column]] <- factor(metadata[[column]])
    }
}
# The primary contrast subsets the data to two levels below. Expanded contrasts
# are computed at the end of the script and need the unsubsetted inputs.
full_metadata <- metadata
full_count_matrix <- count_matrix

if (advanced) {
    model <- model.matrix(model_formula, metadata)
    if (!(coefficient_name %in% colnames(model))) {
        stop(paste(
            "coefficient not found; available coefficients:",
            paste(colnames(model), collapse = ", ")
        ))
    }
    coefficient <- match(coefficient_name, colnames(model))
    group <- factor(metadata[[model_columns[[1]]]])
    contrast_label <- coefficient_name
} else {
    selected <- metadata[[design]] %in% c(reference, test)
    metadata <- droplevels(metadata[selected, , drop = FALSE])
    count_matrix <- count_matrix[, rownames(metadata), drop = FALSE]
    group <- factor(metadata[[design]], levels = c(reference, test))
    if (anyNA(group) || any(table(group) < 2)) {
        stop("reference and test must each contain at least two samples")
    }
    metadata$.txsuite_test <- as.integer(group == test)
    model <- model.matrix(reformulate(c(covariates, ".txsuite_test")), metadata)
    coefficient <- match(".txsuite_test", colnames(model))
    contrast_label <- paste(test, "vs", reference)
}
if (qr(model)$rank < ncol(model)) {
    stop("the design matrix is not full rank; check covariates and groups")
}

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
    metadata[, model_columns, drop = FALSE],
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
palette <- grDevices::rainbow(max(1, nlevels(group)))
pdf(file.path(outdir, "mds.pdf"), width = 7, height = 6)
plotMDS(dge, labels = colnames(dge), col = colors, main = contrast_label)
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
    main = contrast_label
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
    main = contrast_label
)
dev.off()

ranked <- order(result$padj, result$pvalue, na.last = NA)
ranked <- head(ranked, min(top_genes, length(ranked)))
pdf(file.path(outdir, "top-genes-heatmap.pdf"), width = 9, height = 9)
if (length(ranked) >= 2) {
    heatmap(
        log_cpm[result$gene_id[ranked], , drop = FALSE],
        scale = "row",
        ColSideColors = palette[colors],
        margins = c(8, 8)
    )
} else {
    plot.new()
    text(0.5, 0.5, "Too few tested genes for a heatmap")
}
dev.off()

summary <- data.frame(
    method = if (method == "limma") "limma-voom" else "edgeR quasi-likelihood",
    contrast = contrast_label,
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

# Expanded contrasts.
#
# Unlike the DESeq2 script, edgeR and limma cannot read extra contrasts off the
# primary fit: this script models each comparison as a two-level subset, so the
# model matrix only knows about the reference and test levels. Each additional
# contrast therefore repeats the subset, filtering, and fit on the full inputs.
# The upside is that the primary contrast's numbers are unchanged by this
# feature, and every contrast is computed exactly the way a single-contrast run
# would compute it.
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
        for (pair in combn(levels_present, 2, simplify = FALSE)) {
            pairs[[length(pairs) + 1]] <- c(pair[[2]], pair[[1]])
        }
    }
    Filter(function(pair) !(pair[[1]] == test && pair[[2]] == reference), pairs)
}
pairwise_result <- function(ref_level, test_level) {
    selected <- full_metadata[[design]] %in% c(ref_level, test_level)
    pair_metadata <- droplevels(full_metadata[selected, , drop = FALSE])
    pair_counts <- full_count_matrix[, rownames(pair_metadata), drop = FALSE]
    pair_group <- factor(pair_metadata[[design]], levels = c(ref_level, test_level))
    if (anyNA(pair_group) || any(table(pair_group) < 2)) {
        stop(paste(
            "each level in a contrast needs at least two samples;",
            test_level, "vs", ref_level, "does not"
        ))
    }
    pair_metadata$.txsuite_test <- as.integer(pair_group == test_level)
    pair_model <- model.matrix(reformulate(c(covariates, ".txsuite_test")), pair_metadata)
    if (qr(pair_model)$rank < ncol(pair_model)) {
        stop(paste(
            "the design matrix is not full rank for", test_level, "vs", ref_level
        ))
    }
    pair_coefficient <- match(".txsuite_test", colnames(pair_model))
    pair_dge <- DGEList(counts = pair_counts)
    keep_pair <- filterByExpr(pair_dge, design = pair_model)
    if (!any(keep_pair)) {
        stop(paste("no genes pass expression filtering for", test_level, "vs", ref_level))
    }
    pair_dge <- calcNormFactors(pair_dge[keep_pair, , keep.lib.sizes = FALSE])
    pair_normalized <- cpm(pair_dge, normalized.lib.sizes = TRUE, log = FALSE)
    if (method == "edger") {
        pair_dge <- estimateDisp(pair_dge, pair_model)
        pair_fit <- glmQLFit(pair_dge, pair_model)
        pair_table <- topTags(
            glmQLFTest(pair_fit, coef = pair_coefficient), n = Inf, sort.by = "none"
        )$table
        data.frame(
            gene_id = rownames(pair_table),
            baseMean = rowMeans(pair_normalized[rownames(pair_table), , drop = FALSE]),
            log2FoldChange = pair_table$logFC,
            stat = sign(pair_table$logFC) * sqrt(pair_table$F),
            pvalue = pair_table$PValue,
            padj = pair_table$FDR,
            check.names = FALSE,
            row.names = NULL
        )
    } else {
        pair_voom <- voom(pair_dge, pair_model, plot = FALSE)
        pair_fit <- eBayes(lmFit(pair_voom, pair_model))
        pair_table <- topTable(
            pair_fit, coef = pair_coefficient, number = Inf, sort.by = "none"
        )
        data.frame(
            gene_id = rownames(pair_table),
            baseMean = rowMeans(pair_normalized[rownames(pair_table), , drop = FALSE]),
            log2FoldChange = pair_table$logFC,
            stat = pair_table$t,
            pvalue = pair_table$P.Value,
            padj = pair_table$adj.P.Val,
            check.names = FALSE,
            row.names = NULL
        )
    }
}

contrast_index <- data.frame(
    contrast_id = if (advanced) coefficient_name else paste0(test, "_vs_", reference),
    reference = if (advanced) "" else reference,
    test = if (advanced) "" else test,
    genes_tested = nrow(result),
    significant_genes = nrow(significant),
    path = paste0(method, "-results.tsv"),
    stringsAsFactors = FALSE
)
if (!advanced) {
    design_levels <- levels(factor(full_metadata[[design]]))
    extra_contrasts <- expanded_contrasts(contrast_mode, design_levels, reference, test)
    if (length(extra_contrasts) + 1L > 50L) {
        stop(paste0(
            "contrasts='", contrast_mode, "' expands column '", design, "' with ",
            length(design_levels), " levels into ", length(extra_contrasts) + 1L,
            " comparisons; select levels or use 'vs-reference'"
        ))
    }
    if (length(extra_contrasts)) {
        dir.create(file.path(outdir, "contrasts"), recursive = TRUE, showWarnings = FALSE)
    }
    for (pair in extra_contrasts) {
        pair_table <- pairwise_result(pair[[2]], pair[[1]])
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
writeLines(capture.output(sessionInfo()), file.path(outdir, "session-info.txt"))
