suppressPackageStartupMessages(library(clusterProfiler))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 9) {
    stop(paste(
        "usage: enrichment.R MODE DE_RESULTS GENESETS_GMT OUTDIR",
        "PADJ LFC MIN_SIZE MAX_SIZE ADJUST_METHOD"
    ))
}

mode <- args[[1]]
de_path <- args[[2]]
genesets_path <- args[[3]]
outdir <- args[[4]]
padj_threshold <- as.numeric(args[[5]])
lfc_threshold <- as.numeric(args[[6]])
min_size <- as.integer(args[[7]])
max_size <- as.integer(args[[8]])
adjust_method <- args[[9]]

if (!(mode %in% c("ora", "gsea"))) {
    stop("mode must be 'ora' or 'gsea'")
}
if (!is.finite(padj_threshold) || padj_threshold <= 0 || padj_threshold > 1) {
    stop("adjusted p-value threshold must be in (0, 1]")
}
if (!is.finite(lfc_threshold) || lfc_threshold < 0) {
    stop("absolute log2 fold-change threshold must be non-negative")
}
if (is.na(min_size) || is.na(max_size) || min_size < 1 || max_size < min_size) {
    stop("gene-set sizes must satisfy 1 <= min-size <= max-size")
}
if (!(adjust_method %in% p.adjust.methods)) {
    stop(paste("unknown p-value adjustment method:", adjust_method))
}

de <- read.delim(de_path, check.names = FALSE, stringsAsFactors = FALSE)
required <- c("gene_id", "log2FoldChange", "padj")
missing_columns <- setdiff(required, colnames(de))
if (length(missing_columns)) {
    stop(paste("DE results are missing columns:", paste(missing_columns, collapse = ", ")))
}
de$gene_id <- as.character(de$gene_id)
if (anyNA(de$gene_id) || any(!nzchar(de$gene_id)) || anyDuplicated(de$gene_id)) {
    stop("DE results must contain unique, non-empty gene IDs")
}

gmt_lines <- readLines(genesets_path, warn = FALSE)
gmt_parts <- strsplit(gmt_lines[nzchar(gmt_lines)], "\t", fixed = TRUE)
if (!length(gmt_parts) || any(lengths(gmt_parts) < 3)) {
    stop("GMT rows require a term, description, and at least one gene")
}
term_to_gene <- unique(do.call(rbind, lapply(gmt_parts, function(parts) {
    data.frame(term = parts[[1]], gene = parts[-c(1, 2)], stringsAsFactors = FALSE)
})))
term_to_name <- unique(do.call(rbind, lapply(gmt_parts, function(parts) {
    data.frame(term = parts[[1]], name = parts[[2]], stringsAsFactors = FALSE)
})))
overlap <- intersect(de$gene_id, term_to_gene$gene)
if (!length(overlap)) {
    stop("no DE result gene IDs occur in the GMT file")
}

empty_result <- function(mode) {
    if (mode == "ora") {
        return(data.frame(
            ID = character(), Description = character(), GeneRatio = character(),
            BgRatio = character(), pvalue = numeric(), p.adjust = numeric(),
            qvalue = numeric(), geneID = character(), Count = integer()
        ))
    }
    data.frame(
        ID = character(), Description = character(), setSize = integer(),
        enrichmentScore = numeric(), NES = numeric(), pvalue = numeric(),
        p.adjust = numeric(), qvalue = numeric(), rank = integer(),
        leading_edge = character(), core_enrichment = character()
    )
}

if (mode == "ora") {
    selected <- de$gene_id[
        !is.na(de$padj) &
            de$padj <= padj_threshold &
            abs(de$log2FoldChange) >= lfc_threshold
    ]
    selected <- intersect(selected, term_to_gene$gene)
    enrichment <- if (length(selected)) {
        enricher(
            gene = selected,
            universe = de$gene_id,
            pvalueCutoff = 1,
            qvalueCutoff = 1,
            pAdjustMethod = adjust_method,
            minGSSize = min_size,
            maxGSSize = max_size,
            TERM2GENE = term_to_gene,
            TERM2NAME = term_to_name
        )
    } else {
        NULL
    }
    result <- if (is.null(enrichment)) empty_result(mode) else as.data.frame(enrichment)
    output_name <- "ora-results.tsv"
} else {
    rank_column <- if ("stat" %in% colnames(de)) "stat" else "log2FoldChange"
    ranked <- de[[rank_column]]
    names(ranked) <- de$gene_id
    ranked <- sort(ranked[is.finite(ranked) & names(ranked) %in% term_to_gene$gene], decreasing = TRUE)
    if (length(ranked) < min_size) {
        stop("too few finite ranked genes overlap the GMT file")
    }
    enrichment <- GSEA(
        geneList = ranked,
        pvalueCutoff = 1,
        pAdjustMethod = adjust_method,
        minGSSize = min_size,
        maxGSSize = max_size,
        TERM2GENE = term_to_gene,
        TERM2NAME = term_to_name,
        verbose = FALSE,
        seed = TRUE,
        by = "fgsea"
    )
    result <- if (is.null(enrichment)) empty_result(mode) else as.data.frame(enrichment)
    output_name <- "gsea-results.tsv"
}

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
write.table(
    result,
    file.path(outdir, output_name),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)

pdf(file.path(outdir, "enrichment-top.pdf"), width = 10, height = 7)
if (nrow(result) && any(is.finite(result$p.adjust))) {
    ranked_rows <- head(order(result$p.adjust, na.last = NA), 20L)
    labels <- if ("Description" %in% colnames(result)) result$Description[ranked_rows] else result$ID[ranked_rows]
    values <- if (mode == "gsea") result$NES[ranked_rows] else -log10(pmax(result$p.adjust[ranked_rows], .Machine$double.xmin))
    colors <- if (mode == "gsea") ifelse(values >= 0, "firebrick", "steelblue") else "steelblue"
    par(mar = c(5, 16, 3, 2))
    barplot(
        rev(values),
        names.arg = rev(substr(labels, 1, 60)),
        horiz = TRUE,
        las = 1,
        col = rev(colors),
        xlab = if (mode == "gsea") "normalized enrichment score" else "-log10 adjusted p-value",
        main = toupper(mode)
    )
} else {
    plot.new()
    text(0.5, 0.5, "No gene sets passed size and overlap filters")
}
dev.off()

significant_sets <- if (nrow(result) && "p.adjust" %in% colnames(result)) {
    sum(result$p.adjust <= padj_threshold, na.rm = TRUE)
} else {
    0L
}
write.table(
    data.frame(
        metric = c("mode", "overlapping_genes", "tested_gene_sets", "significant_gene_sets"),
        value = c(mode, length(overlap), nrow(result), significant_sets)
    ),
    file.path(outdir, "enrichment-summary.tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)
writeLines(capture.output(sessionInfo()), file.path(outdir, "session-info.txt"))
