process BULK_DE {
    tag "${meta.id}:${meta.method}"
    label 'bulk_r'

    input:
    tuple val(meta), path(pseudobulk_dir)

    output:
    tuple val(meta), path("de/${meta.id}"), emit: results

    script:
    def outputDir = "de/${meta.id}"
    def counts = "${pseudobulk_dir}/pseudobulk-counts.tsv"
    def metadata = "${pseudobulk_dir}/pseudobulk-metadata.tsv"
    def command = meta.method == 'deseq2' \
        ? "Rscript /opt/txsuite/deseq2.R \"${counts}\" \"${metadata}\" ${meta.design} ${meta.reference} ${meta.test}" \
        : "Rscript /opt/txsuite/alternative_de.R ${meta.method} \"${counts}\" \"${metadata}\" ${meta.design} ${meta.reference} ${meta.test}"
    """
    mkdir -p "${outputDir}"
    ${command} "${outputDir}" ${meta.padj} ${meta.lfc} ${meta.top_genes} "${meta.covariates}"
    """

    stub:
    """
    mkdir -p "de/${meta.id}"
    printf 'gene_id\tbaseMean\tlog2FoldChange\tpadj\nGENE1\t15\t1\t0.05\n' > "de/${meta.id}/${meta.method}-results.tsv"
    printf 'gene_id\tbaseMean\tlog2FoldChange\tpadj\nGENE1\t15\t1\t0.05\n' > "de/${meta.id}/significant-genes.tsv"
    printf 'metric\tvalue\nstatus\tstub\n' > "de/${meta.id}/analysis-summary.tsv"
    touch "de/${meta.id}/pca.pdf" "de/${meta.id}/session-info.txt"
    """
}
