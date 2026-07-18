process DESEQ2 {
    tag "${params.test}_vs_${params.reference}"
    label 'bulk_r'
    container params.bulk_image
    publishDir "${params.outdir}/deseq2", mode: 'copy'

    input:
    path counts
    path metadata

    output:
    path '*.tsv', emit: tables
    path '*.pdf', emit: figures
    path 'session-info.txt', emit: session_info

    script:
    """
    Rscript /opt/txsuite/deseq2.R \
        "${counts}" "${metadata}" \
        ${params.design} ${params.reference} ${params.test} . \
        ${params.padj} ${params.lfc} ${params.top_genes} ''
    """

    stub:
    """
    printf 'gene_id\tbaseMean\tlog2FoldChange\tpadj\nGENE1\t15\t1\t0.05\n' > deseq2-results.tsv
    printf 'metric\tvalue\nstatus\tstub\n' > analysis-summary.tsv
    touch pca.pdf session-info.txt
    """
}
