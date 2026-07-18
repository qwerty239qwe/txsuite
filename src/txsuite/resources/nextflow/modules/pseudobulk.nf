process PSEUDOBULK {
    tag 'pseudobulk'
    label 'single_cell'
    container params.single_cell_image
    publishDir "${params.outdir}/pseudobulk", mode: 'copy'

    input:
    path h5ad

    output:
    path 'pseudobulk-counts.tsv', emit: counts
    path 'pseudobulk-metadata.tsv', emit: metadata

    script:
    """
    python /opt/txsuite/single_cell.py pseudobulk \
        "${h5ad}" . \
        --sample-column ${params.sample_column} \
        --design ${params.design}
    """

    stub:
    """
    printf 'gene_id\tsample_A\tsample_B\nGENE1\t10\t20\n' > pseudobulk-counts.tsv
    printf 'sample\t${params.design}\nsample_A\t${params.reference}\nsample_B\t${params.test}\n' > pseudobulk-metadata.tsv
    """
}
