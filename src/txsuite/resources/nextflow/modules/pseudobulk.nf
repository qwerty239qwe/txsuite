process PSEUDOBULK {
    tag "${meta.id}"
    label 'single_cell'

    input:
    tuple val(meta), path(h5ad)

    output:
    tuple val(meta), path("pseudobulk/${meta.id}"), emit: data

    script:
    def outputDir = "pseudobulk/${meta.id}"
    def groupArgs = meta.group_column ? "--group-column ${meta.group_column} --group-value ${meta.group_value}" : ''
    def covariateArgs = meta.covariates ? meta.covariates.split(',').collect { "--covariate ${it}" }.join(' ') : ''
    """
    mkdir -p "${outputDir}"
    python /opt/txsuite/single_cell.py pseudobulk \
        "${h5ad}" "${outputDir}" \
        --sample-column ${meta.sample_column} \
        --design ${meta.design} \
        --reference ${meta.reference} \
        --test ${meta.test} \
        ${groupArgs} ${covariateArgs}
    """

    stub:
    """
    mkdir -p "pseudobulk/${meta.id}"
    printf 'gene_id\tsample_A\tsample_B\nGENE1\t10\t20\n' > "pseudobulk/${meta.id}/pseudobulk-counts.tsv"
    printf 'sample\t${meta.design}\nsample_A\t${meta.reference}\nsample_B\t${meta.test}\n' > "pseudobulk/${meta.id}/pseudobulk-metadata.tsv"
    """
}
