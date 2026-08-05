process CELLRANGER_MKREF {
    tag params.genome_name
    label 'cellranger'
    container params.cellranger_image
    publishDir "${params.outdir}/reference", mode: 'copy'

    input:
    path fasta
    path gtf

    output:
    path "${params.genome_name}", emit: reference

    script:
    """
    cellranger mkref \
        --genome=${params.genome_name} \
        --fasta=${fasta} \
        --genes=${gtf} \
        --nthreads=${params.cellranger_threads} \
        --memgb=${params.cellranger_memory_gb}
    """

    stub:
    """
    mkdir -p ${params.genome_name}/fasta ${params.genome_name}/genes ${params.genome_name}/star
    cat > ${params.genome_name}/reference.json <<STUB
{"genomes": ["${params.genome_name}"], "fasta_hash": "stub", "gtf_hash": "stub"}
STUB
    """
}

process CELLRANGER_COUNT {
    tag params.sample
    label 'cellranger'
    container params.cellranger_image
    publishDir "${params.outdir}/counts", mode: 'copy'

    input:
    path reference
    path fastqs

    output:
    path "${params.sample}/outs", emit: outs
    path "${params.sample}/outs/filtered_feature_bc_matrix", emit: matrix
    path "${params.sample}/outs/web_summary.html", emit: web_summary

    script:
    """
    cellranger count \
        --id=${params.sample} \
        --transcriptome=${reference} \
        --fastqs=${fastqs} \
        --sample=${params.sample} \
        --create-bam=${params.cellranger_create_bam} \
        --localcores=${params.cellranger_threads} \
        --localmem=${params.cellranger_memory_gb}
    """

    stub:
    """
    mkdir -p ${params.sample}/outs/filtered_feature_bc_matrix
    cat > barcodes.tsv <<STUB
AAACCTGAGAAACCAT-1
AAACCTGAGAAACCGC-1
STUB
    cat > features.tsv <<STUB
ENSG00000000001	GENE1	Gene Expression
ENSG00000000002	GENE2	Gene Expression
STUB
    cat > matrix.mtx <<STUB
%%MatrixMarket matrix coordinate integer general
%
2 2 2
1 1 3
2 2 5
STUB
    gzip -c barcodes.tsv > ${params.sample}/outs/filtered_feature_bc_matrix/barcodes.tsv.gz
    gzip -c features.tsv > ${params.sample}/outs/filtered_feature_bc_matrix/features.tsv.gz
    gzip -c matrix.mtx > ${params.sample}/outs/filtered_feature_bc_matrix/matrix.mtx.gz
    rm barcodes.tsv features.tsv matrix.mtx
    touch ${params.sample}/outs/web_summary.html
    """
}
