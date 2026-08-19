process STAR_ALIGN {
    tag "$sample"
    label 'star'
    container params.star_image
    publishDir "${params.outdir}/alignments", mode: 'copy', pattern: '*/*.bam'
    publishDir "${params.outdir}/logs", mode: 'copy', pattern: '*/*.out'

    input:
    path index
    tuple val(sample), path(reads)

    output:
    tuple val(sample), path("${sample}/Aligned.sortedByCoord.out.bam"), emit: bam
    tuple val(sample), path("${sample}"), emit: quant
    path "${sample}/Log.final.out", emit: log

    script:
    def read_list = reads instanceof List ? reads : [reads]
    def gzipped = read_list[0].name.endsWith('.gz') ? '--readFilesCommand zcat' : ''
    def two_pass = params.star_two_pass ? '--twopassMode Basic' : ''
    """
    mkdir -p ${sample}
    STAR --genomeDir ${index} \
        --readFilesIn ${read_list.join(' ')} \
        --runThreadN ${task.cpus} \
        --outFileNamePrefix ${sample}/ \
        --outSAMtype BAM SortedByCoordinate \
        --quantMode GeneCounts ${gzipped} ${two_pass}
    test -s ${sample}/Aligned.sortedByCoord.out.bam
    test -s ${sample}/ReadsPerGene.out.tab
    """

    stub:
    """
    mkdir -p ${sample}
    printf 'stub-bam\\n' > ${sample}/Aligned.sortedByCoord.out.bam
    printf 'stub log\\n' > ${sample}/Log.final.out
    cat > ${sample}/ReadsPerGene.out.tab <<STUB
N_unmapped	100	100	100
N_multimapping	50	50	50
N_noFeature	900	900	900
N_ambiguous	10	10	10
gene1	1000	20	980
gene2	500	10	490
STUB
    """
}

process SAMTOOLS_INDEX {
    tag "$sample"
    label 'star'
    container params.star_image
    // Every sample's BAM is named Aligned.sortedByCoord.out.bam, so the index
    // must publish under the sample directory or samples overwrite each other.
    publishDir "${params.outdir}/alignments/${sample}", mode: 'copy'

    input:
    tuple val(sample), path(bam)

    output:
    tuple val(sample), path("${bam}.bai"), emit: bai

    script:
    """
    samtools index -@ ${task.cpus} ${bam}
    """

    stub:
    """
    printf 'stub-bai\\n' > ${bam}.bai
    """
}

process MERGE_GENE_COUNTS {
    tag 'counts'
    label 'star'
    container params.star_image
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path quant_dirs, stageAs: 'quant/*'

    output:
    path 'counts/gene_counts.tsv', emit: counts
    path 'counts/strandedness.tsv', emit: strandedness
    path 'counts', emit: counts_dir

    script:
    // Strandedness is inferred from the forward/reverse split unless pinned, and
    // the evidence is published beside the matrix either way.
    def pinned = params.star_strandedness == 'auto' ? '' : "--strandedness ${params.star_strandedness}"
    """
    python /opt/txsuite/merge_star_counts.py quant/* --outdir counts ${pinned}
    """

    stub:
    """
    mkdir -p counts
    printf 'gene_id\\tsample1\\ngene1\\t980\\ngene2\\t490\\n' > counts/gene_counts.tsv
    printf 'sample\\tstrandedness\\tsource\\tunstranded\\tforward\\treverse\\tforward_fraction\\tunmapped\\tno_feature\\nsample1\\treverse\\tinferred\\t1500\\t30\\t1470\\t0.0200\\t100\\t900\\n' > counts/strandedness.tsv
    """
}
