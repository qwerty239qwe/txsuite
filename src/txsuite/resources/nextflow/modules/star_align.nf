process STAR_ALIGN {
    tag "$sample"
    label 'star'
    container params.star_image
    // STAR is told to prefix its outputs with the sample, so every file is
    // unique and lives at the top level of the task directory. Writing into a
    // per-sample subdirectory instead made the whole directory the single
    // collected output, and publishDir then matched neither the BAM nor the log.
    publishDir "${params.outdir}/alignments", mode: 'copy', pattern: '*.bam'
    publishDir "${params.outdir}/logs", mode: 'copy', pattern: '*.Log.final.out'

    input:
    path index
    tuple val(sample), path(reads)

    output:
    tuple val(sample), path("${sample}.Aligned.sortedByCoord.out.bam"), emit: bam
    path "${sample}.ReadsPerGene.out.tab", emit: counts
    path "${sample}.Log.final.out", emit: log

    script:
    def read_list = reads instanceof List ? reads : [reads]
    def gzipped = read_list[0].name.endsWith('.gz') ? '--readFilesCommand zcat' : ''
    def two_pass = params.star_two_pass ? '--twopassMode Basic' : ''
    """
    STAR --genomeDir ${index} \
        --readFilesIn ${read_list.join(' ')} \
        --runThreadN ${task.cpus} \
        --outFileNamePrefix ${sample}. \
        --outSAMtype BAM SortedByCoordinate \
        --quantMode GeneCounts ${gzipped} ${two_pass}
    test -s ${sample}.Aligned.sortedByCoord.out.bam
    test -s ${sample}.ReadsPerGene.out.tab
    """

    stub:
    """
    printf 'stub-bam\\n' > ${sample}.Aligned.sortedByCoord.out.bam
    printf 'stub log\\n' > ${sample}.Log.final.out
    cat > ${sample}.ReadsPerGene.out.tab <<STUB
N_unmapped\t100\t100\t100
N_multimapping\t50\t50\t50
N_noFeature\t900\t900\t900
N_ambiguous\t10\t10\t10
gene1\t1000\t20\t980
gene2\t500\t10\t490
STUB
    """
}

process SAMTOOLS_INDEX {
    tag "$sample"
    label 'star'
    container params.star_image
    // The BAM name already carries the sample, so indexes cannot collide and
    // the index publishes flat alongside it.
    publishDir "${params.outdir}/alignments", mode: 'copy'

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
    path counts_files, stageAs: 'counts_in/*'

    output:
    path 'counts/gene_counts.tsv', emit: counts
    path 'counts/strandedness.tsv', emit: strandedness
    path 'counts', emit: counts_dir

    script:
    // Strandedness is inferred from the forward/reverse split unless pinned, and
    // the evidence is published beside the matrix either way.
    def pinned = params.star_strandedness == 'auto' ? '' : "--strandedness ${params.star_strandedness}"
    """
    python /opt/txsuite/merge_star_counts.py counts_in/* --outdir counts ${pinned}
    """

    stub:
    """
    mkdir -p counts
    printf 'gene_id\\tsample1\\ngene1\\t980\\ngene2\\t490\\n' > counts/gene_counts.tsv
    printf 'sample\\tstrandedness\\tsource\\tunstranded\\tforward\\treverse\\tforward_fraction\\tunmapped\\tno_feature\\nsample1\\treverse\\tinferred\\t1500\\t30\\t1470\\t0.0200\\t100\\t900\\n' > counts/strandedness.tsv
    """
}
