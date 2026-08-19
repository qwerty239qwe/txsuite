nextflow.enable.dsl = 2

include { STAR_ALIGN } from './modules/star_align'
include { SAMTOOLS_INDEX } from './modules/star_align'
include { MERGE_GENE_COUNTS } from './modules/star_align'

def readSamplesheet(path) {
    return channel
        .fromPath(path, checkIfExists: true)
        .splitCsv(header: true, strip: true)
        .map { row ->
            def sample = row.sample ?: row.id
            if (!sample) {
                error 'Samplesheet rows require a non-empty sample column'
            }
            if (!(sample ==~ /^[A-Za-z0-9][A-Za-z0-9._-]*$/)) {
                error "Sample name may contain only letters, numbers, ., _ and -: ${sample}"
            }
            if (!row.fastq_1) {
                error "Samplesheet row ${sample} is missing fastq_1"
            }
            def reads = [file(row.fastq_1, checkIfExists: true)]
            if (row.fastq_2) {
                reads << file(row.fastq_2, checkIfExists: true)
            }
            return tuple(sample, reads)
        }
}

workflow {
    if (!params.samplesheet) {
        error 'Missing required parameter: --samplesheet'
    }
    if (!params.star_index) {
        error 'Missing required parameter: --star_index'
    }

    reads_ch = readSamplesheet(params.samplesheet)
    index_ch = channel.fromPath(params.star_index, checkIfExists: true, type: 'dir')

    STAR_ALIGN(index_ch.first(), reads_ch)
    SAMTOOLS_INDEX(STAR_ALIGN.out.bam)
    MERGE_GENE_COUNTS(STAR_ALIGN.out.quant.map { it[1] }.collect())
}
