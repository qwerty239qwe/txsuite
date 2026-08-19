nextflow.enable.dsl = 2

include { SIMPLEAF_INDEX } from './modules/salmon'
include { SIMPLEAF_QUANT } from './modules/salmon'
include { ALEVIN_TO_H5AD } from './modules/salmon'

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
            if (!row.fastq_1 || !row.fastq_2) {
                error "Droplet chemistry needs both fastq_1 and fastq_2 for ${sample}"
            }
            def reads = [
                file(row.fastq_1, checkIfExists: true),
                file(row.fastq_2, checkIfExists: true),
            ]
            return tuple(sample, reads)
        }
}

workflow {
    if (!params.samplesheet) {
        error 'Missing required parameter: --samplesheet'
    }
    if (!params.simpleaf_index && !(params.fasta && params.gtf)) {
        error 'Provide --simpleaf_index, or both --fasta and --gtf to build one'
    }
    if (params.simpleaf_index && (params.fasta || params.gtf)) {
        error 'Pass either --simpleaf_index, or --fasta and --gtf, not both'
    }

    reads_ch = readSamplesheet(params.samplesheet)

    if (params.simpleaf_index) {
        index_ch = channel.fromPath(params.simpleaf_index, checkIfExists: true, type: 'dir')
    } else {
        fasta_ch = channel.fromPath(params.fasta, checkIfExists: true)
        gtf_ch = channel.fromPath(params.gtf, checkIfExists: true)
        SIMPLEAF_INDEX(fasta_ch, gtf_ch)
        index_ch = SIMPLEAF_INDEX.out.index
    }

    SIMPLEAF_QUANT(index_ch.first(), reads_ch)
    ALEVIN_TO_H5AD(SIMPLEAF_QUANT.out.quant.map { it[1] }.collect())
}
