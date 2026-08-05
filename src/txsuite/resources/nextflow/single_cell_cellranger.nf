nextflow.enable.dsl = 2

include { CELLRANGER_MKREF } from './modules/cellranger'
include { CELLRANGER_COUNT } from './modules/cellranger'

workflow {
    if (!params.genome_name) {
        error 'Missing required parameter: --genome_name'
    }
    if (!(params.genome_name ==~ /^[A-Za-z0-9][A-Za-z0-9._-]*$/)) {
        error '--genome_name must contain only letters, numbers, ., _ and -'
    }
    if (!params.fastqs) {
        error 'Missing required parameter: --fastqs'
    }
    if (!params.sample) {
        error 'Missing required parameter: --sample'
    }
    if (!params.cellranger_reference && !(params.fasta && params.gtf)) {
        error 'Provide --cellranger_reference, or both --fasta and --gtf to build one'
    }

    fastqs_ch = channel.fromPath(params.fastqs, checkIfExists: true, type: 'dir')

    if (params.cellranger_reference) {
        reference_ch = channel.fromPath(params.cellranger_reference, checkIfExists: true, type: 'dir')
    } else {
        fasta_ch = channel.fromPath(params.fasta, checkIfExists: true)
        gtf_ch = channel.fromPath(params.gtf, checkIfExists: true)
        CELLRANGER_MKREF(fasta_ch, gtf_ch)
        reference_ch = CELLRANGER_MKREF.out.reference
    }

    CELLRANGER_COUNT(reference_ch, fastqs_ch)
}
