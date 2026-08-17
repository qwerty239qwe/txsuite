nextflow.enable.dsl = 2

include { SALMON_REFERENCE } from './modules/salmon'
include { SALMON_INDEX } from './modules/salmon'
include { SALMON_QUANT } from './modules/salmon'
include { SALMON_MERGE } from './modules/salmon'

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
    if (!params.salmon_index && !(params.fasta && params.gtf)) {
        error 'Provide --salmon_index, or both --fasta and --gtf to build one'
    }
    if (params.salmon_index && (params.fasta || params.gtf)) {
        error 'Pass either --salmon_index, or --fasta and --gtf, not both'
    }
    if (!params.salmon_index && !params.tx2gene && !params.gtf) {
        error 'Gene-level counts need --gtf or an explicit --tx2gene table'
    }
    if (params.salmon_index && !params.tx2gene) {
        error 'A prebuilt --salmon_index also requires --tx2gene for gene-level counts'
    }

    reads_ch = readSamplesheet(params.samplesheet)

    if (params.salmon_index) {
        index_ch = channel.fromPath(params.salmon_index, checkIfExists: true, type: 'dir')
        tx2gene_ch = channel.fromPath(params.tx2gene, checkIfExists: true)
    } else {
        fasta_ch = channel.fromPath(params.fasta, checkIfExists: true)
        gtf_ch = channel.fromPath(params.gtf, checkIfExists: true)
        SALMON_REFERENCE(fasta_ch, gtf_ch)
        SALMON_INDEX(SALMON_REFERENCE.out.gentrome, SALMON_REFERENCE.out.decoys)
        index_ch = SALMON_INDEX.out.index
        if (params.tx2gene) {
            tx2gene_ch = channel.fromPath(params.tx2gene, checkIfExists: true)
        } else {
            tx2gene_ch = SALMON_REFERENCE.out.tx2gene
        }
    }

    SALMON_QUANT(index_ch.first(), reads_ch)
    SALMON_MERGE(SALMON_QUANT.out.quant.map { it[1] }.collect(), tx2gene_ch.first())
}
