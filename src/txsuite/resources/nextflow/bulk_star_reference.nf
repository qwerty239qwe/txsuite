nextflow.enable.dsl = 2

include { STAR_GENOME_GENERATE } from './modules/star_reference'
include { SAMTOOLS_FAIDX } from './modules/star_reference'
include { GATK_SEQUENCE_DICTIONARY } from './modules/star_reference'
include { REFERENCE_MANIFEST } from './modules/star_reference'

workflow {
    if (!params.fasta) {
        error 'Missing required parameter: --fasta'
    }
    if (!params.gtf) {
        error 'Missing required parameter: --gtf'
    }
    if (!params.star_sjdb_overhang || params.star_sjdb_overhang.toString().toInteger() < 1) {
        error '--star_sjdb_overhang must be a positive integer'
    }

    fasta_ch = channel.fromPath(params.fasta, checkIfExists: true)
    gtf_ch = channel.fromPath(params.gtf, checkIfExists: true)

    STAR_GENOME_GENERATE(fasta_ch, gtf_ch)
    SAMTOOLS_FAIDX(fasta_ch)
    GATK_SEQUENCE_DICTIONARY(fasta_ch)
    REFERENCE_MANIFEST(
        STAR_GENOME_GENERATE.out.index,
        SAMTOOLS_FAIDX.out.fai,
        GATK_SEQUENCE_DICTIONARY.out.dict,
    )
}
