nextflow.enable.dsl = 2

include { PSEUDOBULK } from './modules/pseudobulk'
include { DESEQ2 } from './modules/deseq2'

workflow {
    if (!params.input) {
        error 'Missing required parameter: --input'
    }
    if (!(params.sample_column ==~ /^[A-Za-z][A-Za-z0-9_.]*$/)) {
        error '--sample_column must be a simple column name'
    }
    if (!(params.design ==~ /^[A-Za-z][A-Za-z0-9_.]*$/)) {
        error '--design must be a simple column name'
    }
    if (!(params.reference ==~ /^[A-Za-z0-9][A-Za-z0-9_.-]*$/) ||
        !(params.test ==~ /^[A-Za-z0-9][A-Za-z0-9_.-]*$/)) {
        error '--reference and --test must be simple factor levels'
    }
    if (params.reference == params.test) {
        error '--reference and --test must differ'
    }

    input_ch = channel.fromPath(params.input, checkIfExists: true)
    PSEUDOBULK(input_ch)
    DESEQ2(PSEUDOBULK.out.counts, PSEUDOBULK.out.metadata)
}
