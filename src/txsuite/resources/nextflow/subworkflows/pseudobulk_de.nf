include { PSEUDOBULK } from '../modules/pseudobulk'
include { BULK_DE } from '../modules/bulk_de'

workflow PSEUDOBULK_DE {
    take:
    comparisons

    main:
    PSEUDOBULK(comparisons)
    BULK_DE(PSEUDOBULK.out.data)

    emit:
    pseudobulk = PSEUDOBULK.out.data
    results = BULK_DE.out.results
}
