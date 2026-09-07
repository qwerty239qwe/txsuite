nextflow.enable.dsl = 2

include { PSEUDOBULK_DE } from './subworkflows/pseudobulk_de'
include { COLLECT_DE } from './modules/collect_de'

def simpleColumn(value, label) {
    if (!(value ==~ /^[A-Za-z][A-Za-z0-9_.]*$/)) {
        error "${label} must be a simple column name"
    }
    value
}

def simpleLevel(value, label) {
    if (!(value ==~ /^[A-Za-z0-9][A-Za-z0-9_.-]*$/)) {
        error "${label} must be a simple factor level"
    }
    value
}

def comparisonMeta(row, sampleColumn, countsLayer) {
    ['comparison', 'design', 'reference', 'test'].each { key ->
        if (!row[key]?.trim()) error "Comparison manifest has an empty ${key}"
    }
    def groupColumn = row.group_column?.trim() ?: ''
    def groupValue = row.group_value?.trim() ?: ''
    if (!!groupColumn != !!groupValue) error 'Group column and value must be used together'
    def method = row.method?.trim() ?: 'deseq2'
    if (!(method in ['deseq2', 'edger', 'limma'])) error "Unsupported DE method: ${method}"
    def covariates = row.covariates?.trim() ?: ''
    if (covariates && !(covariates ==~ /^[A-Za-z][A-Za-z0-9_.]*(,[A-Za-z][A-Za-z0-9_.]*)*$/)) {
        error 'Covariates must be comma-separated simple column names'
    }
    def padj = row.padj?.trim() ? row.padj as Double : 0.05
    def lfc = row.lfc?.trim() ? row.lfc as Double : 1.0
    def topGenes = row.top_genes?.trim() ? row.top_genes as Integer : 50
    if (!(padj > 0 && padj <= 1) || lfc < 0 || topGenes < 1) {
        error 'padj must be in (0, 1], lfc non-negative, and top_genes positive'
    }
    def reference = simpleLevel(row.reference.trim(), 'Reference')
    def test = simpleLevel(row.test.trim(), 'Test')
    if (reference == test) error 'Reference and test levels must differ'
    [
        id: simpleLevel(row.comparison.trim(), 'Comparison'),
        sample_column: simpleColumn(sampleColumn, 'Sample column'),
        counts_layer: simpleColumn(countsLayer, 'Counts layer'),
        group_column: groupColumn ? simpleColumn(groupColumn, 'Group column') : '',
        group_value: groupValue ? simpleLevel(groupValue, 'Group value') : '',
        design: simpleColumn(row.design.trim(), 'Design'),
        reference: reference,
        test: test,
        method: method,
        covariates: covariates,
        padj: padj,
        lfc: lfc,
        top_genes: topGenes,
    ]
}

workflow {
    if (!params.input) error 'Missing required parameter: --input'
    if (!params.single_cell_image || !params.bulk_image) {
        error 'Both --single_cell_image and --bulk_image are required'
    }
    def h5ad = file(params.input, checkIfExists: true)

    if (params.manifest) {
        def allowed = ['comparison', 'group_column', 'group_value', 'design', 'reference',
                       'test', 'method', 'covariates', 'padj', 'lfc', 'top_genes'] as Set
        comparisons = Channel.fromPath(params.manifest, checkIfExists: true)
            .splitCsv(header: true, sep: '\t', strip: true)
            .toList()
            .flatMap { rows ->
                if (!rows) error 'Comparison manifest has no data rows'
                def unknown = rows[0].keySet() - allowed
                if (unknown) error "Comparison manifest has unknown columns: ${unknown.sort().join(', ')}"
                def metas = rows.collect { row -> comparisonMeta(row, params.sample_column, params.counts_layer) }
                def duplicates = metas.groupBy { it.id.toLowerCase() }.findAll { key, values -> values.size() > 1 }
                if (duplicates) error "Duplicated comparison: ${duplicates.keySet().sort().join(', ')}"
                metas.collect { meta -> tuple(meta, h5ad) }
            }
    } else {
        def row = [
            comparison: params.comparison ?: "${params.test}_vs_${params.reference}",
            group_column: params.group_column ?: '', group_value: params.group_value ?: '',
            design: params.design, reference: params.reference, test: params.test,
            method: params.method ?: 'deseq2', covariates: params.covariates ?: '',
            padj: params.padj as String, lfc: params.lfc as String,
            top_genes: params.top_genes as String,
        ]
        comparisons = Channel.of(tuple(comparisonMeta(row, params.sample_column, params.counts_layer), h5ad))
    }

    expected = comparisons.map { meta, input -> meta }.collect()
    PSEUDOBULK_DE(comparisons)
    // toList emits [] when all comparisons fail, so collection still runs.
    result_inputs = PSEUDOBULK_DE.out.results
        .map { meta, directory -> directory }
        .toList()
        .map { directories -> directories.sort { a, b -> a.name <=> b.name } }
    COLLECT_DE(expected, result_inputs, file("${projectDir}/../single_cell_python/collect_de.py"))
}
