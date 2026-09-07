process COLLECT_DE {
    tag 'summary'
    label 'single_cell'

    input:
    val expected
    path result_inputs
    path collector

    output:
    path 'comparison-index.tsv', emit: index
    path 'combined-results.tsv', emit: combined

    script:
    def rows = expected.collect { meta ->
        [meta.id, meta.group_column, meta.group_value, meta.design, meta.reference,
         meta.test, meta.method].join('\t')
    }.join('\n')
    def directories = result_inputs.collect { "'${it}'" }.join(' ')
    """
    printf 'comparison\tgroup_column\tgroup_value\tdesign\treference\ttest\tmethod\n${rows}\n' > expected.tsv
    python3 "${collector}" expected.tsv . ${directories}
    """

}
