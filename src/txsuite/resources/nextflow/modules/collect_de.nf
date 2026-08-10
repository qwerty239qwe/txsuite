process COLLECT_DE {
    tag 'summary'
    label 'single_cell'

    input:
    val expected
    path result_inputs

    output:
    path 'comparison-index.tsv', emit: index
    path 'combined-results.tsv', emit: combined

    script:
    def rows = expected.collect { meta ->
        [meta.id, meta.group_column, meta.group_value, meta.design, meta.reference,
         meta.test, meta.method].join('\t')
    }.join('\n')
    """
    printf 'comparison\tgroup_column\tgroup_value\tdesign\treference\ttest\tmethod\n${rows}\n' > expected.tsv
    python /opt/txsuite/single_cell.py collect-de expected.tsv . ${result_inputs.join(' ')}
    """

    stub:
    """
    printf 'comparison\tstatus\tmethod\tgroup_column\tgroup_value\tdesign\treference\ttest\tresult\tsignificant\terror\n' > comparison-index.tsv
    printf 'comparison\tgroup_column\tgroup_value\tdesign\treference\ttest\tmethod\n' > combined-results.tsv
    """
}
