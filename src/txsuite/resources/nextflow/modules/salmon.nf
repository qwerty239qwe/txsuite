process SALMON_REFERENCE {
    tag 'gentrome'
    label 'salmon'
    container params.salmon_image
    publishDir "${params.outdir}/reference", mode: 'copy', pattern: 'tx2gene.tsv'

    input:
    path fasta
    path gtf

    output:
    path 'gentrome.fa', emit: gentrome
    path 'decoys.txt', emit: decoys
    path 'tx2gene.tsv', emit: tx2gene

    script:
    """
    gffread -w txome.fa -g ${fasta} ${gtf}
    gffread --table transcript_id,gene_id ${gtf} \
        | awk -F '\\t' 'NF == 2 && \$1 != "" && \$2 != "" && !seen[\$1]++' > tx2gene.tsv
    test -s tx2gene.tsv

    # A decoy-aware index treats the whole genome as background, which is the
    # correct default for bulk selective alignment: it suppresses spurious
    # transcript assignments for reads originating outside the transcriptome.
    grep '^>' ${fasta} | cut -d ' ' -f 1 | sed 's/^>//' > decoys.txt
    cat txome.fa ${fasta} > gentrome.fa
    """

    stub:
    """
    printf '>tx1\\nACGTACGTAC\\n>tx2\\nTGCATGCATG\\n' > gentrome.fa
    printf 'chr1\\n' > decoys.txt
    printf 'tx1\\tgene1\\ntx2\\tgene2\\n' > tx2gene.tsv
    """
}

process SALMON_INDEX {
    tag 'index'
    label 'salmon'
    container params.salmon_image
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path gentrome
    path decoys

    output:
    path 'salmon_index', emit: index

    script:
    def gencode = params.salmon_gencode ? '--gencode' : ''
    """
    salmon index \
        --transcripts ${gentrome} \
        --decoys ${decoys} \
        --index salmon_index \
        --kmerLen ${params.salmon_kmer_len} \
        --threads ${task.cpus} ${gencode}
    """

    stub:
    """
    mkdir -p salmon_index
    printf 'stub\\n' > salmon_index/info.json
    """
}

process SALMON_QUANT {
    tag "$sample"
    label 'salmon'
    container params.salmon_image
    publishDir "${params.outdir}/quant", mode: 'copy'

    input:
    path index
    tuple val(sample), path(reads)

    output:
    tuple val(sample), path(sample), emit: quant

    script:
    def read_list = reads instanceof List ? reads : [reads]
    def input_reads = read_list.size() > 1 ? "-1 ${read_list[0]} -2 ${read_list[1]}" : "-r ${read_list[0]}"
    """
    salmon quant \
        --index ${index} \
        --libType ${params.salmon_libtype} \
        ${input_reads} \
        --output ${sample} \
        --threads ${task.cpus} \
        --validateMappings \
        --seqBias \
        --gcBias
    test -s ${sample}/quant.sf
    """

    stub:
    """
    mkdir -p ${sample}
    cat > ${sample}/quant.sf <<STUB
Name	Length	EffectiveLength	TPM	NumReads
tx1	1000	800.0	600.000000	12.0
tx2	500	300.0	400.000000	8.0
STUB
    """
}

process SALMON_MERGE {
    tag 'merge'
    label 'salmon'
    container params.salmon_image
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path quant_dirs, stageAs: 'quant/*'
    path tx2gene

    output:
    path 'counts/gene_counts.tsv', emit: gene_counts
    path 'counts/transcript_counts.tsv', emit: transcript_counts
    path 'counts', emit: counts

    script:
    """
    python /opt/txsuite/merge_quants.py quant/* --tx2gene ${tx2gene} --outdir counts
    """

    stub:
    """
    mkdir -p counts
    printf 'gene_id\\tsample1\\ngene1\\t12\\ngene2\\t8\\n' > counts/gene_counts.tsv
    printf 'gene_id\\tsample1\\ngene1\\t600.000000\\ngene2\\t400.000000\\n' > counts/gene_tpm.tsv
    printf 'transcript_id\\tsample1\\ntx1\\t12\\ntx2\\t8\\n' > counts/transcript_counts.tsv
    printf 'transcript_id\\tsample1\\ntx1\\t600.000000\\ntx2\\t400.000000\\n' > counts/transcript_tpm.tsv
    """
}

process SIMPLEAF_INDEX {
    tag 'index'
    label 'salmon'
    container params.salmon_image
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path fasta
    path gtf

    output:
    path 'simpleaf_index', emit: index

    script:
    """
    export ALEVIN_FRY_HOME=\${ALEVIN_FRY_HOME:-\$PWD/.alevin-fry-home}
    mkdir -p \$ALEVIN_FRY_HOME
    simpleaf set-paths
    simpleaf index \
        --output simpleaf_index \
        --fasta ${fasta} \
        --gtf ${gtf} \
        --rlen ${params.simpleaf_rlen} \
        --threads ${task.cpus}
    """

    stub:
    """
    mkdir -p simpleaf_index/index simpleaf_index/ref
    printf 'tx1\\tgene1\\tS\\ntx2\\tgene2\\tS\\n' > simpleaf_index/index/t2g_3col.tsv
    """
}

process SIMPLEAF_QUANT {
    tag "$sample"
    label 'salmon'
    container params.salmon_image
    publishDir "${params.outdir}/quant", mode: 'copy'

    input:
    path index
    tuple val(sample), path(reads)

    output:
    tuple val(sample), path(sample), emit: quant

    script:
    def read_list = reads instanceof List ? reads : [reads]
    def barcode_reads = read_list[0]
    def cdna_reads = read_list.size() > 1 ? read_list[1] : read_list[0]
    def permit = params.alevin_whitelist ? "--unfiltered-pl ${params.alevin_whitelist}" : '--unfiltered-pl'
    """
    export ALEVIN_FRY_HOME=\${ALEVIN_FRY_HOME:-\$PWD/.alevin-fry-home}
    mkdir -p \$ALEVIN_FRY_HOME
    simpleaf set-paths
    simpleaf quant \
        --index ${index}/index \
        --reads1 ${barcode_reads} \
        --reads2 ${cdna_reads} \
        --chemistry ${params.alevin_chemistry} \
        --resolution ${params.alevin_resolution} \
        ${permit} \
        --t2g-map ${index}/index/t2g_3col.tsv \
        --threads ${task.cpus} \
        --output ${sample}
    test -s ${sample}/af_quant/alevin/quants_mat.mtx
    """

    stub:
    """
    mkdir -p ${sample}/af_quant/alevin
    cat > ${sample}/af_quant/alevin/quants_mat.mtx <<STUB
%%MatrixMarket matrix coordinate real general
%
2 6 3
1 1 3
2 2 5
1 5 1
STUB
    printf 'AAACCTGAGAAACCAT\\nAAACCTGAGAAACCGC\\n' > ${sample}/af_quant/alevin/quants_mat_rows.txt
    printf 'gene1-S\\ngene2-S\\ngene1-U\\ngene2-U\\ngene1-A\\ngene2-A\\n' > ${sample}/af_quant/alevin/quants_mat_cols.txt
    printf '{"usa_mode": true}\\n' > ${sample}/af_quant/alevin/quants_mat.json
    """
}

process ALEVIN_TO_H5AD {
    tag 'h5ad'
    label 'single_cell'
    container params.single_cell_image
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path quant_dirs, stageAs: 'quant/*'

    output:
    path 'matrix/alevin.h5ad', emit: matrix

    script:
    """
    python /opt/txsuite/alevin_to_h5ad.py quant/* --output matrix/alevin.h5ad
    """

    stub:
    """
    mkdir -p matrix
    touch matrix/alevin.h5ad
    """
}
