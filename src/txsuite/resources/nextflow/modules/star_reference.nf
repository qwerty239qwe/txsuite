process STAR_GENOME_GENERATE {
    tag 'star-index'
    label 'star'
    container params.star_image
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path fasta
    path gtf

    output:
    path 'star_index', emit: index

    script:
    def sa_index = params.star_sa_index_nbases ? "--genomeSAindexNbases ${params.star_sa_index_nbases}" : ''
    """
    mkdir -p star_index
    STAR --runMode genomeGenerate \
        --genomeDir star_index \
        --genomeFastaFiles ${fasta} \
        --sjdbGTFfile ${gtf} \
        --sjdbOverhang ${params.star_sjdb_overhang} \
        --runThreadN ${task.cpus} \
        --limitGenomeGenerateRAM ${params.star_memory_gb.toLong() * 1000000000} ${sa_index}
    test -s star_index/SA
    """

    stub:
    """
    mkdir -p star_index
    printf 'stub\\n' > star_index/SA
    printf 'stub\\n' > star_index/Genome
    printf 'versionGenome\\t2.7.4a\\n' > star_index/genomeParameters.txt
    """
}

process SAMTOOLS_FAIDX {
    tag 'faidx'
    label 'star'
    container params.star_image
    publishDir "${params.outdir}/reference", mode: 'copy'

    // GATK resolves the index and dictionary from the reference basename: X.fa
    // needs X.fa.fai and X.dict beside it. These therefore keep the caller's
    // FASTA name, and the planner derives the same names from the fasta input.
    input:
    path fasta

    output:
    path "${fasta}.fai", emit: fai

    script:
    """
    samtools faidx ${fasta}
    """

    stub:
    """
    printf 'chr1\\t1000\\t6\\t60\\t61\\n' > ${fasta}.fai
    """
}

process GATK_SEQUENCE_DICTIONARY {
    tag 'dict'
    label 'star'
    container params.star_image
    publishDir "${params.outdir}/reference", mode: 'copy'

    input:
    path fasta

    output:
    path "${fasta.baseName}.dict", emit: dict

    script:
    """
    gatk CreateSequenceDictionary --REFERENCE ${fasta} --OUTPUT ${fasta.baseName}.dict
    """

    stub:
    """
    printf '@HD\\tVN:1.6\\n@SQ\\tSN:chr1\\tLN:1000\\n' > ${fasta.baseName}.dict
    """
}

process REFERENCE_MANIFEST {
    tag 'manifest'
    label 'star'
    container params.star_image
    publishDir "${params.outdir}/reference", mode: 'copy'

    input:
    path index
    path fai
    path dict

    output:
    path 'reference-manifest.tsv', emit: manifest

    script:
    // The overhang an index was built for cannot be recovered from the index
    // directory, and reusing an index built for a different read length costs
    // junction sensitivity silently. Record it next to the artifacts.
    """
    cat > reference-manifest.tsv <<MANIFEST
key	value
star_index	${index}
fasta_index	${fai}
sequence_dictionary	${dict}
read_length	${params.star_read_length}
sjdb_overhang	${params.star_sjdb_overhang}
MANIFEST
    """

    stub:
    """
    cat > reference-manifest.tsv <<MANIFEST
key	value
star_index	${index}
fasta_index	${fai}
sequence_dictionary	${dict}
read_length	${params.star_read_length}
sjdb_overhang	${params.star_sjdb_overhang}
MANIFEST
    """
}
