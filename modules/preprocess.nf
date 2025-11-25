nextflow.enable.dsl=2

process ALIGN_MINIMAP2 {
    tag "$sample_id"
    publishDir "${params.outdir}/alignment", mode: 'copy'

    input:
    tuple val(sample_id), path(reads)
    path reference

    output:
    tuple val(sample_id), path("${sample_id}.unsorted.bam")

    script:
    def read_opts = reads instanceof List && reads.size() == 2
        ? "${reads[0]} ${reads[1]}"
        : reads
    """
    minimap2 -ax sr ${reference} ${read_opts} \\
      | samtools view -bS - \\
      > ${sample_id}.unsorted.bam
    """
}

process SAMTOOLS_SORT_INDEX_PRE {
    tag "$sample_id"
    publishDir "${params.outdir}/sorted_pre_trim", mode: 'copy'

    input:
    tuple val(sample_id), path(bam)

    output:
    tuple val(sample_id), path("${sample_id}.sorted.bam"), path("${sample_id}.sorted.bam.bai")

    """
    samtools sort -o ${sample_id}.sorted.bam ${bam}
    samtools index ${sample_id}.sorted.bam
    """
}

process IVAR_TRIM {
    tag "$sample_id"
    publishDir "${params.outdir}/trimmed", mode: 'copy'

    input:
    tuple val(sample_id), path(sorted_bam), path(sorted_bai)
    path reference
    path primers

    output:
    tuple val(sample_id), path("${sample_id}.trimmed.bam")

    """
    ivar trim \\
      -i ${sorted_bam} \\
      -b ${primers} \\
      -p ${sample_id}.trimmed \\
      -m ${params.minreadlen} \\
      -q ${params.slop} \\
      ${params.ivar_options}

    samtools view -b ${sample_id}.trimmed.bam > ${sample_id}.trimmed.bam.tmp
    mv ${sample_id}.trimmed.bam.tmp ${sample_id}.trimmed.bam
    """
}

workflow PREPROCESSING {
    take:
    reads_ch
    reference_ch
    primers_ch

    main:
    aligned   = ALIGN_MINIMAP2(reads_ch, reference_ch)
    presorted = SAMTOOLS_SORT_INDEX_PRE(aligned)
    trimmed   = IVAR_TRIM(presorted, reference_ch, primers_ch)

    emit:
    trimmed
}

