nextflow.enable.dsl=2

include { PREPROCESSING } from './modules/preprocess'

process SAMTOOLS_SORT_INDEX_POST {
    tag "$sample_id"
    publishDir "${params.outdir}/sorted_post_trim", mode: 'copy'

    input:
    tuple val(sample_id), path(trimmed_bam)

    output:
    tuple val(sample_id), path("${sample_id}.final.sorted.bam")
    tuple val(sample_id), path("${sample_id}.final.sorted.bam.bai")

    """
    samtools sort -o ${sample_id}.final.sorted.bam ${trimmed_bam}
    samtools index ${sample_id}.final.sorted.bam
    """
}

workflow {
    reads_ch = Channel
        .fromFilePairs(params.reads, flat: false)
        .map { sample_id, reads -> tuple(sample_id, reads) }

    reference_ch = Channel.value(file(params.reference))
    primers_ch   = Channel.value(file(params.primer_bed))

    trimmed_ch = PREPROCESSING(reads_ch, reference_ch, primers_ch)
    final_bams = SAMTOOLS_SORT_INDEX_POST(trimmed_ch)

    final_bams.view { sample_id, bam, bai -> "Final BAM for ${sample_id}: ${bam}" }
}

