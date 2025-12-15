process COVAR {
    tag "$sample_id"
    publishDir "${params.outdir}/covar", mode: 'copy'

    input:
    tuple val(sample_id), path(bam), path(bai)
    path reference
    path annotation

    output:
    tuple val(sample_id), path("${sample_id}.covar.tsv")

    """
    covar \\
    -i ${bam} \\
    -r ${reference} \\
    -a ${annotation} \\
    -t ${task.cpus} \\
    -o ${sample_id}.covar.tsv \\
    ${params.covar_options}
    """
}

process DETECT_CRYPTIC {
    publishDir "${params.outdir}/detect_cryptic", mode: 'copy'

    input:
    path(covar_files)
    path(metadata_file)
    path(detect_script)

    output:
    path("covar_clinical_detections.tsv")
    path("cryptic_variants.tsv")

    script:
    // Handle both single file and list of files
    def fileList = covar_files instanceof List ? covar_files : [covar_files]
    def copyCommands = fileList.collect { file -> "cp '${file}' covar_dir/" }.join('\n    ')
    """
    mkdir -p covar_dir
    ${copyCommands}

    python ${detect_script} \\
    --covar_dir covar_dir \\
    --metadata ${metadata_file} \\
    --output_dir . | tee detect_cryptic.log >&2
    """
}

workflow CRYPTIC_VARIANT_DETECTION {
    take:
    final_bams
    reference_ch
    annotation_ch
    metadata_file
    detect_script

    main:
    covar_ch = COVAR(final_bams, reference_ch, annotation_ch)
    covar_files_ch = covar_ch
        .map { sample_id, covar_file -> covar_file }
        .collect()
    detect_output = DETECT_CRYPTIC(
        covar_files_ch,
        metadata_file,
        detect_script
    )

    emit:
    clinical_detections = detect_output[0]
    cryptic_variants = detect_output[1]
}