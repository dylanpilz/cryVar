nextflow.enable.dsl=2

include { PREPROCESSING } from './modules/preprocess'

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
    -t 8 \\
    -o ${sample_id}.covar.tsv \\
    ${params.covar_options}
    """
}

workflow {
    if( !params.metadata ) {
        throw new IllegalArgumentException("Provide a metadata file containing sample_id and primer_scheme via --metadata.")
    }

    def metadataFile = file(params.metadata)
    if( !metadataFile.exists() ) {
        throw new FileNotFoundException("Metadata file not found: ${metadataFile}")
    }

    def metadata_ch = Channel
        .fromPath(metadataFile)
        .splitCsv(header: true)
        .map { Map record ->
            def sampleId       = (record.sample_id ?: '').toString().trim()
            def primerScheme   = (record.primer_scheme ?: '').toString().trim()
            def collectionDate = (record.collection_date ?: '').toString().trim()

            def primerFile = file("${params.primer_dir}/${primerScheme}.bed")
            if( !primerFile.exists() ) {
                throw new FileNotFoundException("Primer scheme '${primerScheme}' not found in ${params.primer_dir} (expected ${primerFile}).")
            }

            tuple(sampleId, primerFile)
        }

    // Get all fastq files and group by sample ID
    def all_reads_ch = Channel
        .fromPath(params.reads)
        .map { file ->
            def name = file.getName()
            def sample_id = name
                .replaceAll(/\.fastq\.gz$/, '')
                // Handle Illumina format: _S\d+_L\d+_R[12]_\d+ (e.g., _S8_L001_R1_001)
                .replaceAll(/_S\d+_L\d+_R[12]_\d+$/, '')
                // Handle simple format: _R[12] (e.g., _R1, _R2)
                .replaceAll(/_[R12]$/, '')
            tuple(sample_id, file)
        }
        .groupTuple()
        .map { sample_id, files ->
            def reads
            if (files.size() == 2) {
                def sorted = files.sort { it.getName() }
                // Paired end
                reads = sorted
            } else if (files.size() == 1) {
                // Single-end
                reads = files[0]
            } else {
                throw new IllegalArgumentException("Unexpected number of files (${files.size()}) for sample ${sample_id}: ${files}")
            }
            tuple(sample_id, reads)
        }
    
    // Get sample IDs from metadata and join with reads
    def sample_ids_ch = Channel
        .fromPath(metadataFile)
        .splitCsv(header: true)
        .map { Map record -> (record.sample_id ?: '').toString().trim() }
    
    reads_ch = sample_ids_ch
        .map { sample_id -> tuple(sample_id, null) }
        .join(all_reads_ch, by: 0)
        .map { sample_id, meta, reads -> 
            if (reads == null) {
                throw new FileNotFoundException("No reads found for sample ${sample_id} matching pattern ${params.reads}")
            }
            tuple(sample_id, reads)
        }

    reference_ch = Channel.value(file(params.reference))
    annotation_ch = Channel.value(file(params.annotation))

    final_bams = PREPROCESSING(reads_ch, reference_ch, metadata_ch)
    
    covar_ch = COVAR(final_bams, reference_ch, annotation_ch)
}

