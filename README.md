# cryVar

A Nextflow pipeline for detecting cryptic variants in SARS-CoV-2 sequencing data.

## Overview

The pipeline processes SARS-CoV-2 sequencing data through the following steps:
1. **Alignment**: Aligns reads to reference genome using minimap2
2. **Sorting & Indexing**: Sorts and indexes BAM files
3. **Primer Trimming**: Removes primer sequences using ivar trim
4. **Variant Calling**: Calls physically linked variants using covar
5. **Cryptic Variant Detection**: Separate script to query outbreak.info clinical API for cryptic variants

## Prerequisites

- Nextflow (>= 25.10.0)
- Conda/Mamba (for conda installation) or Docker (for containerized execution)
- GISAID account for accessing outbreak.info clinical data

## Installation

### Option 1: Conda/Mamba Environment

1. Create the conda environment from the provided `environment.yml`:
```bash
conda env create -f environment.yml
conda activate cryvar
```

Or with mamba (faster):
```bash
mamba env create -f environment.yml
conda activate cryvar
```

## Input Files

### Required Files

1. **FASTQ files**: Sequencing reads (single-end or paired-end)
   - Place in `data/fastq/` directory (default)
   - Naming convention: `{sample_id}.fastq.gz` or `{sample_id}_{R1,R2}.fastq.gz` for paired-end

2. **Metadata CSV file**: Required with the following columns:
   - `sample_id`: Sample identifier (must match FASTQ file names)
   - `primer_scheme`: Primer scheme name (must match a `.bed` file in `reference/primer/`)
   
   Additional metadata, such as collection date and location info are optional

   Example `data/metadata/sample_metadata.csv`:
   ```csv
   sample_id,primer_scheme,collection_date
   SRR36112015,ARTICv5.3.2,2023-01-15
   ```

3. **Reference files**:
   - Reference genome: `reference/reference_genome/NC_045512_Hu-1.fasta`
   - Annotation file: `reference/annotation/NC_045512_Hu-1.gff`
   - Primer BED files: `reference/primer/{primer_scheme}.bed`

## Usage

### GISAID Authentication
In order to access the outbreak.info clinical API, you must authenticate with GISAID.

Run authentication:
```bash
python gisaid_authemtication.py
```
This will open a web browser and prompt you to enter your GISAID credentials. Upon doing so, an auth token will be stored in your environment, so this only needs to be done once.

### Basic Usage

Run the pipeline with a metadata file:
```bash
nextflow run main.nf --metadata data/metadata/sample_metadata.csv
```
Detect cryptic variants from linked mutation output

```bash
python detect_cryptic.py --covar_dir results/covar --metadata data/metadata/sample_metadata.csv
```

### With Custom Parameters

```bash
nextflow run main.nf \
  --metadata data/metadata/sample_metadata.csv \
  --reads "data/fastq/*.fastq.gz" \
  --outdir results \
  --minreadlen 80
```


## Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--metadata` | `null` | **Required**. Path to metadata CSV file |
| `--reads` | `"data/fastq/*.fastq.gz"` | Glob pattern for input FASTQ files |
| `--reference` | `"reference/reference_genome/NC_045512_Hu-1.fasta"` | Path to reference genome FASTA |
| `--annotation` | `"reference/annotation/NC_045512_Hu-1.gff"` | Path to annotation GFF file |
| `--primer_dir` | `"reference/primer"` | Directory containing primer BED files |
| `--outdir` | `"results"` | Output directory |
| `--minreadlen` | `80` | Minimum read length after trimming |
| `--ivar_options` | `""` | Additional options for ivar trim |
| `--covar_options` | `""` | Additional options for covar |

## Output Structure

### Pipeline outputs
```
results/
├── sorted_post_trim/       # Final sorted BAM files
│   ├── {sample_id}.final.sorted.bam
│   └── {sample_id}.final.sorted.bam.bai
└── covar/                  # Variant calling results
|   └── {sample_id}.covar.tsv
|__ detect_cryptic/
    |__ covar_clinical_detections.tsv # All physically linked mutations with number of clinical detects
    |__ cryptic_variants.tsv # Same as above, but filtered for cryptic variants (<= 10 clinical detects, additional QC filters)
```

### covar_clinical_detections.tsv
| Column           | Description                                                 |
| ---------------- | ------------------------------------------------------------|
| `nt_mutations`   | Nucleotide mutations for this cluster                       |
| `aa_mutations`   | Corresponding amino acid translations (where possible*)     |
| `cluster_depth`  | Total number of read pairs with this cluster of mutations   |
| `total_depth`    | Total number of reads spanning this cluster                 |
| `frequency`      | Mutation frequency (cluster depth / total depth)            |
| `coverage_start` | Maximum read start site for which this cluster was detected |
| `coverage_end`   | Minimum read end site for which this cluster was detected   |
| `query`          | Raw mutation query submitted to outbreak.info API           |
| `num_clinical_detections` | Number of clinical detections for this mutation cluster |