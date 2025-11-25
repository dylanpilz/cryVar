# cryVar
SARS-CoV-2 Cryptic Variant calling pipeline

## Sample metadata

The pipeline now requires a metadata file supplied via `--metadata /path/to/file`.  
The file must be CSV with a header row and the following columns:

| column          | description                                                                 |
|-----------------|-----------------------------------------------------------------------------|
| `sample_id`     | Prefix shared by the paired FASTQ files (matches `fromFilePairs` output).   |
| `primer_scheme` | Primer scheme name; must match a `.bed` file inside `reference/primer`.     |
| `collection_date` | Sample collection date in `YYYY-MM-DD` format.                           |

Example:

```
sample_id,primer_scheme,collection_date
sampleA,ARTICv5.3.2,2024-03-11
```

At runtime the workflow looks up `<primer_scheme>.bed` under `params.primer_dir`
and feeds the correct scheme into `IVAR_TRIM` for each sample.
