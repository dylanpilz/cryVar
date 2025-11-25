# cryVar
SARS-CoV-2 Cryptic Variant calling pipeline

## Sample metadata

The pipeline requires a metadata file supplied via `--metadata /path/to/file`.  

| column          | description                                                                 |
|-----------------|-----------------------------------------------------------------------------|
| `sample_id`     | Prefix shared by the paired FASTQ files.   |
| `primer_scheme` | Primer scheme name; must match a `.bed` file inside `reference/primer`.     |
| `collection_date` | Sample collection date in `YYYY-MM-DD` format.                           |

Example:

```
sample_id,primer_scheme,collection_date
sampleA,ARTICv5.3.2,2024-03-11
```

