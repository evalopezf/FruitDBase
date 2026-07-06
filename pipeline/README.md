# Pipeline Overview

This folder contains a Snakemake workflow for RNA-seq quantification.

## Workflow

The pipeline follows three main steps:

1. **Fastp** trims and quality-filters the raw FASTQ files.
2. **SortMeRNA** removes rRNA reads using the 18S and 28S references.
3. **Kallisto** quantifies the filtered reads against the prepared transcript index.

Snakemake connects these steps automatically, so each sample moves through the workflow based on the files produced by the previous rule.

## Inputs

- Raw FASTQ files
- A metadata table with sample information
- Reference files and a Kallisto index defined in `config/config.yaml`

## Outputs

- Fastp QC reports
- rRNA-filtered reads
- Kallisto quantification results
- Summary reports for the run

## Run

Update the paths in `config/config.yaml` and then run Snakemake from this directory.

```bash
snakemake --cores 4
```

## Notes

- The workflow supports paired-end and single-end samples.
- Batch size and other options can be controlled from the configuration file.
