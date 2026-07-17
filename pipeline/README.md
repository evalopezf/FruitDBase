# Pipeline Overview

This folder contains a Snakemake workflow for RNA-seq quantification.

## Workflow

The pipeline follows three main steps:

1. **Fastp** trims and quality-filters the raw FASTQ files.
2. **SortMeRNA** removes rRNA reads using the 18S and 28S references.
3. **Kallisto** quantifies the filtered reads against the prepared transcript index.

Snakemake connects these steps automatically, so each sample moves through the workflow based on the files produced by the previous rule.

## Inputs

- Raw FASTQ files downloaded `scritps/download_fastq.sh`
- A metadata table with sample information, columns required: `external_id_sample` (Run), `Layout`
- Reference files and a Kallisto index paths defined in `config/config.yaml`

## Outputs

- Fastp QC reports
- Kallisto quantification results
- Summary reports for the run

## Run

Update the paths in `config/config.yaml` and then run Snakemake from this directory.

```bash
snakemake --cores 4
```

## Notes

- The workflow supports paired-end samples.
- Batch size and other options can be controlled from the configuration file.

## Sample integration (post-Snakemake)

Once `snakemake` has finished, `scripts/integration_pipeline.R` integrates all
samples into a gene-level expression atlas by expression score matrices.

Requires the R packages: `optparse`, `readr`, `dplyr`, `data.table`, `fs`,
`tximport`, `BgeeCall`, `Biostrings`, `GenomicFeatures`, `txdbmaker`,
`AnnotationDbi`, `rtracklayer`.

```bash
Rscript scripts/integration_pipeline.R \
  --metadata newMetadata20260701.tsv \
  --kallisto-report kallisto/results/kallisto_report.rds \
  --kallisto-dir kallisto/results \
  --fasta transcriptome.fasta \
  --gtf annotation.gtf \
  --species-id 3755 \
  --output-dir results/integration
```

Run `Rscript scripts/integration_pipeline.R --help` for the full option list.

### Inputs

- `--metadata`: the same curated metadata table used by the Snakemake run,
  with `BioSample`, `tissue_name`, and `development_stage`
  columns added (used to collapse technical replicates and group samples).
- `--kallisto-report`: `kallisto_report.rds` produced by the `kallisto_report`
  rule (used to collapse technical replicates and group samples).
- `--kallisto-dir`: the Kallisto output directory from Snakemake
  (`{KALLISTO_DIR}` in `config.yaml`) - one subfolder per sample with
  `abundance.tsv`.
- `--fasta` / `--gtf`: the transcriptome FASTA and matching GTF annotation
  used to build the Kallisto index (`kallisto_index` in `config.yaml`). The
  transcript-to-gene map is derived automatically from the GTF at runtime -
  no separate precomputed file is needed.
- `--species-id`: NCBI Taxonomy ID for the species.


### Outputs (under `--output-dir`)

- `presence_absence_df.rds` - per-sample, per-gene presence/absence calls.
- `presence_by_tissue.rds`, `presence_by_dev.rds` - expression atlas grouped
  by tissue, and by tissue x developmental stage.
- `ExpressionScorebyTissue.rds`, `ExpressionScoreByDev.rds` - 
expression score matrices (0-100) per gene and condition.
