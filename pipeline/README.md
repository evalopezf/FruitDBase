# Pipeline Overview

This folder contains a Snakemake workflow for RNA-seq quantification.

## Table of contents

- [Workflow](#workflow)
- [Inputs](#inputs)
- [Outputs](#outputs)
- [Run](#run)
- [Notes](#notes)
- [Metadata fields](#metadata-fields)
- [Sample integration (post-Snakemake)](#sample-integration-post-snakemake)


## Workflow

The pipeline follows three main steps:

1. **Fastp** trims and quality-filters the raw FASTQ files.
2. **SortMeRNA** removes rRNA reads using the 18S and 28S references.
3. **Kallisto** quantifies the filtered reads against the prepared transcript index.

Snakemake connects these steps automatically, so each sample moves through the workflow based on the files produced by the previous rule.

## Inputs

- Raw FASTQ files downloaded `scritps/download_fastq.sh`
- A metadata table with sample information, columns required: `external_id_sample` (Run), `Layout`, `tissue_name`, `development_stage`, `BioSample` (used to collapse technical replicates and group samples), `url`(used to download raw data)
> *Recommendation*: Curate metadata table beforehand.
- Reference files and a Kallisto index paths defined in `config/config.yaml`

> *Advise* 
>
> - To obtain GTF annotation with intergenic regions we used `pipeline/preprocessing/intergenic_gtf_creation.R`
>
> - To obtain transcriptome FASTA `gffread -w <transcriptome_fasta> -g <genome_fasta> <gtf_annotation_with_intergenics>`
>
> - To obtain kallisto index `kallisto index -i <output_index> <transcriptome_fasta>`

## Outputs

- Fastp QC reports
- Kallisto quantification results
- Summary reports of pseudoaligment

## Run

Update the paths in `config/config.yaml` and then run Snakemake from this directory.

```bash
snakemake --rerun-incomplete --cores 1 --cluster-config cluster.yaml --cluster "sbatch -p {cluster.partition} -t {cluster.time} --mem {cluster.mem} --cpus-per-task {cluster.cpus} -o {cluster.output} -e {cluster.error}" --jobs 4 
```

## Notes

- The workflow supports paired-end samples.
- Batch size and other options can be controlled from the configuration file.

## Metadata fields

Full column reference for the metadata table. 
Columns actually read by the pipeline code are marked *(used)*.

| Column | Meaning |
| --- | --- |
| `external_id_sample` | *(used)* Run accession in SRA/NCBI. Sample identifier matching the FASTQ files. |
| `Origin` | Repository/database the raw data was obtained from (e.g. SRA, ENA, in-house sequencing). |
| `IDBiologicalSampleFruitDBase` | *(used)* Like BioSample in SRA/NCBI. Biological sample identifier; runs that share this value are technical replicates of the same biological sample, and the integration script keeps only the one with the most processed reads. |
| `developmental_stage` | *(used)* Developmental stage of the sample. |
| `origin_country` | Country where the biological material was collected or grown. |
| `origin_city` | City/locality where the biological material was collected or grown. |
| `url` | Download URL for the raw sequencing data. |
| `project_title` | Title of the FruitDBase project the sample was curated under. |
| `origin` | Free-text description of the biological origin of the sample (e.g. cultivar, orchard, collection site). |
| `BioSample` | NCBI BioSample accession. |
| `Experiment` | SRA Experiment accession (SRX). |
| `tissue_name` | *(used)* Tissue/organ label. Groups samples for the expression atlas and expression score "by tissue". |
| `Pubmed` | PubMed ID of the associated publication, if any. |
| `external_id` | Bioproject in SRA NCBI / ID of samples' project |
| `Instrument` | Sequencing instrument model. |
| `Layout` | *(used)* Sequencing layout: `SINGLE` or `PAIRED`. |
| `Selection_method` | Library selection method (e.g. cDNA, RANDOM, PCR). |
| `SRA_sample` | SRA Sample accession (SRS). |
| `SRA_study` | SRA Study accession. |
| `treatment` | Experimental treatment applied to the sample, if any. |
| `accession_name` | Name of the cultivar/accession/variety. |
| `Study_abstract` | Abstract text of the associated study. |
| `project_title` | Title of the BioProject in NCBI SRA . |
| `Date` | Sample collection or submission date. |



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
  with `IDBiologicalSampleFruitDBase`, `tissue_name`, and `development_stage`
  columns added (used to collapse technical replicates and group samples —
  see [Metadata fields](#metadata-fields)).
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
