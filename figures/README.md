# Figures

R scripts used to generate the figures of the paper. Each script has its parameters at the top; edit them to point to your own data.

## Figures

| Script | Figure | Description | Input | Output |
|---|---|---|---|---|
| [Figure2_sample_distribution_plot.R](Figure2_sample_distribution_plot.R) | Figure 2 | Geographic and biological characterization of the 305 *Prunus dulcis* RNA-Seq samples: (A) world map of projects, (B) projects by country and number of samples, (C-D) samples per tissue and cultivar. | Curated metadata (xlsx) | `plotDistributionMetadata.svg` |
| [Figure3_tsne_plot.R](Figure3_tsne_plot.R) | Figure 3 | t-SNE of the samples by BioProject and tissue (log-normalized counts, top 5000 HVGs). | Metadata (tsv), gene counts (`counts_gene.rds`, from `run_bgeecall.R`) | `tsne_5000_hvg.png` |
| [Figure4_similatiry_plot.R](Figure4_similatiry_plot.R) | Figure 4 | Overlap (Szymkiewicz–Simpson) coefficient between the Bgee and median TPM approaches, and relation between number of samples and number of specific genes detected. | Output of `utils/specific_genes_detection.R` | `overlap_coefficient_plot.svg`, `relation_num_sample_genes_detected.svg`, `combined_plot.svg` |
| [Figure5_goterms_plot.R](Figure5_goterms_plot.R) | Figure 5 | GO semantic similarity networks of the top 15 significant GO terms from the enrichment of tissue-specific genes, per tissue. | Output of `utils/functional_analysis.R` (`tpm_enrichment.rds`, `bgee_enrichment.rds`) | `figure_5_goclusterplot.svg` |
| [Figure8_usecase2.R](Figure8_usecase2.R) | Figure 8 | Expression of DAM genes across flower bud and wood developmental stages. | Expression score by stage (from `bgeecall.R`), `median_tpm_by_tissue_dev.rds` | `dam_genes_expression_pattern.svg` |

## Utils

Scripts in [utils/](utils) process the pipeline results into the intermediate files used by the figures. Run them before the figures that depend on them.

- [specific_genes_detection.R](utils/specific_genes_detection.R): identifies tissue-specific and developmental stage-specific genes and calculates the median TPM by tissue and developmental stage.
  - Input: metadata (tsv), expression score matrix, presence/absence calling and TPM matrix (from `bgeecall.R`).
  - Output: `highlySpecificGeneBgee.rds`, `highlySpecificGeneMeanTPM.rds`, `median_tpm_by_tissue.rds`, `median_tpm_by_tissue_dev.rds`.
  - Used by: Figures 4, 5 (through `functional_analysis.R`) and 8.
- [functional_analysis.R](utils/functional_analysis.R): GO (Biological Process) enrichment of the tissue-specific genes, mapped to their *Arabidopsis* orthologs.
  - Input: `highlySpecificGeneBgee` and `highlySpecificGeneMeanTPM` (from `specific_genes_detection.R`), atlas (from `bgeecall.R`).
  - Output: `tpm_enrichment.rds`, `bgee_enrichment.rds`.
  - Used by: Figure 5.

## Order of execution

```
integration_gf → utils/specific_genes_detection.R → utils/functional_analysis.R → Figures 4, 5, 8
metadata / bgeecall.R outputs → Figures 2, 3
```
