# Figures

R scripts used to generate the figures of the paper. Each script has its parameters at the top; edit them to point to your own data.

## Figures

<table>
  <thead>
    <tr>
      <th width="20%">Script</th>
      <th width="8%">Figure</th>
      <th width="36%">Description</th>
      <th width="18%">Input</th>
      <th width="18%">Output</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td><a href="Figure2_sample_distribution_plot.R">Figure2_sample_distribution_plot.R</a></td>
      <td>Figure 2</td>
      <td>Geographic and biological characterization of the 305 <i>Prunus dulcis</i> RNA-Seq samples: (A) world map of projects, (B) projects by country and number of samples, (C-D) samples per tissue and cultivar.</td>
      <td>Curated metadata (xlsx)</td>
      <td><code>plotDistributionMetadata.svg</code></td>
    </tr>
    <tr>
      <td><a href="Figure3_tsne_plot.R">Figure3_tsne_plot.R</a></td>
      <td>Figure 3</td>
      <td>t-SNE of the samples by BioProject and tissue (log-normalized counts, top 5000 HVGs).</td>
      <td>Metadata (tsv), gene counts (<code>counts_gene.rds</code>, from <code>snakemake post_processing</code>)</td>
      <td><code>tsne_5000_hvg.png</code></td>
    </tr>
    <tr>
      <td><a href="Figure4_similatiry_plot.R">Figure4_similatiry_plot.R</a></td>
      <td>Figure 4</td>
      <td>Overlap (Szymkiewicz–Simpson) coefficient between the Bgee and median TPM approaches, and relation between number of samples and number of specific genes detected.</td>
      <td>Output of <code>utils/specific_genes_detection.R</code></td>
      <td><code>overlap_coefficient_plot.svg</code>, <code>relation_num_sample_genes_detected.svg</code>, <code>combined_plot.svg</code></td>
    </tr>
    <tr>
      <td><a href="Figure5_goterms_plot.R">Figure5_goterms_plot.R</a></td>
      <td>Figure 5</td>
      <td>GO semantic similarity networks of the top 15 significant GO terms from the enrichment of tissue-specific genes, per tissue.</td>
      <td>Output of <code>utils/functional_analysis.R</code> (<code>tpm_enrichment.rds</code>, <code>bgee_enrichment.rds</code>)</td>
      <td><code>figure_5_goclusterplot.svg</code></td>
    </tr>
    <tr>
      <td><a href="Figure8_usecase2.R">Figure8_usecase2.R</a></td>
      <td>Figure 8</td>
      <td>Expression of DAM genes across flower bud and wood developmental stages.</td>
      <td>Expression score by stage (from <code>snakemake post_processing</code>), <code>median_tpm_by_tissue_dev.rds</code></td>
      <td><code>dam_genes_expression_pattern.svg</code></td>
    </tr>
  </tbody>
</table>

## Utils

Scripts in [utils/](utils) process the pipeline results into the intermediate files used by the figures. Run them before the figures that depend on them.

- [specific_genes_detection.R](utils/specific_genes_detection.R): identifies tissue-specific and developmental stage-specific genes and calculates the median TPM by tissue and developmental stage.
  - Input: metadata (tsv), expression score matrix, presence/absence calling and TPM matrix (from `snakemake post_processing`, see [pipeline/README.md](../pipeline/README.md)).
  - Output: `highlySpecificGeneBgee.rds`, `highlySpecificGeneMeanTPM.rds`, `median_tpm_by_tissue.rds`, `median_tpm_by_tissue_dev.rds`.
  - Used by: Figures 4, 5 (through `functional_analysis.R`) and 8.
- [functional_analysis.R](utils/functional_analysis.R): GO (Biological Process) enrichment of the tissue-specific genes, mapped to their *Arabidopsis* orthologs.
  - Input: `highlySpecificGeneBgee` and `highlySpecificGeneMeanTPM` (from `specific_genes_detection.R`), atlas (from `snakemake post_processing`).
  - Output: `tpm_enrichment.rds`, `bgee_enrichment.rds`.
  - Used by: Figure 5.

