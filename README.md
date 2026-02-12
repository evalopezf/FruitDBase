# FruitDBase: Genomic and Transcriptomic Data Repository for *Prunus* Species

A web-based thematic repository for *Prunus* species that enables visualization and download of transcriptomic and genomic data, of *Prunus dulcis*, *Prunus persica*, *Prunus armeniaca*, and *Prunus domestica*.

![Home view](images/home.png)

**Figure 1.** 

---
## Downloader
The Downloader module provides an interactive table interface designed for seamless exploration and retrieval of transcriptomic metadata and datasets. It serves as a central hub for users to filter, select, and export specific subsets of data, including almond atlas expression data.

**Features:**

- Interactive DataTable
- Column-based sorting and global search
- Multi-field filtering using the filter box (organism, tissue, development stage, etc.)
- Direct links to original data sources
- Export options: Copy, CSV, Excel

![Search view](images/searcher.png)

**Figure 2.** Search interface querying species (*Prunus dulcis*) and cultivar (Nonpareil).

---

## Almond Expression Atlas v1.0
Interactive heatmap visualization of gene expression data from Almond Expression Atlas v1.0

**Features:**
- Multi-tissue expression heatmaps (genes × tissues)
- Expression profiles by tissue and developmental stage for flower bud, fruit, and wood
- Bar plots showing expression levels per tissue with color-coded thresholds
- Dynamic filtering by tissue selection, expression score ranges, and top N genes
- Interactive legend for filtering by expression thresholds
- Data export functionality
-  Gene Ontology (GO) term integration for functional characterization of genes using both *Arabidopsis thaliana* and *Prunus dulcis* annotations.

**Example visualization of expression using mean TPM and the Bgee aproach**

![Atlas visualization example](images/6).png)

**Figure 3.** Given a set of selected genes:
- (1.A) Heatmap of general tissues for the selected specific genes  
- (1.B) Heatmap of tissues with developmental stage for floral bud for the selected specific genes  
- (2) Information about the selected gene  
- (3,4) Barplot and table with expression by general tissue  
- (5) Expression data for the samples used to calculate the expression of the selected tissue
- (6) t-SNE plot by tissue and expression level to complement the barplot
- (7) Barplot of expression by development stage and tissue, with an option to visualize expresion by tissue 


---

## SNP Markers

Access to SNP marker data from Axiom 60K SNP array for *Prunus dulcis*.

**Features:**
- Marker information and genomic positions
- Genotype data visualization
- Filtering by chromosome, gene name or marker name
- Export functionality for downstream analysis

**Example of visualization**
![snpExample](images/snpExample.png)
**Figure 4.** Given a selected gene:
- (1.A) Complete catalog of markers according to the performed query
- (1.B) The alleles associated with each marker are shown for the accessions available in the repository
- (2) Information about the selected gene, including GO terms, position, other markers associated with the gene, and gene expression by tissue

---

## Data Sources

- **Expression data:** 
Expression profiles were generated using two approaches: (1) Bgee methodology, and (2) median TPM per tissue. A total of 205 RNA-Seq samples representing all publicly available data for Prunus dulcis were retrieved from the Sequence Read Archive (SRA) obtained using bears R library(Almeida-Silva et al., 2023).  After manual quality filtering and selecting only control conditions, 85 samples were retained for the baseline dataset. An additional 35 samples from an unpublished project were included, resulting in a total of 120 samples for the first almond expression atlas.
​
- **SNP markers:** Data obtained from : Mas-Gómez, J., Gómez-López, F. J., Rubio, M., Gómez-Abajo, M. del M., Dicenta, F., & Martínez-García, P. J. (2025). Integration of linkage mapping, QTL analysis, RNA-Seq data, and Genome-Wide Association Studies (GWAS) to explore relative flowering traits in almond. Horticultural Plant Journal. https://doi.org/10.1016/j.hpj.2025.04.013

- **Functional annotations:** GO terms for *Prunus dulcis* were obtained from the Genome Database for Rosaceae (GDR), and annotations for *Arabidopsis thaliana* were obtained from EnsemblPlants.

---

## Technology

**Backend:** Flask, PostgreSQL  
**Frontend:** D3.js v7, DataTables, Tailwind CSS, JavaScript ES6+

---

## Development Status

**Under active development** - Repository in preparation for public release.

---

## Contact
Department of Plant Breeding, Fruit Breeding Group, CEBAS-CSIC, Espinardo, Murcia, Spain.

For inquiries, please contact: emlopez@cebas.csic.es


---
