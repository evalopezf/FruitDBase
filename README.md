# FruitDBase: Genomic Atlas for *Prunus* Species

Web platform for genomic data visualization and breeding management in almond, peach, apricot, and plum.

---

## Gene Expression Atlas

Interactive visualization of tissue-specific gene expression data from the Bgee database. Features include:

- **Multi-tissue heatmaps** with expression profiling across developmental stages
- **Dynamic filtering** by expression levels, tissue types, and gene identifiers  
- **Gene-centric views** with bar plots and statistical summaries
- **Integration with GO annotations** for functional characterization

![Expression Atlas Screenshot](docs/images/expression_atlas.png)

---

## Breeding Tools

### Pedigree Visualization
Force-directed network graphs displaying family relationships, ancestry, and breeding crosses across germplasm collections.

![Pedigree Network](docs/images/pedigree_network.png)

### Germplasm Browser
Searchable database of accessions with phenological traits, S-allele compatibility data, and breeding metadata.

---

## Genetic Markers Explorer *(In Development)*

Browse and filter SNPs, QTLs, and SSR markers with exportable datasets and detailed marker annotations.

---

## Technical Stack

**Backend:** Flask, SQLAlchemy, PostgreSQL  
**Frontend:** D3.js, DataTables, Tailwind CSS  
**External Data:** Bgee Expression Database, QuickGO

---

## Citation

*Under review - citation information will be available upon publication.*

---

## Contact

For questions or collaboration inquiries, please contact the development team.
- **Presence**: Raw expression measurements from RNA-seq
- **Sample**: RNA-seq and WGS sample metadata
- **Organism**: Prunus species taxonomy

## Key Visualizations

### Expression Heatmap
- **Layout**: Tissues (rows) × Genes (columns)
- **Bar Plots**: Expression level per tissue with color coding
- **Interactivity**: 
  - Drag legend to filter by expression score
  - Click tissue names to expand/collapse
  - Toggle scientific notation
- **Filters**: Tissue selection, expression range (min/max), top N genes

### Pedigree Network
- **Node Types**:
  - Purple (large) = Breeding families
  - Blue = Sires (♂)
  - Pink = Dams (♀)
  - Colored = Offspring
- **Interactions**:
  - Drag nodes to rearrange
  - Click nodes to view complete pedigree tree
  - Hover for family statistics
- **Tree View**: Recursive rendering of ancestors (above) and descendants (below)

## User Roles

- **Admin**: Full database access and data management
- **Group Member**: Access to group-specific datasets
- **Public**: Read-only access to published data

## Development Status

🚧 **In Active Development** - This project is currently under development for internal use and version control.

## Contact

**CEBAS-CSIC Almond Breeding Program**  
Murcia, Spain

---

*Built for advancing Prunus genomics research* 🌰
