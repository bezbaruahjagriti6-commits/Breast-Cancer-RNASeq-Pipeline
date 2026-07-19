#  Comparative Transcriptomic Analyis of Breast Cancer Subtypes: Basal vs. Luminal A

An end-to-end RNA-Seq bioinformatics pipeline to identify differentially expressed genes (DEGs) and enriched biological pathways distinguishing Basal (Triple-Negative) and Luminal A breast cancer subtypes.

## Project Overview

Breast cancer is a heterogeneous disease comprising molecular subtypes with distinct gene expression patterns, prognosis, and treatment responses. This project compares the aggressive Basal/TNBC subtype with Luminal A to identify subtype-specific biomarkers and biological pathways using RNA-Seq analysis.

## Project Goal

This project processes raw RNA-Seq reads to perform quality control, alignment, quantification, and differential expression analysis. The final output identifies key genes and biological pathways that distinguish the aggressive Basal subtype from the hormone-responsive Luminal A subtype.

-   **Data Source:** NCBI GEO (`GSE58135`/ BioProject: PRJNA251383)
-   **Comparison:** 5 Basal vs. 5 Luminal A primary tumor samples
-   **Reference Genome:** GRCh38


## Pipeline Architecture and Toolkits

Raw FASTQ Reads -> FastQC (QC Assessment) -> Trimmomatic (Adapter/Quality Trim) -> HISAT2 (Splice-Aware GRCh38 Alignment) -> featureCounts (Abundance Quantification) -> Raw Count Matrix -> DESeq2 (Wald-Test Modeling) -> Significance Subsetting (Log2FC > 1, padj < 0.05) -> Visualizations (PCA, MA, Volcano, Biclustered Heatmap) -> clusterProfiler (GO and KEGG Pathway Enrichment)

## Core Analytical Results

* **Principal Component Analysis (PCA):** PC1 explained 32% of the total variance, clearly separating Basal and Luminal A samples, indicating distinct transcriptomic profiles.

* **MA Diagnostic Plot:** Symmetrical distribution of log2 fold-changes centered exactly across the baseline 
axis, validating baseline count matrix normalization across all sequencing library depths.

# Key Results

| **Initial Count Features** | 61,552 |
| **Low-Count Filtered Features** | 27,875 |
| **Statistically Confirmed DEGs (Absolute Log2FC >= 1, padj < 0.05)** | **3,906** |

### Key validated biomarkers identified include:

  **Up-regulated in Basal (Aggressive subtype):**
      `FOXC1`: Master transcription factor driving EMT.
      `MSLN`: Cell-surface oncogene; a clear target for Antibody-Drug Conjugates (ADCs).
      `CT83`: Cancer/Testis antigen; an ideal target for CAR-T therapy.
      `ART3`: Enzyme that hyPer-activates Akt and ERK proliferative pathways in TNBC.
      `TLX3`: Oncogenic transcription factor implicated in cell proliferation.
      `HORMAD1`: Cancer-testis antigen linked to genomic instability.
      `CRABP1`: Retinoic acid signaling protein associated with poor prognosis.
  **Down-regulated in Basal (Luminal A Enriched):**
      `GATA3`: Pioneer factor essential for Estrogen Receptor binding.
      `KRT8 / KRT18`: Definitive cytokeratin markers for luminal architecture.
      `DSCAM-AS1`: Estrogen-Receptor-targeted lncRNA driving endocrine resistance.
      `SPAG6`: Structural protein implicated in cilia function and cell signaling.
      `DOK7`: Downstream kinase anchor protein involved in signaling pathways.

### Downstream Visualizations

### Heatmap
Hierarchical clustering successfully separated Basal and Luminal samples using subtype-specific marker genes, demonstrating clear expression difference. The **Basal block** exhibits uniform activation of aggressive stem-like markers (CT83, TLX3, ART3, FOXC1, HORMAD1, CRABP1), while the **Luminal block** shows co-expression of pioneer elements and structural intermediate filaments(GATA3, KRT8, KRT18, DSCAM-AS1, SPAG6)

<img src="Heatmap.png" width="600" title="Bi-Clustered Heatmap">

### Volcano Plot

The volcano plot identified highly significant subtype-specific DEGs, highlighting aggressive Basal markers (e.g., ART3, CT83, TLX3 ,FOXC1, HORMAD1, CRABP1) and Luminal markers (e.g., GATA3, DSCAM-AS1, KRT18, SPAG6, DOK7)


<img src="Volcano plot.png" width="700" title="Volcano Plot Distribution">

### 3. Pathway Mapping (KEGG and GO Enrichment)

3. Pathway Mapping (KEGG and GO Enrichment)
* **Basal Pathways:** Shows specific functional enrichment for **Neuroactive ligand-receptor interactions**, structural cellular re-engineering tracks (**Keratinization**, **Keratinocyte differentiation**), and dense cytokeratin structural changes governing mechanical tissue invasion.
* **Luminal Pathways:** Driven primarily by the **PI3K-Akt signaling pathway**, intracellular macromolecule localization, and vesicle-mediated transport loops reflecting active metabolic synthesis over migration.

  <img src="KEGG_Pathways_Basal.png" width="45%" title="Basal Pathways">
  <img src="KEGG_Pathways_Luminal.png" width="45%" title="Luminal Pathways">
  <img src="GO_Pathways_Basal.png" width="45%" title="Basal GO Terms">
  <img src="GO_Pathways_Luminal.png" width="45%" title="Luminal GO Terms">


## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.


## How to Run

1.  **Clone the repository:**
    ```bash
    git clone <repository-url>
    cd Breast-Cancer-RNASeq-Pipeline
    ```

2.  **Install Dependencies:**
    -   SRA Toolkit (`fasterq-dump`)
    -   FastQC
    -   Trimmomatic
    -   HISAT2 & SAMtools
    -   Subread (`featureCounts`)
    -   R with packages: `DESeq2`, `tidyverse`, `pheatmap`, `EnhancedVolcano`, `clusterProfiler`

3.  **Download Data & References:**
    -   Download SRA files for the 10 samples listed in the `README.md`.
    -   Download the HISAT2 `grch38` genome index and the Ensembl `GRCh38` GTF annotation file.

4.  **Execute the Pipeline:**
    -   Modify and run the `pipeline.sh` script for each sample to generate a count matrix.
    -   Run the R scripts in the `/Scripts` directory to perform the downstream analysis.

    


*Created by [Jagriti Bezbaruah]*
    