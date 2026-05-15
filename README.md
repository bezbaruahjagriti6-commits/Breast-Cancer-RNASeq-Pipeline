# Breast Cancer Subtype Analysis (Luminal vs. Basal)

**Author:** Jagriti Bezbaruah

##  Overview
This repository hosts a complete end-to-end bioinformatics workflow for the transcriptomic profiling of Breast Cancer subtypes. Using the GSE209998 dataset, the project identifies key molecular signatures that distinguish Luminal from Basal-like subtypes.


## Pipeline Architecture

RA (Raw Archive) ──> [prefetch / fasterq-dump] ──> FASTQ (Raw Reads) ──> [HISAT2] ──> SAM (Alignment Text) ──> [SAMtools view & sort] ──> BAM (Sorted Binary) ──> [featureCounts] ──> Raw Counts (Matrix) ──> [DESeq2 in R] ──> DGE Results (Biomarkers)

##  Pipeline & Analysis
The project is divided into two main stages:

### 1. The Processing Pipeline (`Scripts/pipeline.sh`)
Automates high-throughput sequence alignment:
* **Alignment:** HISAT2 (GRCh38 Reference).
* **Post-processing:** SAMtools for sorting and indexing.
* **Quantification:** featureCounts for gene-level expression matrices.

### 2. Statistical Analysis (`Scripts/dge_analysis.R`)
Performs advanced differential gene expression (DGE) modeling:
* **Normalization:** Median-of-ratios (DESeq2).
* **DGE Modeling:** Wald test ($p\text{-adj} < 0.05$, $|\text{log}_2\text{FC}| > 2$).
* **Enrichment:** KEGG Pathway analysis via clusterProfiler.

##  Results & Detailed Analysis
The analysis revealed a sharp transcriptomic divergence between Luminal and Basal samples.

**Key Biological Findings:**
* **Upregulated in Luminal:** Significant enrichment of hormone receptor-related pathways (e.g., $ESR1$, $PGR$).
* **Upregulated in Basal:** Identified a robust immune-cell infiltration signature, specifically within the $IGHV$ (Immunoglobulin Heavy Variable) gene cluster.
* **Pathway Enrichment:** Top pathways included Cell Cycle Control and DNA Replication, with Basal samples showing higher expression of proliferation markers ($MKI67$).

**Data Visualizations:**
* **Volcano Plot:** Highlights the most significantly shifted genes ($p$-value vs. Fold Change).
* **Clustered Heatmap:** Shows clear separation of sample groups based on the top 20 DEGs.
* **KEGG Dot Plot:** Visualizes the most impacted biological pathways.

*(See `z_Figures/` for the original plots).*

##  How to Run

1. **Clone the repository:**
   ```bash
   git clone https://github.com/bezbaruahjagriti6-commits/Breast-Cancer-RNASeq-Pipeline.git
   cd Breast-Cancer-RNASeq-Pipeline
   ```

2. **Add your raw data:**
   Place your raw `.fastq` files, the HISAT2 reference genome index, and your `.gtf` annotation file in your working directory. 
   *(Note: Large sequencing files like `.fastq`, `.sam`, and `.bam` are ignored by git to keep the repository lightweight).*

3. **Execute Shell Script:** 
   ```bash
   chmod +x Scripts/pipeline.sh
   ./Scripts/pipeline.sh
   ```
   *(Requires raw data).*

4. **Execute R Script:**
   ```bash
   Rscript Scripts/dge_analysis.R
   ```
   *(Generates plots and CSV results).*
