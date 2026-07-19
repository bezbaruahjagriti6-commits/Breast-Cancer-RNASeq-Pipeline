library(DESeq2)
library(tidyverse)
library(EnhancedVolcano)
library(pheatmap)
library(clusterProfiler)
library(enrichplot)
library(org.Hs.eg.db)
library(RColorBrewer)
library(biomaRt)
library(ggplot2)

#load data
counts_data <- read.csv("Breastcancer_counts.csv", row.names = 1)
coldata <- read.csv("coldata.csv", row.names = 1)

#Run DESeq2
dds <- DESeqDataSetFromMatrix(countData = counts_data,
                              colData = coldata, 
                              design = ~ Subtype)
dds <- DESeq(dds)
res <- results(dds, contrast=c("Subtype", "Basal", "LuminalA"))
res_df<- as.data.frame(res)

#PCA plot
rld <- rlog(dds)
pca_plot <- plotPCA(rld, intgroup = "Subtype") +
  theme_bw() 
  labs(title = "PCA Plot:",
       subtitle = "Breast Cancer Subtypes")

print(pca_plot)

#MA plot
plotMA(res, alpha = 0.05, main = "MA Plot: Basal vs Luminal A")

#convert ensembl id into gene names
mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

gene_map <- getBM(attributes = c("ensembl_gene_id", "external_gene_name"),
                  filters = "ensembl_gene_id",
                  values = rownames(res_df),
                  mart = mart)

res_annotated <- res_df%>%
  rownames_to_column(var = "ensembl_gene_id") %>%
  left_join(gene_map, by = "ensembl_gene_id")

res_clean <- res_annotated %>% 
  filter(!is.na(external_gene_name) & external_gene_name != "" & !is.na(padj))

#significant degs
sig_degs<- res_clean%>% 
  filter(padj < 0.05& abs(log2FoldChange) >= 1)


upregulated_genes <- sig_degs %>% filter(log2FoldChange > 1) %>% arrange(desc(log2FoldChange))
downregulated_genes <- sig_degs %>% filter(log2FoldChange < 1) %>% arrange(log2FoldChange)


#VOLCANO PLOT

EnhancedVolcano(
  res_clean,
  lab = res_clean$external_gene_name,
  x = 'log2FoldChange',
  y = 'padj',
  pCutoff = 0.05,                     
  FCcutoff = 2.0,                     
  pointSize = 1.5,
  labSize = 3.0,
  labCol = 'black',
  labFace = 'bold',
  col = c('grey', 'green', 'blue', 'red'),
  colAlpha = 0.6,
  title = 'DEG: TNBC vs. Luminal A',
  subtitle = "Volcano plot",
  caption = paste0('Total analyzed features = ', nrow(res_clean), ' | Thresholds: Log2FC > 2.0, padj < 0.05'),
  legendPosition = 'right',
  legendLabSize = 10,
  legendIconSize = 4.0
)
#HEATMAP 
plot_matrix <- assay(rld)
rownames(plot_matrix) <- res_annotated$external_gene_name[match(rownames(plot_matrix), res_annotated$ensembl_gene_id)]

custom_labels <- c(
  "CT83", "ART3", "TLX3", "KRT16", "FOXC1","WNT16", "ITGA10","MSLN","SLC26A9","SOX10","NLRP2","HORMAD1","PRSS41","CRABP1","LINC02437",
  "SPAG6", "TFF3", "FOXA1", "DSCAM-AS1", "PGR", "GATA3", "CCND1", "AR", "KRT18", "KRT8", "BCL2", "AGR2","DOK7","ABCC12"
)

heatmap_genes <- custom_labels[custom_labels %in% rownames(plot_matrix)]
plot_matrix_subset <- plot_matrix[heatmap_genes, ]

sample_info <- data.frame(Subtype = coldata$Subtype)
rownames(sample_info) <- colnames(plot_matrix_subset)
sample_info$Subtype <- as.factor(sample_info$Subtype)

ann_colors <- list(Subtype = c(Basal = "coral2", LuminalA = "darkcyan"))

pheatmap(
  plot_matrix_subset,
  scale = "row",            
  clustering_distance_rows = "correlation",
  cluster_cols = FALSE,           
  cluster_rows = TRUE,            
  annotation_col = sample_info,    
  annotation_colors = ann_colors,  
  # Color gradient defined using standard text names instead of hex numbers:
  color = colorRampPalette(c("navyblue", "royalblue3", "lightskyblue", "lightyellow", "orange", "red2", "darkred"))(100),
  border_color = "white",          
  fontsize_row = 10,          
  fontsize_col = 9,                
  main = "HEATMAP: Basal vs Luminal A Breast Cancer"
  )
  
#KEGG Enrichment

sig_kegg <- sig_degs[sig_degs$padj < 0.05 & !is.na(sig_degs$padj),]
 
#map gene symbols to ENTERZID ID
entrez_mapped <- suppressWarnings(
  bitr(
    sig_kegg$external_gene_name,
    fromType = "SYMBOL",
    toType   = "ENTREZID",
    OrgDb    = org.Hs.eg.db
  )
)

# Merge back to keep the Log2FoldChange data tied to the new Entrez IDs
merged_kegg <- merge(entrez_mapped, sig_kegg, by.x = "SYMBOL", by.y = "external_gene_name")

# Separate into distinct cohort vectors based on the Log2 Fold Change scale
basal_entrez   <- merged_kegg$ENTREZID[merged_kegg$log2FoldChange > 0]
luminalA_entrez <- merged_kegg$ENTREZID[merged_kegg$log2FoldChange < 0]


# 3. RUN KEGG PATHWAY ENRICHMENT

#  Basal Pathway Analysis 
kk_basal <- enrichKEGG(
  gene          = basal_entrez,
  organism      = "hsa",               
  pvalueCutoff  = 0.05,
  pAdjustMethod = "BH"                 
)

kk_basal <- setReadable(kk_basal, OrgDb = org.Hs.eg.db, keyType="ENTREZID")

# LuminalA Pathway Analysis 
kk_luminalA <- enrichKEGG(
  gene          = luminal_entrez,
  organism      = "hsa",
  pvalueCutoff  = 0.05,
  pAdjustMethod = "BH"
)
kk_luminalA <- setReadable(kk_luminalA, OrgDb = org.Hs.eg.db, keyType="ENTREZID")


# Generate publication-grade dotplots showing the top 10 enriched pathways
dotplot(kk_basal, showCategory = 10, title = "Top KEGG Pathways in Basal")
dotplot(kk_luminalA, showCategory = 10, title = "Top KEGG Pathways in Luminal A")


#GO Enrichment Analyis

sig_df <- sig_degs[sig_degs$padj < 0.05 & !is.na(sig_degs$padj), ]
             
# Split into explicit lineage vectors based on the Log2 Fold Change scale
basal_genes <- sig_df$external_gene_name[sig_degs$log2FoldChange > 1]
luminalA_genes <- sig_df$external_gene_name[sig_degs$log2FoldChange < 1]

# Define the entire dataset profile as your reference genomic baseline
background_genes <- res_df$external_gene_name

print(paste("Basal Genes:", length(basal_genes), "| LuminalA Genes:", length(luminalA_genes)))


# 2. RUN GENE ONTOLOGY (GO) BIOLOGICAL PROCESS ANALYSIS

#  Basal Enrichment Pipeline 
go_basal <- enrichGO(
  gene          = basal_genes,
  universe      = background_genes,
  OrgDb         = org.Hs.eg.db,
  keyType       = "SYMBOL",
  ont           = "BP",       
  pAdjustMethod = "BH",           
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  readable      = TRUE            
)

# Luminal A Enrichment Pipeline 
go_luminalA <- enrichGO(
  gene          = luminal_genes,
  universe      = background_genes,
  OrgDb         = org.Hs.eg.db,
  keyType       = "SYMBOL",
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  readable      = TRUE
)


# Generate horizontal dot plots (ideal for checking gene counts vs. p-adjustments)
dotplot(go_basal, showCategory = 15, title = "Top Biological Processes in Basal")
dotplot(go_luminalA, showCategory = 15, title = "Top Biological Processes in Luminal A")


