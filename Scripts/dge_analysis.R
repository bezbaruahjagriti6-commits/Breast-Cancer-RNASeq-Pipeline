# Step 1: Install and Load
if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install("DESeq2")
library(DESeq2)

# Step 2: Load your file (Skip the first line of metadata)
data <- read.table("final_counts.txt", header=TRUE, skip=1, row.names=1)

# Step 3: Pick your two count columns
# Assuming columns 6 and 7 are your actual counts
counts <- data[, c(6, 7)] 
colnames(counts) <- c("Luminal", "Basal")

# Step 4: Create a simple metadata table
info <- data.frame(condition=factor(c("Luminal", "Basal")))

# Step 5: The "Hack" for No Replicates
dds <- DESeqDataSetFromMatrix(countData = counts, colData = info, design = ~condition)
dds <- estimateSizeFactors(dds)

# We manually tell R the "noise" level is 0.1 so it doesn't crash
dispersions(dds) <- 0.1 
dds <- nbinomWaldTest(dds)

# Step 6: Get your results!
res <- results(dds)
res_sorted <- res[order(res$log2FoldChange, decreasing=TRUE), ]
write.csv(as.data.frame(res_sorted), "Final_Project_Results.csv")
# Install the translation tools
if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("AnnotationDbi", "org.Hs.eg.db"))

library(AnnotationDbi)
library(org.Hs.eg.db)

# Convert IDs (This assumes your row names are ENSG IDs)
res_df <- as.data.frame(res_sorted)
res_df$symbol <- mapIds(org.Hs.eg.db,
                        keys=rownames(res_df),
                        column="SYMBOL",
                        keytype="ENSEMBL",
                        multiVals="first")

# Save the version with names!
write.csv(res_df, "Final_Results_with_Names.csv")

# Basic Volcano Plot
plot(res_df$log2FoldChange, -log10(res_df$pvalue), 
     pch=20, main="Volcano Plot: Basal vs Luminal", 
     xlab="log2 Fold Change", ylab="-log10 P-value",
     col=ifelse(abs(res_df$log2FoldChange) > 2, "red", "black"))

# Add horizontal/vertical threshold lines
abline(h=-log10(0.05), col="blue", lty=2) # Significance threshold
abline(v=c(-2, 2), col="blue", lty=2)    # Fold-change threshold

# Sort by log2FoldChange to find the most significant movers
top_up <- head(res_df[order(res_df$log2FoldChange, decreasing = TRUE), ], 5)
top_down <- head(res_df[order(res_df$log2FoldChange, decreasing = FALSE), ], 5)

print(top_up[, c("symbol", "log2FoldChange", "pvalue")])
print(top_down[, c("symbol", "log2FoldChange", "pvalue")])


if (!require("pheatmap", quietly = TRUE)) install.packages("pheatmap")
library(pheatmap)

# Select the top 20 genes (10 up, 10 down)
top_genes <- rbind(head(res_df[order(res_df$log2FoldChange, decreasing = TRUE), ], 10),
                   head(res_df[order(res_df$log2FoldChange, decreasing = FALSE), ], 10))

# Get the normalized counts for these genes
# Note: counts(dds, normalized=TRUE) gives normalized counts
mat <- counts(dds, normalized=TRUE)[rownames(top_genes), ]

# Use symbols as row names for the plot
rownames(mat) <- ifelse(is.na(top_genes$symbol), rownames(top_genes), top_genes$symbol)

# Log transform for better colors
mat_log <- log2(mat + 1)

# Draw the Heatmap
pheatmap(mat_log, 
         cluster_rows = TRUE, 
         cluster_cols = TRUE, 
         show_colnames = TRUE,
         main = "Top Differentially Expressed Genes")


# 1. Install the necessary packages
if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install("clusterProfiler")
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)

# 2. Convert your Gene Symbols to Entrez IDs
# Assuming 'res_df' is your results dataframe
genes_to_test <- res_df$symbol[!is.na(res_df$symbol)]

entrez_ids <- bitr(genes_to_test, 
                   fromType = "SYMBOL", 
                   toType = "ENTREZID", 
                   OrgDb = org.Hs.eg.db)

# 3. Run KEGG Enrichment
kegg_results <- enrichKEGG(gene = entrez_ids$ENTREZID,
                           organism = 'hsa', # 'hsa' is for Human
                           pvalueCutoff = 0.05)

# 4. Visualize the results
dotplot(kegg_results, showCategory=15) + 
  ggtitle("KEGG Pathway Enrichment")

