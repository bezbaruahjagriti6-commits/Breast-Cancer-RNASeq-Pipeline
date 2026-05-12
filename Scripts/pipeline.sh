#!/bin/bash
# Project: Breast Cancer Subtype Analysis (Luminal vs Basal)
# Author: Jagriti Bezbaruah

echo "Starting Bioinformatics Pipeline..."

# 1. Alignment with HISAT2
# hisat2 -x grch38_index -1 basal_R1.fastq -2 basal_R2.fastq -S basal.sam

# 2. Convert to BAM and Sort
# samtools view -bS basal.sam | samtools sort -o basal_sorted.bam

# 3. Quantification with featureCounts
# featureCounts -a annotation.gtf -o counts.txt basal_sorted.bam luminal_sorted.bam

echo "Pipeline script completed."
