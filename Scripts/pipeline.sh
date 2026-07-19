
#!/bin/bash
# Project: Breast Cancer Subtype Analysis (Basal vs Luminal)

SECONDS=0

echo "Starting Bioinformatics Pipeline..."

# Convert .sra to .fastq
fasterq-dump --split-files --progress

# Run fastqc
fastqc data/demo.fastq -o data/

# Run trimmomatic 
java -jar ~/Desktop/demo/tools/Trimmomatic-0.39/trimmomatic-0.39.jar SE -threads 4 data/demo.fastq data/demo_trimmed.fastq LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:36
echo "Trimmomatic finished running!"

# Run HISAT2
# mkdir HISAT2
# get the genome indices
# wget https://genome-idx.s3.amazonaws.com/hisat/grch38_genome.tar.gz

# Run alignment and output sorted BAM directly
hisat2 -q --rna-strandness R -x HISAT2/grch38/genome -U data/demo_trimmed.fastq | samtools sort -o HISAT2/demo_sorted.bam
echo "HISAT2 finished running!"

# Run featureCounts - Quantification
# get gtf
# wget http://ftp.ensembl.org/pub/release-106/gtf/homo_sapiens/Homo_sapiens.GRCh38.106.gtf.gz

#featureCounts command 
featureCounts -T 4 -t exon -g gene_id -a ../annotation/Homo_sapiens.GRCh38.106.gtf -o quants/demo_featurecounts.txt HISAT2/demo_sorted.bam
echo "featureCounts finished running!"

echo "Pipeline script completed."

duration=$SECONDS
echo "$(($duration / 60)) minutes and $(($duration % 60)) seconds elapsed."





