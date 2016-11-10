## Mapping of the Kc167 RNA-seq reads.

# Load the needed packages
setwd('/RnaSeq')
library(Rsubread)
library(limma)
library(edgeR)
library(GenomicFeatures)

# Index the dmel Genome
buildindex(basename="dmel_5.53",reference="/data/dmel-all-chromosome-r5.53.fasta")

# Map the reads
align(index="dmel_5.53",readfile1='/data/Marc_RNA_Kc167_7407_allreads.fastq.gz',input_format="gzFASTQ",output_format="BAM", output_file="./Marc_RNA_Kc167_7407_all_mapped.bam",unique=TRUE,indels=4)

## Not finished. Not a priority
