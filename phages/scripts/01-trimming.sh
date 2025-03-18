#!/bin/bash
# This script trims the raw FAST reads.

# This takes raw FASTQ files, and trims using Trimmomatic.

SAMPLE="$1"
CPU="$2"
STAGING="$3"

# Illumina adapter to use for trimming:
wget https://github.com/Microbial-Ecology-Group/AMRplusplus/raw/refs/heads/master/data/adapters/nextera.fa

#echo "Step 1: trim sequences"

trimmomatic \
      PE \
      -threads ${CPU} \
      ${SAMPLE}_R1_001.fastq.gz ${SAMPLE}_R2_001.fastq.gz \
      ${SAMPLE}.1P.fastq.gz ${SAMPLE}.1U.fastq.gz \
      ${SAMPLE}.2P.fastq.gz ${SAMPLE}.2U.fastq.gz \
      ILLUMINACLIP:nextera.fa:2:30:10:3:TRUE \
      LEADING:3 \
      TRAILING:3 \
      SLIDINGWINDOW:"4:15" \
      MINLEN:36