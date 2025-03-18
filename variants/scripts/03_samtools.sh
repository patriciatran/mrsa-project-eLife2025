#!/bin/bash

FOLDER="$1"
SAMPLE="$2"
REF="$3"

echo "convert sam to bam for $FILE"
samtools view -Sb ${SAMPLE}_vs_${REF}.sam > ${SAMPLE}_vs_${REF}.bam

echo "filter bam files for mapped reads only"
samtools view -b -F 4 ${SAMPLE}_vs_${REF}.bam > ${SAMPLE}_vs_${REF}.mapped.bam

echo "sort BAM, mapped"
samtools sort -o ${SAMPLE}_vs_${REF}.mapped.sorted.bam ${SAMPLE}_vs_${REF}.mapped.bam 

echo "sort BAM"
samtools sort -o ${SAMPLE}_vs_${REF}.sorted.bam ${SAMPLE}_vs_${REF}.bam 

echo "index sorted BAM FILE"
# note order: input then output:
samtools index ${SAMPLE}_vs_${REF}.sorted.bam ${SAMPLE}_vs_${REF}.sorted.bai

echo "create coverage tsv file"
samtools depth ${SAMPLE}_vs_${REF}.sorted.bam > ${SAMPLE}_vs_${REF}.coverage.tsv

echo "done"

