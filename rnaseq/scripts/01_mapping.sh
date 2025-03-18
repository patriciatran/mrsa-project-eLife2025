#/bin/bash

SAMPLE="$1"
REF="$2"

ls -lht

bowtie2-build ${REF}_assembly.fasta ${REF}
ls -lht

bowtie2 -x ${REF} -1 ${SAMPLE}_R1_001.fastq.gz -2 ${SAMPLE}_R2_001.fastq.gz -S ${SAMPLE}_vs_${REF}.sam --very-sensitive -p $3

ls -lht
samtools view -bS ${SAMPLE}_vs_${REF}.sam | samtools sort -o ${SAMPLE}_vs_${REF}.sorted.bam

ls -lht

