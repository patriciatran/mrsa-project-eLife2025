#!/bin/bash
# This scripts takes trimmed fastq illumina data and removes the host (Staph aureus), uses the non-host sequences for assembly.

#SAMPLE="EDN_SA_Ph1_S203"
#HOST="SaureusRN4220.fna"
#STAGING="/staging/ptran5"
#CPUS=8


SAMPLE="$1"
CPUS="$2"
STAGING="$3"
HOST="$4"

cp ${STAGING}/hosts/${HOST}* .

minimap2 --split-prefix=tmp$$ -a -xsr ${HOST} ${SAMPLE}.1P.fastq.gz ${SAMPLE}.2P.fastq.gz | samtools view -bh | samtools sort -o output.bam
samtools index output.bam

mkdir host not_host
echo "Host:"
samtools fastq -F 3588 -f 65 output.bam | gzip -c > host/output_${SAMPLE}_R1.fastq.gz
echo "R2 matching host genome:"
samtools fastq -F 3588 -f 129 output.bam | gzip -c > host/output_${SAMPLE}_R2.fastq.gz

echo "Not Host:"
samtools fastq -F 3584 -f 77 output.bam  | gzip -c > not_host/output_${SAMPLE}_R1.fastq.gz
samtools fastq -F 3584 -f 141 output.bam | gzip -c > not_host/output_${SAMPLE}_R2.fastq.gz
samtools fastq -f 4 -F 1 output.bam | gzip -c > not_host/output_S_Singletons.fastq.gz

echo "Using non-host for SPADES assembly:"
spades.py -1 not_host/output_${SAMPLE}_R1.fastq.gz -2 not_host/output_${SAMPLE}_R2.fastq.gz -o ${SAMPLE}_output_folder --only-assembler -t ${CPUS}

tar -zcvf ${SAMPLE}_output_folder.tar.gz ${SAMPLE}_output_folder
