#/bin/bash

SAMPLE="$1"
FOLDER="$2"

ls -lht

filtlong -p 95 ${SAMPLE}.fastq.gz | gzip > ${SAMPLE}_filtlong.fastq.gz

ls -lh

mkdir -p ${FOLDER}/filtlong

echo "done"
