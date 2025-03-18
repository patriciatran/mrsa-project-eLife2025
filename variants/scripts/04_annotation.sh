#!/bin/bash
SAMPLE="$1"
CPU="$2"

echo "SAMPLE : ${SAMPLE}"
echo "CPUS: ${CPU}"

MPLCONFIGDIR="./tmp"

bakta --db /projects/bacteriology_tran_data/bakta/db/ -v --tmp-dir ./tmp --keep-contig-headers --output ${SAMPLE}_annotation --threads ${CPU} ${SAMPLE}_assembly.fasta

tar -czvf bakta_${SAMPLE}.tar.gz ${SAMPLE}_annotation
