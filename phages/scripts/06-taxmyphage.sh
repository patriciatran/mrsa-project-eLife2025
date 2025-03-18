#!/bin/bash

SAMPLE="$1"
CPU="$2"
DB="$3"

#mkdir genomes
#cp /staging/ptran5/MoLab/charlie/phage_genomes/* genomes/.

#ls genomes/*

#taxmyphage run -i ./genomes/ -t ${CPU} -o taxmyphage_${SAMPLE} -db ${DB}
taxmyphage run -i ${SAMPLE}.fasta -t ${CPU} -o taxmyphage_${SAMPLE} -db ${DB}

tar -czvf taxmyphage_${SAMPLE}.tar.gz taxmyphage_${SAMPLE}

