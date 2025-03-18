#!/bin/bash
SAMPLE="$1"
CPUS="$2"

tar -xvzf ${SAMPLE}_spades_output_folder.tar.gz

blastn -db nr \
	-query ${SAMPLE}_output_folder/scaffolds.fasta \
	-remote \
	-max_target_seqs 5 \
	-out ${SAMPLE}.result.blast.tsv \
	-evalue 1e-10 \
	-outfmt 6
#	-num_threads ${CPUS}


