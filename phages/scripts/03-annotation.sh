#!/bin/bash

cp /staging/ptran5/MoLab/cleaned_data/$1_spades_output_folder.tar.gz
tar -xvzf $1_spades_output_folder.tar.gz

bakta --db /projects/bacteriology_tran_data/bakta/db/ \
	$1_output_folder/scaffolds.fasta \
	--output bakta_$1 \
	--threads $2

tar -czvf bakta_$1.tar.gz bakta_$1
