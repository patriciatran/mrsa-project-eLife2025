#!/bin/bash

SAMPLE="$1"

ls -lh

tar -xvzf bakta_${SAMPLE}_1.10.3.tar.gz

ls -lh 
ls -lh ${SAMPLE}_annotation/*

mkdir ${SAMPLE}_out
ls

export EGGNOG_DATA_DIR="/projects/bacteriology_tran_data/eggnog"
echo $EGGNOG_DATA_DIR

emapper.py -i ${SAMPLE}_annotation/${SAMPLE}_assembly.faa \
	--itype proteins \
	--data_dir $EGGNOG_DATA_DIR \
	--temp_dir ./tmp \
	--scratch_dir ./scratch \
	--cpu 20 \
	--output ${SAMPLE}_ \
	--excel \
	--report_orthologs \
	--decorate_gff yes \
	--pfam_realign realign \
	--output_dir ${SAMPLE}_out

ls -lh

tar -czvf ${SAMPLE}_eggnog.tar.gz ${SAMPLE}_out*

ls -lh


