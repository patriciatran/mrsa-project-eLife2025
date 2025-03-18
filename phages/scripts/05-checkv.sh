#!/bin/bash

SAMPLE="$1"
CPU="$2"

#CHECKVDB=/projects/bacteriology_tran_data/checkV_db

checkv end_to_end ${SAMPLE}.fasta ${SAMPLE}_checkv_out -t ${CPU} -d /projects/bacteriology_tran_data/checkV_db/checkv-db-v1.5/

tar -czvf ${SAMPLE}_checkv_out.tar.gz ${SAMPLE}_checkv_out
