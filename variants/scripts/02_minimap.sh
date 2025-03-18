#/bin/bash

set -e 

FOLDER="$1"
SAMPLE="$2"
REF="$3"

echo "FOLDER:"${FOLDER}
echo "SAMPLE:"${SAMPLE}
echo "REF":${REF}

#cp ${FOLDER}/filtlong/${SAMPLE}_filtlong.fastq.gz .
#cp ${FOLDER}/refs/*${REF}.fna .

ls .

minimap2 -x map-ont -a ${REF}_assembly.fasta ${SAMPLE}_filtlong.fastq.gz > ${SAMPLE}_vs_${REF}.sam

mkdir -p ${FOLDER}/mapping

echo "done"
