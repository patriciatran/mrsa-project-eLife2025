#/bin/bash

SAMPLE="$1"
REF="$2"

featureCounts -T 8 -t gene \
    -g gene \
    -a ${REF}_assembly.gtf \
    -o {SAMPLE}_vs_${REF}_counts.txt \
    ${SAMPLE}_vs_${REF}.sorted.bam
