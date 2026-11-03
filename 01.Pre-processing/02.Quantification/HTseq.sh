#!/bin/bash
# Usage: bash 01.HTseq.sh <SampleID> settings.sh
ID=$1; source $2

make_dir ${HTSEQ_DIR}

TMP="${HTSEQ_DIR}/${ID}.namesorted.bam"
samtools sort -n -@ 4 -o ${TMP} ${ALIGN_ALL}/${ID}.bam

htseq-count -f bam -r name -s no -a 10 -t exon -m intersection-strict \
    ${TMP} ${GTF} > ${HTSEQ_DIR}/${ID}.htseq.count.txt

rm -f ${TMP}
