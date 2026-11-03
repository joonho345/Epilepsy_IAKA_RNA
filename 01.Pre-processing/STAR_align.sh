#!/bin/bash
# Usage: bash 03.STAR_align.sh <SampleID> settings.sh
ID=$1; source $2

make_dir ${ALIGN_DIR}
make_dir ${ALIGN_ALL}

STAR \
    --runThreadN 10 \
    --runMode alignReads \
    --genomeDir ${INDEX} \
    --readFilesIn ${FASTP_DIR}/${FASTP1} ${FASTP_DIR}/${FASTP2} \
    --sjdbOverhang 100 \
    --sjdbGTFfile ${GTF} \
    --readFilesCommand gunzip -c \
    --genomeLoad NoSharedMemory \
    --twopassMode Basic \
    --outFileNamePrefix ${ALIGN_DIR}/${ID}_ \
    --outSAMtype BAM SortedByCoordinate \
    --outSAMattributes All

samtools index ${ALIGN_DIR}/${ID}_Aligned.sortedByCoord.out.bam

ln -sf ${ALIGN_DIR}/${ID}_Aligned.sortedByCoord.out.bam     ${ALIGN_ALL}/${ID}.bam
ln -sf ${ALIGN_DIR}/${ID}_Aligned.sortedByCoord.out.bam.bai ${ALIGN_ALL}/${ID}.bam.bai
