#!/bin/bash
# Usage: bash 01.fastp_paired.sh <SampleID> settings.sh
ID=$1; source $2

make_dir ${FASTP_DIR}

fastp \
    --in1  ${FASTQ_DIR}/${FASTQ1} \
    --in2  ${FASTQ_DIR}/${FASTQ2} \
    --out1 ${FASTP_DIR}/${FASTP1} \
    --out2 ${FASTP_DIR}/${FASTP2} \
    --html ${FASTP_DIR}/${ID}.report.html \
    --json ${FASTP_DIR}/${ID}.report.json \
    --length_required 50 \
    --low_complexity_filter \
    --detect_adapter_for_pe \
    --overrepresentation_analysis -P 20 -w 4
