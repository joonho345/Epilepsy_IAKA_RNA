#!/bin/bash
# settings.sh — edit the two lines below, then source before running shell scripts
#   source settings.sh && bash 01.Pre-processing/01.fastp_paired.sh <SampleID> settings.sh

DATA_DIR="/path/to/data"       # pipeline data output root (FASTQ, BAM, HTSeq)
REF_DIR="/path/to/reference"   # GRCm39 reference files

GTF="${REF_DIR}/Mus_musculus.GRCm39.112.gtf"
REF="${REF_DIR}/Mus_musculus.GRCm39.dna.primary_assembly.fa"
INDEX="${DATA_DIR}/index/Mouse_101"

FASTQ_DIR="${DATA_DIR}/01.raw"
FASTP_DIR="${DATA_DIR}/03.fastp/${ID}"
ALIGN_DIR="${DATA_DIR}/04.aligned/${ID}"
ALIGN_ALL="${DATA_DIR}/04.aligned/ALL"
HTSEQ_DIR="${DATA_DIR}/02.HTseq"

FASTQ1="${ID}_1.fastq.gz"
FASTQ2="${ID}_2.fastq.gz"
FASTP1="${ID}.FP_R1.fq.gz"
FASTP2="${ID}.FP_R2.fq.gz"

make_dir() { [ -d "$1" ] || mkdir -p "$1"; }
