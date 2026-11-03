#!/bin/bash
# Usage: bash 02.STAR_index.sh settings.sh
# sjdbOverhang = read_length - 1 = 100 for 101 bp reads
source $1; make_dir ${INDEX}

STAR \
    --runThreadN 10 \
    --runMode genomeGenerate \
    --genomeDir ${INDEX} \
    --genomeFastaFiles ${REF} \
    --sjdbGTFfile ${GTF} \
    --sjdbOverhang 100
