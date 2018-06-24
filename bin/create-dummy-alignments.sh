#!/bin/bash

aln_dir="../alignments"

if [ ! -e "$aln_dir" ]
then
    mkdir -p "$aln_dir"
fi

NSPECIES=2
NGENOMES=40
NCHARS=500

i=1
while [ "$i" -lt 6 ]
do
    prefix="c${i}sp"
    outfile="${aln_dir}/comp0${i}-${NSPECIES}species-${NGENOMES}genomes-00${NCHARS}chars.nex"
    ./generate-dummy-biallelic-alignment.py \
        --nspecies "$NSPECIES" \
        --ngenomes "$NGENOMES" \
        --ncharacters "$NCHARS" \
        --prefix "$prefix" \
        > "$outfile"
    i=`expr $i + 1`
done

NCHARS=1000

i=1
while [ "$i" -lt 6 ]
do
    prefix="c${i}sp"
    outfile="${aln_dir}/comp0${i}-${NSPECIES}species-${NGENOMES}genomes-0${NCHARS}chars.nex"
    ./generate-dummy-biallelic-alignment.py \
        --nspecies "$NSPECIES" \
        --ngenomes "$NGENOMES" \
        --ncharacters "$NCHARS" \
        --prefix "$prefix" \
        > "$outfile"
    i=`expr $i + 1`
done

NCHARS=10000

i=1
while [ "$i" -lt 6 ]
do
    prefix="c${i}sp"
    outfile="${aln_dir}/comp0${i}-${NSPECIES}species-${NGENOMES}genomes-${NCHARS}chars.nex"
    ./generate-dummy-biallelic-alignment.py \
        --nspecies "$NSPECIES" \
        --ngenomes "$NGENOMES" \
        --ncharacters "$NCHARS" \
        --prefix "$prefix" \
        > "$outfile"
    i=`expr $i + 1`
done
