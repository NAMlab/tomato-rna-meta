#!/bin/bash
# This script saves FASTA extracts of the ITAG protein sequences of the core response genes.
# It requires `seqtk` to be installed.

for stress_type in `tail -n+2 "output/1.4_core_response_genes.csv" | cut -f1 -d"," | sort | uniq`; do
  for direction in `tail -n+2 "output/1.4_core_response_genes.csv" | cut -f2 -d"," | sort | uniq`; do
    echo $stress_type,$direction
    grep $stress_type,$direction "output/1.4_core_response_genes.csv" | cut -f3 -d"," | sort > id.lst
    seqtk subseq ../data/ITAG4.1_proteins.lfs.fasta.gz id.lst | awk '{print $1}' | sed 's/\*//g' > output/sequences/${stress_type}_${direction}.lfs.fasta
  done
done

# Get a random 500 sequences for baseline
#zgrep ">" ../data/ITAG4.1_proteins.lfs.fasta.gz | cut -f1 -d" " | cut -c2- | shuf | head -500 | sort > random.lst
seqtk subseq ../data/ITAG4.1_proteins.lfs.fasta.gz id.lst | awk '{print $1}' | sed 's/\*//g' > output/sequences/random_500.lfs.fasta

rm id.lst random.lst

