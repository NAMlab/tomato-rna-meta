#!/bin/bash

# Save a random 500 protein IDs from the genomes in random.lst
grep ">" ITAG4.1_proteins.fasta | cut -f1 -d" " | cut -c2- | shuf | head -500 | sort > random.lst

# Get those sequences
seqtk subseq ITAG4.1_proteins.fasta random.lst | grep -v ">" > seqs

# Save them together
echo "target,peptide_sequence" > random.csv
paste -d"," random.lst seqs >> random.csv

# Now do the same for the HS core genes (this list needs to be sorted!)
seqtk subseq ITAG4.1_proteins.fasta hs_core_genes.lst | grep -v ">" > seqs
echo "target,peptide_sequence" > hs_core_genes.csv
paste -d"," hs_core_genes.lst seqs >> hs_core_genes.csv

# Clean up
rm seqs random.lst
