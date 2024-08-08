#!/bin/bash
# For each gene in our heat stress core set and the random baseline set:
# align the sequences (both CDS and protein) using mafft and calculate
# position-specific conservation scores using al2co ( https://doi.org/10.1093/bioinformatics/17.8.700 )
# For the protein sequences we're using the BLOSUM80 substitution matrix, for the 
# CDS I spoofed one based on that (see NUC file)
# Saves alignments to output/alignments and conservation scores to output/al2co_scores

module load seqkit mafft

grep ">" "analyses/1_core_response/output/sequences/random_500.lfs.fasta" | cut -c2- > random_baseline_ids.txt
grep ">" "analyses/1_core_response/output/sequences/heat_upregulated.fasta" | cut -c2- > heat_core_ids.txt
grep ">" "analyses/1_core_response/output/sequences/heat_downregulated.fasta" | cut -c2- >> heat_core_ids.txt

cat hs_transcript_ids.txt random_baseline_ids.txt | shuf > all_ids.txt
rm random_baseline_ids.txt heat_core_ids.txt

mkdir -p output/proteins/alignments output/proteins/al2co_scores output/cds/alignments output/cds/al2co_scores

while read t; do
  echo "Protein $t"
  echo "  extracting sequences"
  bash 2_get_seqs.sh "$t" "proteins" > sequences.fa
  echo "  running alignment"
  mafft --auto --clustalout --reorder "sequences.fa" > "output/proteins/alignments/$t.aln" 2>/dev/null
  # Remove header line
  tail -n +2 output/proteins/alignments/$t.aln > alignment.aln
  echo "  calculating conservation scores"
  al2co/al2co -i alignment.aln -n F -c 2 -s input/BLOSUM80 -o scores.tsv
  # Remove last 12 lines
  head -n -14 scores.tsv > output/proteins/al2co_scores/$t.txt
done <all_ids.txt

while read t; do
  echo "CDS $t"
  echo "  extracting sequences"
  bash 2_get_seqs.sh "$t" "cds" > sequences.fa
  echo "  running alignment"
  mafft --auto --clustalout --reorder "sequences.fa" > "output/cds/alignments/$t.aln" 2>/dev/null
  # Remove header line
  tail -n +2 output/cds/alignments/$t.aln > alignment.aln
  echo "  calculating conservation scores"
  al2co/al2co -i alignment.aln -n F -c 2 -s input/NUC -o scores.tsv
  # Remove last 12 lines
  head -n -14 scores.tsv > output/cds/al2co_scores/$t.txt
done <all_ids.txt

rm all_ids.txt scores.tsv alignment.aln sequences.fa
