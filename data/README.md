## data
This directory contains all data needed to reproduce the analyses in this repository. 

| filename                            | description          | source                                                                 |
|-------------------------------------|----------------------|------------------------------------------------------------------------|
| GOMAP_annotation.map.lfs.gz         | GO terms for all genes in the tomato genome, reformated for the R package `topGO`.  | https://doi.org/10.25739/zh2v-4p15 |
| ITAG4.1_PrDs.csv                    | predicted prion-like protein domains | [PLAAC](https://github.com/whitehead/plaac ) on the ITAG4.1 proteome with an alpha value of 0.5. |
| ITAG4.1_descriptions.lfs.txt.gz     | names of the proteins in ITAG4.1 | https://solgenomics.net/ftp/tomato_genome/annotation/ITAG4.1_release/ITAG4.1_descriptions.txt |
| ITAG4.1_proteins.lfs.fasta.gz       | amino acid sequences for proteins in ITAG4.1| https://solgenomics.net/ftp/tomato_genome/annotation/ITAG4.1_release/ITAG4.1_proteins.fasta |
| TAIR10_functional_descriptions.lfs.tsv.gz | description of Arabidopsis genes from TAIR 10 | https://www.arabidopsis.org/download_files/Genes/TAIR10_genome_release/TAIR10_functional_descriptions |
| combined_abundance.lfs.tsv.gz       | kallisto mapping results for all samples (counts and TPMs) | created with https://github.com/NAMlab/rnaseq-mapper |
| display_datasource_mapping.csv      | mapping from original publication_id in the samples list to the datasource id (DS-X) for the manuscript | manually created |
| fimo_hsfs_matches.tsv               | matches of a heat shock factor consensus sequence in the tomato genome | see below |
| id_mapping.csv.gz                   | mapping of Solyc identifiers to NCBI, UniProt and UniParc IDs | https://ensembl.gramene.org/biomart/martview , Plant Genes 68 |
| oma_orthologues.tsv.gz              | Arabidopsis orthologs for each tomato gene | https://omabrowser.org/cgi-bin/gateway.pl?f=PairwiseOrthologs&p1=SOLLC&p2=ARATH&p3=Source |
| samples_annotation.csv              | samples used in this analysis and their metadata | manually collected and curated |
| uniprot_descriptions.tsv.gz         | protein and gene names for UniProt accessions | UniProt |


### fimo_hsfs_matches.tsv
FIMO from meme version 4.12. was run with default parameters using the motif below on the 500bp promoters of all protein coding tomato genes (extracted using https://github.com/RimGubaev/extract_promoters ).
Background frequencies (for bgfile) were estimated by running fasta-get-markov on the chromosome sequences of SL4.0.

```
MEME version 4.4

ALPHABET= ACGT

strands: + -

Background letter frequencies (from /mnt/thumper-e1/home/shhuang/projects/dap/analysis.v4/gem07_rep_memechip03/ABI3VP1_tnt/AT5G18090_col_a/background):
A 0.30000 C 0.20000 G 0.20000 T 0.30000 

MOTIF HSF_manual.tnt.HSFC1_col_a_m1 HSF_manual

letter-probability matrix: alength= 4 w= 8 nsites= 600 E= 6.6e-1430
  0.000000	  0.005000	  0.995000	  0.000000	
  0.980000	  0.003333	  0.000000	  0.016667	
  0.896667	  0.001667	  0.030000	  0.071667	
  0.205000	  0.163333	  0.401667	  0.230000	
  0.116667	  0.605000	  0.211667	  0.066667	
  0.060000	  0.001667	  0.000000	  0.938333	
  0.010000	  0.000000	  0.000000	  0.990000	
  0.000000	  1.000000	  0.000000	  0.000000
  ```