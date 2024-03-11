## data
This directory contains all data needed to reproduce the analyses in this repository. 

| filename                            | description          | source                                                                 |
|-------------------------------------|----------------------|------------------------------------------------------------------------|
| GOMAP_annotation.map.lfs.gz         | GO terms for all genes in the tomato genome, reformated for the R package `topGO`.  | https://doi.org/10.25739/zh2v-4p15 |
| ITAG4.1_PrDs.csv                    |                      |                                                                        |
| ITAG4.1_descriptions.lfs.txt.gz     |                      |                                                                        |
| ITAG4.1_proteins.lfs.fasta.gz       |                      |                                                                        |
| TAIR10_functional_descriptions.lfs.tsv.gz |                      |                                                                        |
| combined_abundance.lfs.tsv.gz       |                      |                                                                        |
| display_datasource_mapping.csv      |                      |                                                                        |
| fimo_hsfs_matches.tsv               |                      |                                                                        |
| id_mapping.csv.gz                   |                      |                                                                        |
| oma_orthologues.tsv.gz              |                      |                                                                        |
| samples_annotation.csv              |                      |                                                                        |
| uniprot_descriptions.tsv.gz         |                      |                                                                        |


- https://solgenomics.net/ftp//tomato_genome/annotation/ITAG4.1_release/ITAG4.1_descriptions.txt
- https://omabrowser.org/cgi-bin/gateway.pl?f=PairwiseOrthologs&p1=SOLLC&p2=ARATH&p3=Source
- https://www.arabidopsis.org/download_files/Genes/TAIR10_genome_release/TAIR10_functional_descriptions

PrDs were predicted using PLAAC (https://github.com/whitehead/plaac ) on the ITAG4.1 proteome with an alpha value of 0.5.
