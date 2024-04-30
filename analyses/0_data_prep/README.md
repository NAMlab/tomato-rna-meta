## 0_data_prep
Data preparation for other analyses, currently only generates a file with descriptions and orthologs for all protein coding genes in the ITAG4.1 tomato genome.
The first script in this folder take all genes in the ITAG4.1 tomato genome and add Arabidopsis orthologes found on [OMA](https://omabrowser.org/) as well as the descriptions of these orthologues from [TAIR](https://www.arabidopsis.org/).
The second script then also adds any protein and gene descriptions from [UniProt](https://www.uniprot.org/) where it can find any.
None of the scripts perform a life search on the respective platforms, they use the downloaded files from the `data` folder.
