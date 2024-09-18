# Make another separate file only containing the entries of the hs core genes for convenience.
cut -f1 -d"," ../analyses/1_core_response/input/hs_core_genes_internal_names.csv > hs_core_ids.txt
zgrep -f hs_core_ids.txt output/gene_family_annotations.csv.gz | sed 's/,/ /' > output/hs_core_gene_family_annotations.csv
rm hs_core_ids.txt
