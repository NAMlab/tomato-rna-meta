# This script takes the produced protein descriptions and adds UniProtKB descriptions to it.
library(stringr)

id_mapping = read.csv("input/id_mapping.csv.gz")
uniprot_descriptions = read.csv("input/uniprot_descriptions.tsv.gz", sep="\t")
desc = read.csv("output/protein_descriptions.lfs.csv.gz")

id_mapping$Gramene.ID = str_split_fixed(id_mapping$Gramene.ID, "\\.", 2)[,1]
m = merge(id_mapping, uniprot_descriptions, by.x="UniProtKB.ID", by.y="Entry")[c("Gramene.ID","Reviewed","Protein.names", "Gene.Names")]
names(m) = c("gene", "UniProtKB.reviewed", "UniProtKB.protein.names", "UniProtKB.gene.names")
d = merge(desc, m, all.x=T)
write.csv(d, "output/protein_descriptions.lfs.csv", row.names = F, na="")
system("pigz -11 output/protein_descriptions.lfs.csv")
