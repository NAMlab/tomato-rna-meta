library(stringr)

id.mapping = read.csv("input/id_mapping.csv.gz")
uniprot.descriptions = read.csv("input/uniprot_descriptions.tsv.gz", sep="\t")
desc = read.csv("output/protein_descriptions.lfs.csv.gz")

id.mapping$Gramene.ID = str_split_fixed(id.mapping$Gramene.ID, "\\.", 2)[,1]
m = merge(id.mapping, uniprot.descriptions, by.x="UniProtKB.ID", by.y="Entry")[c("Gramene.ID","Reviewed","Protein.names", "Gene.Names")]
names(m) = c("gene", "UniProtKB.reviewed", "UniProtKB.protein.names", "UniProtKB.gene.names")
d = merge(desc, m, all.x=T)
write.csv(d, "output/protein_descriptions.lfs.csv", row.names = F, na="")
system("pigz -11 output/protein_descriptions.lfs.csv")
