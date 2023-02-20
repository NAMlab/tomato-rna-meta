library(stringr)

d = read.csv("output/combined_deepest_levels.csv")
d = d[d$set == "heat",]
d$level = d$level + 1

labels = scan("clades_hierarchy.txt", character())[(min(d$level)):(max(d$level))]
labels[length(labels)] = "S. lycopersicum"

labels = data.frame(level = 1:length(labels), origin=labels)
d = merge(d, labels, all.x=T)

protein.descriptions = read.csv("../0_data_prep/output/protein_descriptions.lfs.csv.gz")[c(1,5)]
protein.descriptions = protein.descriptions[!duplicated(protein.descriptions[,c('gene')]),]
protein.descriptions$ITAG4.1_description = str_remove(protein.descriptions$ITAG4.1_description, " \\(AHRD V.*$")
d$gene.short = str_remove(d$target, "\\..*")
d = merge(d, protein.descriptions, by.x="gene.short", by.y="gene", all.x=T)

d = d[order(d$level),c("target", "origin", "ITAG4.1_description")]
write.csv(d, "output/origins_of_heat_genes.csv", row.names=F)
