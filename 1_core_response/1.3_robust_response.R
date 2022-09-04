# Select genes which are robustly changed through the treatments
library(stringr)

diffexp = read.csv("output/differential_expression.lfs.csv.gz")

responsive.genes = data.frame(gene=str_split_fixed(diffexp$gene, "\\.", 2)[,1])
responsive.genes$median=apply(diffexp[grep("_PValue", names(diffexp))], 1, median, na.rm=T)
responsive.genes$median.adj = p.adjust(responsive.genes$median, method="BH")
responsive.genes$sigCount = apply(diffexp[grep("_FDR", names(diffexp))], 1, FUN=function(x) { sum(x < 0.05) })

protein.descriptions = read.csv("input/protein_descriptions.csv.gz")
responsive.genes$gene = str_split_fixed(responsive.genes$gene, "\\.", 2)[,1]
responsive.genes = merge(responsive.genes, protein.descriptions)

