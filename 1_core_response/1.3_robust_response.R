# Select genes which are robustly changed through the treatments
library(stringr)

diffexp = read.csv("output/differential_expression.lfs.csv.gz")
diffexp$median = apply(diffexp[grep("_PValue", names(diffexp))], 1, median, na.rm=T)
diffexp$median = p.adjust(diffexp$median, method="BY")
diffexp = diffexp[diffexp$median < 0.05,]


protein.descriptions = read.csv("input/protein_descriptions.csv.gz")
diffexp$gene = str_split_fixed(diffexp$gene, "\\.", 2)[,1]
diffexp = merge(diffexp, protein.descriptions)
