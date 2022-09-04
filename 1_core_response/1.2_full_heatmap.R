library(ComplexHeatmap)
library(stringr)

diffexp = read.csv("output/differential_expression.lfs.csv.gz")

# Filter out genes which are not differentially expressed in any sample
min_fdr = apply(diffexp[grep("_FDR", names(diffexp))], 1, FUN = min) 
diffexp = diffexp[min_fdr < 0.01,]

logfc = as.matrix(diffexp[grep("_logFC", names(diffexp))])
rownames(logfc) = diffexp$gene
colnames(logfc) = str_split_fixed(colnames(logfc), "_", 2)[,1]

png("output/plots/full_heatmap.png")
heatmap(logfc)
dev.off()