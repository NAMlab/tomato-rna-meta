library(chromoMap)

a = read.csv("input/core_genes.csv", header=F)
q = read.csv("input/qtls.csv", header=F)
# Transform MB positions to b
q$V3 = q$V3 * 1000000
q$V4 = q$V4 * 1000000
a = rbind(a, q)
write.table(a, "input/out.txt", sep="\t", row.names=F, col.names = F, quote=F)

chromoMap("input/chromosomes.txt", "input/out.txt", labels=T, segment_annotation = T,
          data_based_color_map = T, interactivity = T, 
          data_type = "categorical", export.options = T,
          data_colors = list(c("red","lightblue", "green")))
