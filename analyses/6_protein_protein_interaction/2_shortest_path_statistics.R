paths = read.csv("output/shortest_paths.csv")
paths$length = paths$length - 1 # remove the starting node
print("Overall difference in path lengths:")
print(t.test(paths$length ~ paths$set))

print("Looking only at those genes outside of QLTs:")
qtl.locs = read.csv("../1_core_response/input/qtls.csv", header=F)
genes.loc = read.csv("../1_core_response/input/core_genes.csv", header=F)
# Transform MB positions to b
qtl.locs$V3 = qtl.locs$V3 * 1000000
qtl.locs$V4 = qtl.locs$V4 * 1000000

inside_qtls = c()
for(i in 1:nrow(genes.loc)) {
  x = genes.loc[i,]
  if(any(qtl.locs$V2 == c(x[2]) & qtl.locs$V3 <= c(x[3]) & qtl.locs$V4 >= c(x[4])))
    inside_qtls = c(inside_qtls, unlist(x[1]))
}

names = read.csv("../1_core_response/input/hs_core_genes_internal_names.csv")
inside_qtl_ids = names[names$internal.name %in% inside_qtls,]$target

print("Genes inside QTLs:")
print(t.test(paths[paths$source %in% inside_qtl_ids,]$length, paths[paths$set == "random",]$length))

print("Genes outside QTLs:")
print(t.test(paths[!paths$source %in% inside_qtl_ids & paths$set == "heat-core",]$length, paths[paths$set == "random",]$length))
