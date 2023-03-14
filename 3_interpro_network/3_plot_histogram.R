library(rjson)

edges = data.frame(gene = character(), family = character(), ipr.id = character())
nodes.genes = data.frame(gene = character(), direction = character(), heat = logical(), drought = logical(), salt = logical())

fams = data.frame(gene = character(), direction = character(), stress.type = character(), family = character())
for(stress.type in c("heat", "drought", "salt")) {
  for(direction in c("upregulated", "downregulated")) {
    f = fromJSON(file=paste0("input/", stress.type, "_", direction, ".json"))
    # Edges
    for(g in names(f)) {
      for(e in f[[g]]) {
        fams = rbind(fams, data.frame(gene=g, family=e[[3]], stress.type=stress.type, direction=direction))
      }
    }
  }
}

# Filter out families with only one member
a = table(fams$family)
fams = fams[fams$family %in% names(a[a > 2]), ]

a = table(fams$stress.type, fams$family)
col.order = names(sort(colSums(a), decreasing = T))
a = table(fams$stress.type, fams$family, fams$direction)
# order by frequency
a = a[, col.order,]
pdf("output/plots/3.3_families_barplot.pdf", 16, 16)
par(mar=c(8, 2, 1, 1), mfrow=c(2,1))
barplot(a[,,"upregulated"], col=c("#0000FF", "#FF0000", "#00FF00"), las=2,
        legend.text = rownames(a), # Legend values
        args.legend = list(x = "topright", inset = c(.1, 0)))
par(mar=c(0,2,1,1))
barplot(-a[,,"downregulated"], col=c("#0000FF", "#FF0000", "#00FF00"), names.arg = "")
dev.off()
