library(rjson)
edges = data.frame(gene = character(), family = character())
nodes.genes = data.frame(gene = character(), direction = character(), heat = logical(), drought = logical(), salt = logical())

for(stress.type in c("heat", "drought", "salt")) {
  for(direction in c("upregulated", "downregulated")) {
    f = fromJSON(file=paste0("input/", stress.type, "_", direction, ".json"))
    # Nodes
    # add new nodes to the table
    if(length(setdiff(names(f), nodes.genes$gene)) > 0) {
      nodes.genes = rbind(nodes.genes, data.frame(gene = setdiff(names(f), nodes.genes$gene),
                                                  direction = direction, heat=F, drought=F, salt=F))
    }
    # set stress type to T for all nodes in the respective JSON
    nodes.genes[nodes.genes$gene %in% names(f),][[stress.type]] = T
    # Edges
    for(g in names(f)) {
      for(e in f[[g]]) {
        edges = rbind(edges, data.frame(gene=g, family=e[[3]]))
      }
    }
  }
}
edges = unique(edges)

nodes.genes$color = rgb(nodes.genes$heat*0.5, nodes.genes$salt*0.5, nodes.genes$drought*0.5)

sink('network.dot')
cat("digraph G {\n")
cat(paste0('  "', nodes.genes$gene, '"[fillcolor="',nodes.genes$color,'",fontcolor="white",style="filled",label="',nodes.genes$gene,'"];', collapse="\n"))
cat(paste0('  "', edges$gene, '" -> "', edges$family, '";', collapse="\n"))
cat("\n}")
sink()
system("circo -Tpdf network.dot > network.pdf")

# Now let's prune it a bit
# Remove families with only 1 gene
a = table(edges$family)
edges = edges[edges$family %in% names(a[a>1]),]
# Remove genes without edges
nodes.genes = nodes.genes[nodes.genes$gene %in% unique(edges$gene),]

sink('network.dot')
cat("digraph G {\n")
cat(paste0('  "', nodes.genes$gene, '"[fillcolor="',nodes.genes$color,'",fontcolor="white",style="filled",label="',nodes.genes$gene,'"];', collapse="\n"))
cat(paste0('  "', edges$gene, '" -> "', edges$family, '";', collapse="\n"))
cat("\n}")
sink()
system("circo -Tpdf network.dot > network_pruned.pdf")
