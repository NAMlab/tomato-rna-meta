library(rjson)
library(grid)
library(circlize)
library(stringr)
edges = data.frame(gene = character(), family = character(), ipr.id = character())
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
        edges = rbind(edges, data.frame(gene=g, family=e[[3]], ipr.id=e[[1]]))
      }
    }
  }
}
edges = unique(edges)

nodes.families = unique(edges[2:3])
nodes.families = merge(nodes.families, read.csv("output/family_sizes.csv"), by.x="ipr.id", by.y="interpro.id")

nodes.genes$color = rgb(nodes.genes$heat*0.5, nodes.genes$salt*0.5, nodes.genes$drought*0.5)
protein.descriptions = read.csv("../0_data_prep/output/protein_descriptions.lfs.csv.gz")[c(1,5)]
protein.descriptions = protein.descriptions[!duplicated(protein.descriptions[,c('gene')]),]
protein.descriptions$ITAG4.1_description = str_remove(protein.descriptions$ITAG4.1_description, " \\(AHRD V.*$")
nodes.genes$gene.short = str_remove(nodes.genes$gene, "\\..*")
nodes.genes = merge(nodes.genes, protein.descriptions, by.x="gene.short", by.y="gene", all.x=T)
nodes.genes$shape = ifelse(nodes.genes$direction == "upregulated", 'oval', 'box')

sink('network.dot')
cat("digraph G {\n")
cat(paste0('  "', nodes.genes$gene, '"[fillcolor="',nodes.genes$color,'",fontcolor="white",style="filled",shape="',nodes.genes$shape,'",label="',paste0(nodes.genes$ITAG4.1_description, '\n(', nodes.genes$gene),')"];', collapse="\n"))
cat(paste0('  "', nodes.families$ipr.id, '"[label="',paste0(nodes.families$family, '\n(', nodes.families$n.tomato.proteins),' members in tomato) "];', collapse="\n"))
cat(paste0('  "', edges$gene, '" -> "', edges$ipr.id, '";', collapse="\n"))
cat("\n}")
sink()
system("circo -Tpdf network.dot > output/plots/3.2_network.pdf")

# Now let's prune it a bit
# Remove families with only 1 gene
a = table(edges$family)
edges = edges[edges$family %in% names(a[a>1]),]
# Remove genes without edges
nodes.genes = nodes.genes[nodes.genes$gene %in% unique(edges$gene),]
nodes.families = nodes.families[nodes.families$ipr.id %in% unique(edges$ipr.id),]

sink('network.dot')
cat("digraph G {\n")
cat(paste0('  "', nodes.genes$gene, '"[fillcolor="',nodes.genes$color,'",fontcolor="white",style="filled",shape="',nodes.genes$shape,'",label="',paste0(nodes.genes$ITAG4.1_description, '\n(', nodes.genes$gene),')"];', collapse="\n"))
cat(paste0('  "', nodes.families$ipr.id, '"[label="',paste0(nodes.families$family, '\n(', nodes.families$n.tomato.proteins),' members in tomato) "];', collapse="\n"))
cat(paste0('  "', edges$gene, '" -> "', edges$ipr.id, '";', collapse="\n"))
cat("\n}")
sink()
system("circo -Tpdf network.dot > output/plots/3.2_network_pruned.pdf")

# Try visualizing it in a "kind of heatmap" situation
library(ComplexHeatmap)

pdf("output/plots/3.2_heatmap.pdf", 17, 30)
nodes.genes = nodes.genes[order(nodes.genes$gene),]
row_ha = rowAnnotation(heat = nodes.genes$heat, drought = nodes.genes$drought, salt = nodes.genes$salt, direction=nodes.genes$direction,
                       col = list(heat = c("TRUE" = "#880000", "FALSE" = "white"), drought = c("TRUE" = "#000088", "FALSE" = "white"), salt = c("TRUE" = "#008800", "FALSE" = "white"),
                                  direction = c("upregulated" = "red", "downregulated" = "blue")))
draw(Heatmap(table(edges$gene, edges$family),
             right_annotation = row_ha,
             #column_names_gp = grid::gpar(fontsize = 10),
             column_names_max_height = max_text_width(edges$family),
             col=c("grey90", "black")
))
dev.off()
