# Visualize our data set as a sankey diagram

library(networkD3)
library(htmlwidgets)
sankey.keys = read.csv("input/display_datasource_mapping.csv")
samples.annotation = read.csv("input/samples_annotation.csv")
samples.annotation$sample.group = make.names(samples.annotation$sample.group)
samples.annotation = merge(samples.annotation, sankey.keys, by.x="publication_id", by.y="old_id")

# Only look at contrasts
samples.annotation = unique(samples.annotation[samples.annotation$stress_type != "none",c("display_id", "tissue", "stress_type", "sample.group")])

links = data.frame(source = character(), target=character(), value=numeric())

a = as.data.frame(table(samples.annotation$display_id, samples.annotation$stress_type))
names(a) = c("source", "target", "value")
links = rbind(links, a)

a = as.data.frame(table(samples.annotation$stress_type, samples.annotation$tissue))
names(a) = c("source", "target", "value")
links = rbind(links, a)

links = links[links$value > 0,]

nodes <- data.frame(
  name=unique(c(as.character(links$source), 
         as.character(links$target)))
)

links$IDsource <- match(links$source, nodes$name)-1 
links$IDtarget <- match(links$target, nodes$name)-1


node.names = sankey.keys$display_id
node.cols = rep("black", length(node.names))
node.names = c(node.names, "drought", "salt", "heat", "anther", "fruit", "leaf", "pollen", "seed", "seedling", "ovaries", "root")
node.cols = c(node.cols, "#0072B2", "#009E73", "#D55E00", "#CC79A7", "#D55E00", "#009E73", "#F0E442", "#E69F00", "#56B4E9", "pink", "brown")

my_color <- paste0('d3.scaleOrdinal() .domain([', '"', paste(node.names, collapse='", "'), '"]) .range([', '"', paste(node.cols, collapse='", "'), '"])')
p <- sankeyNetwork(Links = links, Nodes = nodes,
                   Source = "IDsource", Target = "IDtarget",
                   Value = "value", NodeID = "name", fontSize = 18, nodeWidth = 30,
                   sinksRight=FALSE, colourScale = my_color, width=600, height=800)
p
saveWidget(p, "output/plots/1_sankey_plot.html")
