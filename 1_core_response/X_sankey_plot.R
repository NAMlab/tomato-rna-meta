library(networkD3)

samples.annotation = read.csv("input/samples_annotation.csv")
samples.annotation$sample.group = make.names(samples.annotation$sample.group)

# Only look at contrasts
samples.annotation = unique(samples.annotation[samples.annotation$stress_type != "none",c("publication_id", "tissue", "stress_type", "sample.group")])

links = data.frame(source = character(), target=character(), value=numeric())

a = as.data.frame(table(samples.annotation$publication_id, samples.annotation$stress_type))
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

p <- sankeyNetwork(Links = links, Nodes = nodes,
                   Source = "IDsource", Target = "IDtarget",
                   Value = "value", NodeID = "name", fontSize = 18, nodeWidth = 30,
                   sinksRight=FALSE)
p

