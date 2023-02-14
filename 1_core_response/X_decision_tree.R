# Select genes which are robustly changed through the treatments
library(stringr)
library(rpart)
library(rpart.plot)

diffexp = read.csv("output/differential_expression.lfs.csv.gz")

# Drop Yunlong columns and e6
diffexp = diffexp[,-grep("yunlong|e6", names(diffexp))]

# Build table for training: Classify each gene x sample in upregulated, downregulated, or non-significant (this is what we want to predict)
gene.effect = data.frame(sapply(unique(str_split_fixed(colnames(diffexp)[2:ncol(diffexp)], "_", 2)[,1]), FUN=function(s) {
  sapply(1:nrow(diffexp), FUN=function(g) {
    # @TODO not sure why there are NAs but there are (e.g. in the column "e5.2")
    if(is.na(diffexp[g,paste0(s,"_PValue")]) || diffexp[g,paste0(s,"_PValue")] >= 0.05)
      "ns"
    else {
      if(diffexp[g,paste0(s,"_logFC")] < 0)
        "down"
      else
        "up"
    }
  })
}))
row.names(gene.effect) = str_split_fixed(diffexp$gene, "\\.", 2)[,1]
gene.effect = data.frame(t(gene.effect))

# Add x vals for prediction
sample.annotation = unique(read.csv("input/samples_annotation.csv")[c("sample.group", "tissue", "temperature", "stress.duration")])
sample.annotation$sample.group = make.names(sample.annotation$sample.group)

# Merge them
gene.effect$sample.group = row.names(gene.effect)
train.table = merge(sample.annotation, gene.effect, by="sample.group", all=F)

# @TODO for Robert's data we will just assume 72 hours as duration
train.table[is.na(train.table$stress.duration),]$stress.duration = 72

train.table$tissue = as.factor(train.table$tissue)
#train.table$temperature = as.factor(train.table$temperature)

# Now we'll try predicting the effect on a specific gene by all our input factors
gene.predictability = data.frame(gene=character(), up=numeric(), down=numeric(), accuracy=numeric())
for(g in sample(names(gene.effect), 500)) {
  print(g)
  train.table[[g]] = factor(train.table[[g]], levels=c("down","ns","up"))
  r_tree =  rpart(get(g) ~ tissue + temperature + stress.duration, data = train.table, method="class")
  root.node.error = as.numeric(na.omit(str_match(capture.output(printcp(r_tree)), "Root node error: .* = (.*)$"))[1,2])
  cp = printcp(r_tree)
  accuracy = 1 - cp[nrow(cp),4] * root.node.error
  up = table(train.table[[g]])[3] / nrow(train.table)
  down = table(train.table[[g]])[1] / nrow(train.table)
  gene.predictability = rbind(gene.predictability, data.frame(gene=g, up=up, down=down, accuracy=accuracy))
}
gene.predictability$non.ns = gene.predictability$up + gene.predictability$down
write.csv(gene.predictability, "output/gene_predictability_decision_tree.lfs.csv", row.names=F)
system("pigz -11 output/gene_predictability_decision_tree.lfs.csv")

show_tree <- function(g) {
  train.table[[g]] = factor(train.table[[g]], levels=c("down","ns","up"))
  r_tree =  rpart(get(g) ~ tissue + temperature + stress.duration, data = train.table, method="class", control = rpart.control(cp = 0, maxdepth=3, minsplit=3))
  print(r_tree); printcp(r_tree); rpart.plot(r_tree)
}
# Example:
show_tree("Solyc02g005140")
# Example Solyc01g005985

for(g in a$gene) {
  tryCatch(
    {show_tree(g)},
    error=function(cond){}
  )
}


# Find an interesting tree
for(g in gene.predictability$gene) {
  train.table[[g]] = factor(train.table[[g]], levels=c("down","ns","up"))
  r_tree =  rpart(get(g) ~ tissue + temperature + stress.duration, data = train.table, method="class", control = rpart.control(cp = 0))
  if(length(capture.output(print(r_tree))) > 8)
    print(g)
}


# protein.descriptions = read.csv("input/protein_descriptions.csv.gz")
# responsive.genes$gene = str_split_fixed(responsive.genes$gene, "\\.", 2)[,1]
# responsive.genes = merge(responsive.genes, protein.descriptions)

