library(UpSetR)
library(data.table)
library(stringr)
diffexp = read.csv("output/differential_expression.lfs.csv.gz")

samples.annotation = read.csv("input/samples_annotation.csv")
samples.annotation[samples.annotation$sra_run_id == "",]$sra_run_id = samples.annotation[samples.annotation$sra_run_id == "",]$sample.name
samples.annotation$sample.group = make.names(samples.annotation$sample.group)
samples.annotation = samples.annotation[samples.annotation$sample.group %in% unique(str_split_fixed(names(diffexp), "_", 2)[,1]),]

degs = list()
for(contrast in unique(str_split_fixed(names(diffexp), "_", 2)[,1])[-1]) {
  print(contrast)
  if (any(is.na(diffexp[[paste0(contrast, "_FDR")]]))) {
    next
  }
  degs[[paste0(contrast, "_up")]] <- diffexp[diffexp[[paste0(contrast, "_FDR")]] < 0.05 & diffexp[[paste0(contrast, "_logFC")]] > 0,]$gene
  degs[[paste0(contrast, "_down")]] <- diffexp[diffexp[[paste0(contrast, "_FDR")]] < 0.05 & diffexp[[paste0(contrast, "_logFC")]] < 0,]$gene
}

metadata = unique(samples.annotation[c("sample.group", "tissue", "genotype.name", "temperature", "stress_type", "stress.duration")])
metadata$stress.duration = log(metadata$stress.duration)
names(metadata)[1] = "sets"
names(metadata)[3] = "genotype"
metadata = metadata[metadata$sets %in% unique(str_split_fixed(names(degs), "_", 2)[,1]),]
metadata_up = metadata
metadata_up$sets = paste0(metadata_up$sets, "_up")
metadata_up$direction = "up"
metadata_down = metadata
metadata_down$sets = paste0(metadata_down$sets, "_down")
metadata_down$direction = "down"
metadata = rbind(metadata_up, metadata_down)
set.metadata = list(data = metadata, plots = list(
    list(type = "heat", column = "stress_type", assign = 3,
         colors = c(drought="#0072B2", heat="#009E73", salt="#D55E00")),
    list(type = "heat", colors=c("anther" = "#CC79A7", "fruit" = "#D55E00", "leaf" = "#009E73", "pollen" = "#F0E442", "seed" = "#E69F00", "seedling" = "#56B4E9", "ovaries"="pink", "root"="brown"),
         column = "tissue", assign=3),
    list(type = "text", column = "genotype", assign = 8),
    list(type = "heat", column = "stress.duration", assign = 3),
    list(type="bool", column="direction", assign=3)
))

pdf("../../plots/1.X_big_upset.pdf", 16, 20)
print(upset(fromList(degs), nsets=120, order.by="freq", nintersects = 80, set.metadata = set.metadata,
            mb.ratio = c(0.3, 0.7)))
dev.off()
