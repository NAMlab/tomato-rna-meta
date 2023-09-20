library(stringr)
library(htmlwidgets)
library(plotly)
library(ggfortify)
library(gridExtra)

source("../config.R")

samples.annotation = read.csv("input/samples_annotation.csv")
samples.annotation[samples.annotation$sra_run_id == "",]$sra_run_id = samples.annotation[samples.annotation$sra_run_id == "",]$sample.name

abundance = read.csv("input/combined_abundance.tsv.gz", sep="\t", check.names=F)
row.names(abundance) = abundance$target_id

plot_and_save_pca <- function(pca_res, point.annotations, outfile_name, x = 1, y = 2) {
  scores = as.data.frame(pca_res$x)
  
  ## 2D
  pdf(paste0("output/plots/", outfile_name, ".pdf"))
  print(autoplot(pca_res, data = point.annotations, colour='tissue', shape='stress_type', x = x, y = y) +
          scale_color_manual(values = colors[["tissue"]]) + theme_minimal())
  dev.off()
}

### ----- log(TPMs) -----
log.tpms = log(abundance[grep("_tpm", names(abundance))])
# Kick out rows which contain -Inf (formerly 0 values)
log.tpms = log.tpms[rowSums(log.tpms == -Inf) == 0,]
log.tpms = t(as.matrix(log.tpms))
row.names(log.tpms) = str_remove(row.names(log.tpms), "_tpm")

pca_res <- prcomp(log.tpms)
scores = as.data.frame(pca_res$x)

point.annotations = data.frame(sample = row.names(scores))
point.annotations = merge(point.annotations, samples.annotation, by.x = "sample", by.y="sra_run_id", all.x=T)

plot.log.tpm <- autoplot(pca_res, data = point.annotations, colour='tissue', shape='stress_type', x = 1, y = 2) +
  scale_color_manual(values = colors[["tissue"]]) + ggtitle("PCA of logTPMs") + theme_minimal()+ theme(legend.position = "none")

pdf("output/plots/1.2_supplement_PCA_all_logTPMs.lfs.pdf", 14, 8.5)
for(i in seq(3, 337, 2)) {
  message(i)
  print(autoplot(pca_res, data = point.annotations, colour='tissue', shape='stress_type', x = i, y = i+1) +
          scale_color_manual(values = colors[["tissue"]]) + ggtitle("PCA of logTPMs") + theme_minimal())
}
dev.off()

library(ape)
pdf("output/plots/1.2_supplement_clustering_all_logTPMs.pdf", 15, 60)
sa = samples.annotation[order(samples.annotation$sra_run_id),]
sa$color <- colors[["tissue"]][match(sa$tissue, names(colors[["tissue"]]))]
tree = hclust(d = dist(x = log.tpms, method = "euclidean"))
tree.phylo = as.phylo(tree)
plot(tree.phylo)
dev.off()

rm(abundance, scores, pca_res, log.tpms, point.annotations)

#### ---- fold-changes -------
diffexp = read.csv("output/differential_expression.lfs.csv.gz")
row.names(diffexp) = diffexp$gene

fcs = t(as.matrix(diffexp[grep("_logFC", names(diffexp))]))
row.names(fcs) = str_remove(row.names(fcs), "_logFC")

pca_res <- prcomp(fcs)
scores = as.data.frame(pca_res$x)

point.annotations = data.frame(sample = row.names(scores))
samples.annotation$sample.group = str_replace_all(samples.annotation$sample.group, "-", ".")
samples.annotation = unique(samples.annotation[c("genotype.name", "tissue", "stress.duration", "temperature", "sample.group", "stress_type")])
point.annotations = merge(point.annotations, samples.annotation, by.x = "sample", by.y="sample.group", all.x=T, all.y=F)

plot.fcs <- autoplot(pca_res, data = point.annotations, colour='tissue', shape='stress_type', x = 1, y = 2) +
  scale_color_manual(values = colors[["tissue"]]) + ggtitle("PCA of log Fold Changes") + theme_minimal()

pdf("output/plots/1.2_PCA.pdf", 14, 8.5)
print(grid.arrange(plot.log.tpm, plot.fcs, ncol = 2))
dev.off()


## Finally, OPLS-DA
# Finally, we want to see if we can split the contrasts according to stress types, tissues etc.
# This could be very interesting to see which genes are used here.
library(ropls)

fcs.heat = fcs[which(point.annotations$stress_type == "heat"),]
points.heat = point.annotations[point.annotations$stress_type == "heat",]

heat.psda = opls(fcs.heat, points.heat$tissue, predI=1, subset="odd")
trainVi <- getSubsetVi(heat.psda)
confusion_train.tb <- table(points.heat$tissue[trainVi], fitted(heat.psda))
confusion_train.tb

confusion_test.tb <- table(points.heat$tissue[-trainVi],
                           predict(heat.psda, fcs.heat[-trainVi,]))
confusion_test.tb

# For a combination of multiple y's (online possible with numeric values):
heat.psda <- opls(fcs.heat, as.matrix(points.heat[c("stress.duration", "temperature")]))

#### (snippet for playing with UMAP) ###
# umap_res = umap(tpms)
# row.names(point.annotations) = point.annotations$sample
# p = merge(point.annotations, umap_res$layout, by="row.names")
# ggplot(data=p, aes(x=V1, y=V2))+
#        geom_point(aes(color=publication_id, shape=stress_type)) + xlab("UMAP 1") + ylab("UMAP 2")
