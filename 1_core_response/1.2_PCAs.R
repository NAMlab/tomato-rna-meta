library(stringr)
library(htmlwidgets)
library(plotly)
library(ggfortify)
library(gridExtra)
library(VCA)
library(ape)

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
  scale_color_manual(name="Tissue", values = colors[["tissue"]]) + ggtitle("PCA of logTPMs") + 
  scale_shape_manual(name="Treatment", values=c("drought"=15,"heat"=16,"salt"=17,"none"=3), drop=F) + 
  theme_minimal(base_size = 18) + theme(legend.position = "left") + guides(color = guide_legend(override.aes = list(size = 3)), shape = guide_legend(override.aes = list(size = 3)))
plot.log.tpm

pdf("output/plots/1.2_supplement_PCA_all_logTPMs.lfs.pdf", 14, 8.5)
for(i in seq(3, 337, 2)) {
  message(i)
  print(autoplot(pca_res, data = point.annotations, colour='tissue', shape='stress_type', x = i, y = i+1) +
          scale_color_manual(values = colors[["tissue"]]) + ggtitle("PCA of logTPMs") + theme_minimal())
}
dev.off()

# Now do VCA
vca.res = data.frame(gene = character(), percent.stress_type = numeric(), percent.tissue = numeric(),
                     percent.datasource_id = numeric(), percent.genotype = numeric(), percent.error = numeric())
samples.annotation$tissue = as.factor(samples.annotation$tissue)
samples.annotation$datasource_id = as.factor(samples.annotation$datasource_id)
samples.annotation$stress_type = as.factor(samples.annotation$stress_type)
samples.annotation$genotype.name = as.factor(samples.annotation$genotype.name)
for(i in 1:ncol(log.tpms)) {
  print(i)
  df.vca = cbind(data.frame(log.tpm = log.tpms[,i]), samples.annotation)
  tryCatch({
    a = fitVCA(log.tpm ~ stress_type + tissue + datasource_id + genotype.name, df.vca, method="reml")
    vca.res = rbind(vca.res, list(gene = colnames(log.tpms)[i], 
                                  percent.stress_type = a$aov.tab[2,3],
                                  percent.tissue = a$aov.tab[3,3],
                                  percent.datasource_id = a$aov.tab[4,3],
                                  percent.genotype = a$aov.tab[5,3],
                                  percent.error = a$aov.tab[6,3]
    ))
  }, error=function(cond) {
    
  })
}
write.csv(vca.res, "output/1.2_VCA_logTPMs.csv", row.names=F)
system("pigz -11 output/1.2_VCA_logTPMs.csv")

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
point.annotations$row_number = 1:nrow(point.annotations)
samples.annotation$sample.group = str_replace_all(samples.annotation$sample.group, "-", ".")
samples.annotation = unique(samples.annotation[c("genotype.name", "tissue", "stress.duration", "temperature", "sample.group", "stress_type", "datasource_id")])
point.annotations = merge(point.annotations, samples.annotation, by.x = "sample", by.y="sample.group", all.x=T, all.y=F)
point.annotations = point.annotations[order(point.annotations$row_number),]

plot.fcs <- autoplot(pca_res, data = point.annotations, colour='tissue', shape='stress_type', x = 1, y = 2) +
  scale_color_manual(values = colors[["tissue"]]) + 
  scale_shape_manual(values=c("drought"=15,"heat"=16,"salt"=17,"none"=3), drop=F) + 
  theme_minimal(base_size = 16) + theme(legend.position = "none") +
  ggtitle("PCA of log Fold Changes")

pdf("output/plots/1.2_PCA.pdf", 14, 8.5)
plot_grid(plot.log.tpm, plot.fcs, align="h", ncol=2, rel_widths = c(0.58, 0.42))
dev.off()

# Now do VCA
point.annotations$stress_type = droplevels(point.annotations$stress_type)
vca.res = data.frame(gene = character(), percent.stress_type = numeric(), percent.tissue = numeric(),
                     percent.datasource_id = numeric(), percent.genotype = numeric(), percent.error = numeric())
for(i in 1:ncol(fcs)) {
  print(i)
  df.vca = cbind(data.frame(log.fc = fcs[,i]), point.annotations)
  tryCatch({
    a = fitVCA(log.fc ~ stress_type + tissue + datasource_id + genotype.name, df.vca, method="reml")
    vca.res = rbind(vca.res, list(gene = colnames(fcs)[i], 
                                  percent.stress_type = a$aov.tab[2,3],
                                  percent.tissue = a$aov.tab[3,3],
                                  percent.datasource_id = a$aov.tab[4,3],
                                  percent.genotype = a$aov.tab[5,3],
                                  percent.error = a$aov.tab[6,3]
    ))
  }, error=function(cond) {
    
  })
}
write.csv(vca.res, "output/1.2_VCA_logFC.csv", row.names=F)
system("pigz -11 output/1.2_VCA_logFC.csv")

## Finally, OPLS-DA
# Finally, we want to see if we can split the contrasts according to stress types, tissues etc.
# This could be very interesting to see which genes are used here.
# library(ropls)
# 
# fcs.heat = fcs[which(point.annotations$stress_type == "heat"),]
# points.heat = point.annotations[point.annotations$stress_type == "heat",]
# 
# heat.psda = opls(fcs.heat, points.heat$tissue, predI=1, subset="odd")
# trainVi <- getSubsetVi(heat.psda)
# confusion_train.tb <- table(points.heat$tissue[trainVi], fitted(heat.psda))
# confusion_train.tb
# 
# confusion_test.tb <- table(points.heat$tissue[-trainVi],
#                            predict(heat.psda, fcs.heat[-trainVi,]))
# confusion_test.tb
# 
# # For a combination of multiple y's (online possible with numeric values):
# heat.psda <- opls(fcs.heat, as.matrix(points.heat[c("stress.duration", "temperature")]))
