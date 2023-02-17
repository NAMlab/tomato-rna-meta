library(stringr)
library(htmlwidgets)
library(plotly)
library(ggfortify)

samples.annotation = read.csv("input/samples_annotation.csv")
samples.annotation[samples.annotation$sra_run_id == "",]$sra_run_id = samples.annotation[samples.annotation$sra_run_id == "",]$sample.name

abundance = read.csv("input/combined_abundance.tsv.gz", sep="\t", check.names=F)
row.names(abundance) = abundance$target_id

plot_and_save_pca <- function(pca_res, point.annotations, outfile_name) {
  scores = as.data.frame(pca_res$x)
  
  ## 2D
  pdf(paste0("output/plots/", outfile_name, ".pdf"))
  print(autoplot(pca_res, data = point.annotations, colour='tissue', shape='stress_type'))
  dev.off()
  
  ## 3D
  fig <- plot_ly(scores, x = ~PC1, y = ~PC2, z = ~PC3,
                 text=paste0(point.annotations$sample, "\n",
                             "Sample Group: ", point.annotations$sample.group, "\n",
                             "Temperature: ", point.annotations$temperature, "\n",
                             "Duration: ", point.annotations$stress.duration, "\n",
                             "Genotype: ", point.annotations$genotype.name, "\n"
                 ), 
                 hoverinfo="text", mode="markers", type="scatter3d", color = as.factor(point.annotations$tissue)) %>%
    layout(scene = list(xaxis = list(title = paste0('PC1 (',summary(pca_res)$importance[2,1]*100,"%)")),
                        yaxis = list(title = paste0('PC2 (',summary(pca_res)$importance[2,2]*100,"%)")),
                        zaxis = list(title = paste0('PC3 (',summary(pca_res)$importance[2,3]*100,"%)"))))
  saveWidget(fig, paste0("output/plots/",outfile_name,".html"), selfcontained = F, libdir = "lib")
  
  rm(pca_res, scores, point.annotations)
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

plot_and_save_pca(pca_res, point.annotations, "1.2_log_tpms_PCA")

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

plot_and_save_pca(pca_res, point.annotations, "1.2_FC_PCA")

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
