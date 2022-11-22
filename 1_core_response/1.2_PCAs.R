library(stringr)
library(htmlwidgets)
library(plotly)
library(ggfortify)

samples.annotation = read.csv("input/samples_annotation.csv")
samples.annotation[samples.annotation$sra_run_id == "",]$sra_run_id = samples.annotation[samples.annotation$sra_run_id == "",]$sample.name

abundance = read.csv("input/combined_abundance.tsv.gz", sep="\t", check.names=F)
row.names(abundance) = abundance$target_id
quantified.samples = unique(str_split_fixed(names(abundance)[3:ncol(abundance)], "_", 2)[,1])

samples.annotation = samples.annotation[samples.annotation$sra_run_id %in% quantified.samples,]
samples.annotation$temperature = as.numeric(samples.annotation$temperature)

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

### ----- TPMs -----

tpms = t(as.matrix(abundance[grep("_tpm", names(abundance))]))
row.names(tpms) = str_remove(row.names(tpms), "_tpm")

pca_res <- prcomp(tpms)
scores = as.data.frame(pca_res$x)

point.annotations = data.frame(sample = row.names(scores))
point.annotations = merge(point.annotations, samples.annotation, by.x = "sample", by.y="sra_run_id", all.x=T)

plot_and_save_pca(pca_res, point.annotations, "1.2_tpms_PCA")

rm(abundance, scores, pca_res, tpms, point.annotations)

#### ---- fold-changes -------
diffexp = read.csv("output/differential_expression.lfs.csv.gz")
row.names(diffexp) = diffexp$gene

fcs = t(as.matrix(diffexp[grep("_logFC", names(diffexp))]))
row.names(fcs) = str_remove(row.names(fcs), "_logFC")
# fill in NAs with 0
fcs[is.na(fcs)] <- 0

pca_res <- prcomp(fcs)
scores = as.data.frame(pca_res$x)

point.annotations = data.frame(sample = row.names(scores))
samples.annotation$sample.group = str_replace_all(samples.annotation$sample.group, "-", ".")
samples.annotation = unique(samples.annotation[c("genotype.name", "tissue", "stress.duration", "temperature", "sample.group", "stress_type")])
point.annotations = merge(point.annotations, samples.annotation, by.x = "sample", by.y="sample.group", all.x=T, all.y=F)

plot_and_save_pca(pca_res, point.annotations, "1.2_FC_PCA")
