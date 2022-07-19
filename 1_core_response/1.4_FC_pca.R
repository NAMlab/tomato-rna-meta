library(stringr)
library(htmlwidgets)
library(plotly)

diffexp = read.csv("output/differential_expression.lfs.csv.gz")
row.names(diffexp) = diffexp$gene

samples.annotation = read.csv("input/samples_annotation.csv")
samples.annotation[samples.annotation$sra_run_id == "",]$sra_run_id = samples.annotation[samples.annotation$sra_run_id == "",]$sample.name

fcs = t(as.matrix(diffexp[grep("_logFC", names(diffexp))]))
row.names(fcs) = str_remove(row.names(fcs), "_logFC")
# fill in NAs with 0
fcs[is.na(fcs)] <- 0

pca_res <- prcomp(fcs)
scores = as.data.frame(pca_res$x)

point.annotations = data.frame(sample = row.names(scores))
samples.annotation$sample.group = str_replace_all(samples.annotation$sample.group, "-", ".")
samples.annotation = unique(samples.annotation[c("genotype.name", "tissue", "stress.duration", "temperature", "sample.group")])
point.annotations = merge(point.annotations, samples.annotation, by.x = "sample", by.y="sample.group", all.x=T, all.y=F)



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
saveWidget(fig, "output/plots/1.5_pca_fc_all.html", selfcontained = F, libdir = "lib")

# Set all non-significant fold-changes to 0
fdrs = t(as.matrix(diffexp[grep("_FDR", names(diffexp))]))
fcs[fdrs >= 0.05] <- 0

pca_res <- prcomp(fcs)
scores = as.data.frame(pca_res$x)

point.annotations = data.frame(sample = row.names(scores))
samples.annotation$sample.group = str_replace_all(samples.annotation$sample.group, "-", ".")
samples.annotation = unique(samples.annotation[c("genotype.name", "tissue", "stress.duration", "temperature", "sample.group")])
point.annotations = merge(point.annotations, samples.annotation, by.x = "sample", by.y="sample.group", all.x=T, all.y=F)



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
saveWidget(fig, "output/plots/1.5_pca_fc_only_significant.html", selfcontained = F, libdir = "lib")
