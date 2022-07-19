library(stringr)
library(htmlwidgets)
library(plotly)

samples.annotation = read.csv("input/samples_annotation.csv")
samples.annotation[samples.annotation$sra_run_id == "",]$sra_run_id = samples.annotation[samples.annotation$sra_run_id == "",]$sample.name
abundance = read.csv("input/combined_abundance.tsv.gz", sep="\t", check.names=F)
row.names(abundance) = abundance$target_id
quantified.samples = unique(str_split_fixed(names(abundance)[3:ncol(abundance)], "_", 2)[,1])

samples.annotation = samples.annotation[samples.annotation$sra_run_id %in% quantified.samples,]

tpms = t(as.matrix(abundance[grep("_tpm", names(abundance))]))
row.names(tpms) = str_remove(row.names(tpms), "_tpm")

pca_res <- prcomp(tpms)
scores = as.data.frame(pca_res$x)


point.annotations = data.frame(sample = row.names(scores))
point.annotations = merge(point.annotations, samples.annotation, by.x = "sample", by.y="sra_run_id", all.x=T)



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
saveWidget(fig, "output/plots/1.4_pca_tpm_full.html", selfcontained = F, libdir = "lib")


# Remove the outliers
remove.points = samples.annotation[grep("e13-", samples.annotation$sample.group),]$sra_run_id
remove.points = c(remove.points, samples.annotation[grep("e2-", samples.annotation$sample.group),]$sra_run_id)
remove.points = c(remove.points, samples.annotation[grep("e14-", samples.annotation$sample.group),]$sra_run_id)

tpms = tpms[!row.names(tpms) %in% remove.points,]

pca_res <- prcomp(tpms)
scores = as.data.frame(pca_res$x)

point.annotations = data.frame(sample = row.names(scores))
point.annotations = merge(point.annotations, samples.annotation, by.x = "sample", by.y="sra_run_id", all.x=T)

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
saveWidget(fig, "output/plots/1.4_pca_tpm_zoomed_in.html", selfcontained = F, libdir = "lib")
