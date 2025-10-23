# This script determines the core heat stress response for each stress type (heat, drought, salt)
# by identifying genes that are differentially expressed in at least 80% of the samples for that
# stress type. It then performs GO term enrichment analysis on these genes and generates a heatmap
# of the core response genes across all samples. Finally, it creates Upset plots to visualize
# the overlap of core response genes and GO terms between different stress types.

library(stringr)
library(data.table)
library(ComplexHeatmap)
library(circlize)
library(grid)

diffexp = read.csv("output/differential_expression.lfs.csv.gz")
row.names(diffexp) = diffexp$gene
diffexp = diffexp[-c(1)]

samples.annotation = read.csv("input/samples_annotation.csv")
samples.annotation[samples.annotation$sra_run_id == "",]$sra_run_id = samples.annotation[samples.annotation$sra_run_id == "",]$sample.name
samples.annotation$sample.group = make.names(samples.annotation$sample.group)
samples.annotation = samples.annotation[samples.annotation$sample.group %in% unique(str_split_fixed(names(diffexp), "_", 2)[,1]),]
samples.annotation = samples.annotation[samples.annotation$stress_type == "heat",]

diffexp = diffexp[grep(paste0(c("gene", samples.annotation$sample.group), collapse="|"), names(diffexp))]

hsfs = read.csv("input/hsfs.csv")

# Make a heatmap of HS core genes across all the samples
cells = diffexp[row.names(diffexp) %in% hsfs$target_id,
                       grep("_logFC", names(diffexp))]
names(cells) = str_remove(names(cells), "_logFC")

cells$gene = row.names(cells)
cells = merge(cells, hsfs, sort=F, by.x="row.names", by.y="target_id")
row.names(cells) = cells$transcription_factor
cells = cells[ , -which(names(cells) %in% c("Row.names","gene", "gene_id", "transcription_factor"))]

column_data = unique(samples.annotation[c("sample.group", "tissue", "temperature", "stress.duration", "genotype.name")])

# Order sample groups by stress.duration
ordered_groups <- column_data$sample.group[order(as.numeric(column_data$stress.duration))]
column_data <- column_data[match(ordered_groups, column_data$sample.group), ]
cells <- cells[, ordered_groups]

column_ha = HeatmapAnnotation(
  stress_duration = anno_text(
    ifelse(is.na(column_data$stress.duration), "?", as.character(as.numeric(column_data$stress.duration))),
    rot = 90,
    gp = gpar(fontsize = 10),
    show_name = TRUE
  ),
  tissue = as.factor(column_data$tissue),
  temperature = anno_text(
    ifelse(is.na(column_data$temperature), "?", as.character(as.numeric(column_data$temperature))),
    rot = 90,
    gp = gpar(fontsize = 10),
    show_name = TRUE
  ),
  genotype = column_data$genotype,
  col = list(
    tissue = c(
      "anther" = "#CC79A7", "fruit" = "#D55E00", "leaf" = "#009E73",
      "pollen" = "#F0E442", "seed" = "#E69F00", "seedling" = "#56B4E9",
      "ovaries" = "pink", "root" = "brown"
    )
  )
)

names(cells) = str_replace_all(names(cells), "\\.", "-")
pdf("output/plots/1.7_known_hsfs_heatmap.pdf", 12, 6)
draw(Heatmap(as.matrix(cells), col = colorRamp2(c(-10, 0, 10), c("blue", "white", "red")),
             row_names_max_width = max_text_width(rownames(cells), gp = gpar(fontsize = 11)), row_names_gp = grid::gpar(fontsize = 11),
             bottom_annotation = column_ha, cluster_columns = FALSE,
             heatmap_legend_param = list(title = "logFC")),
     heatmap_legend_side="right", annotation_legend_side = 'right')
dev.off()
