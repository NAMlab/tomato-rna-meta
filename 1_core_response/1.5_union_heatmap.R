# This script generates a heatmap visualizing all genes that are differentially expressed
# in *any* sample (that is, the union of all of them). To ensure we don't just plot the whole genome,
# we are applying a somewhat stricter p-value threshold: Instead of just adjusting the p-value within each
# experiment separately, we use the raw p values and then adjust them across *all* the contrasts.
# Only genes differentially expressed in any contrast using that p value are included in the heatmap.

library(stringr)
library(data.table)
library(ComplexHeatmap)
library(circlize)
library(grid)
diffexp = read.csv("output/differential_expression.lfs.csv.gz")

samples.annotation = read.csv("input/samples_annotation.csv")
samples.annotation[samples.annotation$sra_run_id == "",]$sra_run_id = samples.annotation[samples.annotation$sra_run_id == "",]$sample.name
samples.annotation$sample.group = make.names(samples.annotation$sample.group)
samples.annotation = samples.annotation[samples.annotation$sample.group %in% unique(str_split_fixed(names(diffexp), "_", 2)[,1]),]

# Remove all samples from diffexp that are not part of our samples.annotation
diffexp = diffexp[grep(paste0(c("gene", unique(samples.annotation$sample.group)), collapse="|"), names(diffexp))]

# Reshape the table to long format and calculate adjusted p-value across all samples
long = melt.data.table(data.table(diffexp), measure=patterns("_logFC$", "_PValue$", "_FDR$"), value.name=c('logFC', 'PValue', 'FDR'), variable.name="contrast", na.rm=TRUE, id='gene')
long$p.adj = p.adjust(long$PValue, method="BH")

# Now, if we take all genes with an adjusted p-value < 0.05 at least once, we get 25k genes, so still
# the majority of the whole genome. We'll do that for now and see what Andrea and JJ say.
genes = unique(long[long$p.adj < 0.05,]$gene)
#genes = names(table(genes)[table(genes) > 2])

row.names(diffexp) = diffexp$gene
diffexp = diffexp[-c(1)]

g = diffexp[row.names(diffexp) %in% genes, grep("_logFC", names(diffexp))]
meta = unique(samples.annotation[c("sample.group", "tissue", "temperature", "stress.duration", "genotype.name")])
g = g[paste0(meta$sample.group, "_logFC")]

g.p = diffexp[row.names(diffexp) %in% genes, grep("_FDR", names(diffexp))]
g.p = g.p[paste0(meta$sample.group, "_FDR")]

fimo_matches = str_remove(read.csv("../data/fimo_hsfs_matches.tsv", sep="\t")$sequence_name, "\\([+-]\\)")
genes_matched = str_remove(genes, "\\.[0-9]+$") %in% fimo_matches

col_fun = colorRamp2(c(-10, 0, 10), c("blue", "white", "red"))
column_ha = HeatmapAnnotation(tissue = as.factor(meta$tissue), temperature = as.numeric(meta$temperature), log.duration = log(as.numeric(meta$stress.duration)),
                              genotype = meta$genotype)
row_ha = rowAnnotation(hsf_binding = genes_matched, col = list(hsf_binding = c("TRUE" = "black", "FALSE" = "white")))

pdf("output/plots/union_heatmap.lfs.pdf", 10, 200)
Heatmap(as.matrix(g), show_row_names = FALSE, col = col_fun, left_annotation = row_ha, bottom_annotation = column_ha)
dev.off()
