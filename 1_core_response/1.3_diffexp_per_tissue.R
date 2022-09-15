library(stringr)
library(grid)
library(gridExtra)
library(UpSetR)
diffexp = read.csv("output/differential_expression.lfs.csv.gz")
row.names(diffexp) = diffexp$gene
diffexp = diffexp[-c(1)]

samples.annotation = read.csv("input/samples_annotation.csv")
samples.annotation[samples.annotation$sra_run_id == "",]$sra_run_id = samples.annotation[samples.annotation$sra_run_id == "",]$sample.name
samples.annotation$sample.group = make.names(samples.annotation$sample.group)
samples.annotation = samples.annotation[samples.annotation$sample.group %in% unique(str_split_fixed(names(diffexp), "_", 2)[,1]),]

# Remove recovery samples
samples.annotation = samples.annotation[samples.annotation$temperature > 30,]
# Remove e-9 pollen (their 'heat stress' is 32 degrees)
samples.annotation = samples.annotation[-grep("e9.", samples.annotation$sample.group),]

upregulated.genes = list()
downregulated.genes = list()
for(t in unique(samples.annotation$tissue)) {
  print(t)
  samples = unique(samples.annotation[samples.annotation$tissue == t,]$sample.group)
  print(samples)
  b = rowSums(diffexp[paste0(samples, "_FDR")] < 0.05 & diffexp[paste0(samples, "_logFC")] > 0)
  upregulated.genes[[t]] = names(b[b >= (0.8 * length(samples))])
  b = rowSums(diffexp[paste0(samples, "_FDR")] < 0.05 & diffexp[paste0(samples, "_logFC")] < 0)
  downregulated.genes[[t]] = names(b[b >= (0.8 * length(samples))])
}

pdf("output/plots/diffexp_genes_per_tissue.pdf")
# Number of sample groups per tissue
a = sapply(unique(samples.annotation$tissue), FUN=function(t){length(unique(samples.annotation[samples.annotation$tissue == t,]$sample.group))})
grid.table(a, rows=names(a))
upset(fromList(upregulated.genes), order.by = "freq")
grid.text("Upregulated Genes",x = 0.65, y=0.95, gp=gpar(fontsize=16))
upset(fromList(downregulated.genes), order.by = "freq")
grid.text("Downregulated Genes",x = 0.65, y=0.95, gp=gpar(fontsize=16))
dev.off()

protein.descriptions = read.csv("input/protein_descriptions.csv.gz")

up.leaves = protein.descriptions[protein.descriptions$gene %in% str_split_fixed(upregulated.genes[["leaf"]], "\\.", 2)[,1],]
write.csv(up.leaves, "output/leaf_upregulated_genes.csv", row.names = F)

# For paintomics and such @TODO (currently only quick, dirty and prelim)
# gene.names = paste(str_split_fixed(upregulated.genes$leaf, "\\.", 3)[,1], str_split_fixed(upregulated.genes$leaf, "\\.", 3)[,2], sep=".")
# write.table(gene.names, "paintomics_relevant.tsv", row.names=F, quote=F)
# leaf.samples = diffexp[grep(paste0(unique(samples.annotation[samples.annotation$tissue == "leaf",]$sample.group), collapse="|"), names(diffexp))]
# leaf.samples = leaf.samples[grep("_logFC", names(leaf.samples))]
# row.names(leaf.samples) = paste(str_split_fixed(row.names(leaf.samples), "\\.", 3)[,1], str_split_fixed(row.names(leaf.samples), "\\.", 3)[,2], sep=".")
# write.table(leaf.samples, "paintomics_vals.tsv", sep="\t", quote=F, col.names = F)

g = diffexp[row.names(diffexp) %in% upregulated.genes[["leaf"]],grep("_logFC", names(diffexp))]
meta = unique(samples.annotation[c("sample.group", "tissue", "temperature", "stress.duration", "genotype.name")])
g = g[paste0(meta$sample.group, "_logFC")]
g$gene = str_split_fixed(row.names(g), "\\.", 2)[,1]
g = merge(g, unique(protein.descriptions[c("gene", "ITAG4.1_description")]), all.y=F)
row.names(g) = paste0(gsub(" \\(AHRD.*\\)", "", g$ITAG4.1_description), " (", g$gene, ")")
library(ComplexHeatmap)
library(circlize)
library(grid)
col_fun = colorRamp2(c(-10, 0, 10), c("blue", "white", "red"))
column_ha = HeatmapAnnotation(tissue = as.factor(meta$tissue), temperature = as.numeric(meta$temperature), log.duration = log(as.numeric(meta$stress.duration)),
                              genotype = meta$genotype)
#cols = str_split_fixed(names(g), "_", 2)[,1]

pdf("output/plots/leaf_genes_heatmap.pdf", 10, 20)
Heatmap(as.matrix(g[2:18]), row_names_gp = grid::gpar(fontsize = 4), col = col_fun, bottom_annotation = column_ha)
dev.off()
