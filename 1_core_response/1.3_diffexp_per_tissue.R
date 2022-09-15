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


protein.descriptions = read.csv("input/protein_descriptions.csv.gz")
protein.descriptions$ITAG4.1_description = sapply(lapply(protein.descriptions$ITAG4.1_description, strwrap, width=50), paste, collapse="\n")
protein.descriptions$OMA_orthologues = sapply(lapply(protein.descriptions$OMA_orthologues, strwrap, width=50), paste, collapse="\n")

#plot.new()
#always.up = Reduce(intersect, upregulated.genes)
#grid.table(protein.descriptions[protein.descriptions$gene %in% str_split_fixed(always.up, "\\.", 2)[,1],], theme=ttheme_minimal(base_size = 6))

plot.new()
grid.table(protein.descriptions[protein.descriptions$gene %in% str_split_fixed(upregulated.genes[["leaf"]], "\\.", 2)[,1],], theme=ttheme_minimal(base_size = 6))
dev.off()

up.leaves = protein.descriptions[protein.descriptions$gene %in% str_split_fixed(upregulated.genes[["leaf"]], "\\.", 2)[,1],]
write.csv(up.leaves, "output/leaf_upregulated_genes.csv", row.names = F)

# For paintomics and such @TODO (currently only quick, dirty and prelim)
# gene.names = paste(str_split_fixed(upregulated.genes$leaf, "\\.", 3)[,1], str_split_fixed(upregulated.genes$leaf, "\\.", 3)[,2], sep=".")
# write.table(gene.names, "paintomics_relevant.tsv", row.names=F, quote=F)
# leaf.samples = diffexp[grep(paste0(unique(samples.annotation[samples.annotation$tissue == "leaf",]$sample.group), collapse="|"), names(diffexp))]
# leaf.samples = leaf.samples[grep("_logFC", names(leaf.samples))]
# row.names(leaf.samples) = paste(str_split_fixed(row.names(leaf.samples), "\\.", 3)[,1], str_split_fixed(row.names(leaf.samples), "\\.", 3)[,2], sep=".")
# write.table(leaf.samples, "paintomics_vals.tsv", sep="\t", quote=F, col.names = F)
