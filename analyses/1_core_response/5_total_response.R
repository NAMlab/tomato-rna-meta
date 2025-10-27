# This script generates (for each stress type) a heatmap visualizing all genes that are differentially expressed
# in *any* sample (that is, the union of all of them). To ensure we don't just plot the whole genome,
# we are applying a somewhat stricter p-value threshold: Instead of just adjusting the p-value within each
# experiment separately, we use the raw p values and then adjust them across *all* the contrasts.
# Only genes differentially expressed in any contrast using that p value are included in the heatmap.
plot.heatmaps = F

library(stringr)
library(data.table)
library(ComplexHeatmap)
library(circlize)
library(grid)
library(UpSetR)
source("lib_GO_enrichment.R")
diffexp = read.csv("output/differential_expression.lfs.csv.gz")

samples.annotation = read.csv("input/samples_annotation.csv")
samples.annotation[samples.annotation$sra_run_id == "",]$sra_run_id = samples.annotation[samples.annotation$sra_run_id == "",]$sample.name
samples.annotation$sample.group = make.names(samples.annotation$sample.group)
samples.annotation = samples.annotation[samples.annotation$sample.group %in% unique(str_split_fixed(names(diffexp), "_", 2)[,1]),]

protein.descriptions = read.csv("input/protein_descriptions.csv.gz")

union.genes = list(upregulated=list(), downregulated=list())
union.GOs   = list(upregulated=list(), downregulated=list())
GO.names = data.frame(GO.ID = character(), Term = character())

for(direction in c("upregulated", "downregulated")) {
  for(stress.type in c("heat", "drought", "salt")) {
    samples.stress = unique(samples.annotation[samples.annotation$stress_type == stress.type,c("sample.group", "tissue", "temperature", "stress.duration", "genotype.name")])
    diffexp.stress = diffexp[grep(paste0(c("gene", samples.stress$sample.group), collapse="|"), names(diffexp))]
   
    # Reshape the table to long format and calculate adjusted p-value across all samples
    long = melt.data.table(data.table(diffexp.stress), measure=patterns("_logFC$", "_PValue$", "_FDR$"), value.name=c('logFC', 'PValue', 'FDR'), variable.name="contrast", na.rm=TRUE, id='gene')
    long$p.adj = p.adjust(long$PValue, method="BH")
    
    if(direction == "upregulated")
      genes = unique(long[long$p.adj < 0.05 & long$logFC > 0,]$gene)
    else
      genes = unique(long[long$p.adj < 0.05 & long$logFC < 0,]$gene)
    union.genes[[direction]][[stress.type]] = genes
    
    # GO term enrichment
    enriched.GOs = go_enrichment(str_remove(genes, "\\.[0-9]+\\.[0-9]+$"))
    write.csv(enriched.GOs, paste0("output/1.5_union_",stress.type,"_",direction,"_GOs.csv"), row.names=F)
    union.GOs[[direction]][[stress.type]] = enriched.GOs$GO.ID
    GO.names = unique(rbind(GO.names, enriched.GOs[1:2]))
    
    row.names(diffexp.stress) = diffexp.stress$gene
    diffexp.stress = diffexp.stress[-c(1)]
    
    # HSF binding sites
    fimo_matches = str_remove(read.csv("input/fimo_hsfs_matches.tsv", sep="\t")$sequence_name, "\\([+-]\\)")
    genes_matched = str_remove(genes, "\\.[0-9]+$") %in% fimo_matches
    
    # Heatmap
    if(plot.heatmaps) {
      cells = diffexp.stress[row.names(diffexp.stress) %in% genes, grep("_logFC", names(diffexp.stress))]
      ## Add protein descriptions to the name
      cells$gene = str_split_fixed(row.names(cells), "\\.", 2)[,1]
      cells = merge(cells, unique(protein.descriptions[c("gene", "ITAG4.1_description")]), all.y=F)
      row.names(cells) = paste0(gsub(" \\(AHRD.*\\)", "", cells$ITAG4.1_description), " (", cells$gene, ")")
      cells = cells[ , -which(names(cells) %in% c("gene","ITAG4.1_description"))]
      ## Heatmap Annotations
      column_ha = HeatmapAnnotation(tissue = as.factor(samples.stress$tissue), temperature = as.numeric(samples.stress$temperature), log.duration = log(as.numeric(samples.stress$stress.duration)),
                                    genotype = samples.stress$genotype, col = list(tissue = c("anther" = "#CC79A7", "fruit" = "#D55E00", "leaf" = "#009E73", "pollen" = "#F0E442", "seed" = "#E69F00", "seedling" = "#56B4E9")))
      row_ha = rowAnnotation(hsf_binding = genes_matched, col = list(hsf_binding = c("TRUE" = "black", "FALSE" = "white")))
      
      pdf(paste0("output/plots/1.5_union_",stress.type,"_",direction,".lfs.pdf"), 10, 450)
      print(Heatmap(as.matrix(cells), col = colorRamp2(c(-10, 0, 10), c("blue", "white", "red")), bottom_annotation = column_ha, left_annotation = row_ha, 
                    row_names_gp = gpar(fontsize = 1)))
      dev.off()
    }
  }
}

# Prettify names
unlisted.genes = unlist(union.genes, recursive=F)
names(unlisted.genes) = c("Heat (up)", "Drought (up)", "Salt (up)", "Heat (down)", "Drought (down)", "Salt (down)")
unlisted.GOs = unlist(union.GOs, recursive=F)
names(unlisted.GOs) = c("Heat (up)", "Drought (up)", "Salt (up)", "Heat (down)", "Drought (down)", "Salt (down)")

# Make Upset Plots for genes and GO terms
pdf("output/plots/1.5_union_upsets.pdf", 5, 3)
print(upset(fromList(unlisted.genes), order.by = "freq", nsets=6, point.size=1.5, nintersects = 24, text.scale = 1, show.numbers = "no"))
grid.text("Overlapping Genes",x = 0.65, y=0.95, gp=gpar(fontsize=12))
print(upset(fromList(unlisted.GOs), order.by = "freq", nsets=6, point.size=1.5, nintersects = 24, text.scale = 1, show.numbers = "no"))
grid.text("Overlapping Gene\n Ontology Terms",x = 0.65, y=0.9, gp=gpar(fontsize=12))
dev.off()