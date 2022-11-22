# This script generates (for each stress type) a heatmap visualizing all genes that are differentially expressed
# in *any* sample (that is, the union of all of them). To ensure we don't just plot the whole genome,
# we are applying a somewhat stricter p-value threshold: Instead of just adjusting the p-value within each
# experiment separately, we use the raw p values and then adjust them across *all* the contrasts.
# Only genes differentially expressed in any contrast using that p value are included in the heatmap.

library(stringr)
library(data.table)
library(ComplexHeatmap)
library(circlize)
library(grid)
library(VennDiagram)
source("lib_GO_enrichment.R")
diffexp = read.csv("output/differential_expression.lfs.csv.gz")

samples.annotation = read.csv("input/samples_annotation.csv")
samples.annotation[samples.annotation$sra_run_id == "",]$sra_run_id = samples.annotation[samples.annotation$sra_run_id == "",]$sample.name
samples.annotation$sample.group = make.names(samples.annotation$sample.group)
samples.annotation = samples.annotation[samples.annotation$sample.group %in% unique(str_split_fixed(names(diffexp), "_", 2)[,1]),]

protein.descriptions = read.csv("input/protein_descriptions.csv.gz")

union.genes = list()
union.GOs   = list()
GO.names = data.frame(GO.ID = character(), Term = character())

for(stress.type in c("heat", "drought", "salt")) {
  samples.stress = unique(samples.annotation[samples.annotation$stress_type == stress.type,c("sample.group", "tissue", "temperature", "stress.duration", "genotype.name")])
  diffexp.stress = diffexp[grep(paste0(c("gene", samples.stress$sample.group), collapse="|"), names(diffexp))]
 
  # Reshape the table to long format and calculate adjusted p-value across all samples
  long = melt.data.table(data.table(diffexp.stress), measure=patterns("_logFC$", "_PValue$", "_FDR$"), value.name=c('logFC', 'PValue', 'FDR'), variable.name="contrast", na.rm=TRUE, id='gene')
  long$p.adj = p.adjust(long$PValue, method="BH")
  
  genes = unique(long[long$p.adj < 0.05,]$gene)
  union.genes[[stress.type]] = genes
  
  # GO term enrichment
  enriched.GOs = go_enrichment(str_remove(genes, "\\.[0-9]+\\.[0-9]+$"))
  write.csv(enriched.GOs, paste0("output/1.3_union_",stress.type,"_GOs.csv"), row.names=F)
  union.GOs[[stress.type]] = enriched.GOs$GO.ID
  GO.names = unique(rbind(GO.names, enriched.GOs[1:2]))
  
  row.names(diffexp.stress) = diffexp.stress$gene
  diffexp.stress = diffexp.stress[-c(1)]
  
  # HSF binding sites
  fimo_matches = str_remove(read.csv("../data/fimo_hsfs_matches.tsv", sep="\t")$sequence_name, "\\([+-]\\)")
  genes_matched = str_remove(genes, "\\.[0-9]+$") %in% fimo_matches
  
  # Heatmap
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
  
  pdf(paste0("output/plots/1.3_union_",stress.type,".lfs.pdf"), 10, 450)
  print(Heatmap(as.matrix(cells), col = colorRamp2(c(-10, 0, 10), c("blue", "white", "red")), bottom_annotation = column_ha, left_annotation = row_ha, 
          row_names_gp = gpar(fontsize = 1)))
  dev.off()
}

# Make Venn Diagrams for genes and GO terms
pdf("output/plots/1.3_union_venns.pdf")
grid.draw(venn.diagram(union.genes, filename=NULL, main="Gene Overlaps", disable.logging = T))
grid.newpage()
grid.draw(venn.diagram(union.GOs, filename=NULL, main="GO Overlaps", disable.logging = T))
dev.off()

# @TODO output lists of which genes & GO terms are shared/unique for each stress type
