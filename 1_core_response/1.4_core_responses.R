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
library(UpSetR)
source("lib_GO_enrichment.R")
diffexp = read.csv("output/differential_expression.lfs.csv.gz")
row.names(diffexp) = diffexp$gene
diffexp = diffexp[-c(1)]

samples.annotation = read.csv("input/samples_annotation.csv")
samples.annotation[samples.annotation$sra_run_id == "",]$sra_run_id = samples.annotation[samples.annotation$sra_run_id == "",]$sample.name
samples.annotation$sample.group = make.names(samples.annotation$sample.group)
samples.annotation = samples.annotation[samples.annotation$sample.group %in% unique(str_split_fixed(names(diffexp), "_", 2)[,1]),]

protein.descriptions = read.csv("input/protein_descriptions.csv.gz")

core.genes = list(upregulated=list(), downregulated=list())
core.GOs   = list(upregulated=list(), downregulated=list())
GO.names = data.frame(GO.ID = character(), Term = character())

for(direction in c("upregulated", "downregulated")) {
  for(stress.type in c("heat", "drought", "salt")) {
    samples.stress = unique(samples.annotation[samples.annotation$stress_type == stress.type,c("sample.group", "tissue", "temperature", "stress.duration", "genotype.name")])
    diffexp.stress = diffexp[grep(paste0(c("gene", samples.stress$sample.group), collapse="|"), names(diffexp))]
    
    samples = unique(samples.stress$sample.group)
    n.without.p.val = sum(is.na(diffexp.stress[1,paste0(samples, "_FDR")])) # simply count NAs in first row
    if(direction == "upregulated")
      b = rowSums(diffexp.stress[paste0(samples, "_FDR")] < 0.05 & diffexp.stress[paste0(samples, "_logFC")] > 0, na.rm = T)
    else
      b = rowSums(diffexp.stress[paste0(samples, "_FDR")] < 0.05 & diffexp.stress[paste0(samples, "_logFC")] < 0, na.rm = T)
    genes = names(b[b >= (0.8 * (length(samples) - n.without.p.val))])
    
    core.genes[[direction]][[stress.type]] = genes
    
    # GO term enrichment
    enriched.GOs = go_enrichment(str_remove(genes, "\\.[0-9]+\\.[0-9]+$"))
    write.csv(enriched.GOs, paste0("output/1.4_core_response_",stress.type,"_",direction,"_GOs.csv"), row.names=F)
    core.GOs[[direction]][[stress.type]] = enriched.GOs$GO.ID
    GO.names = unique(rbind(GO.names, enriched.GOs[1:2]))
    
    # HSF binding sites
    fimo_matches = str_remove(read.csv("../data/fimo_hsfs_matches.tsv", sep="\t")$sequence_name, "\\([+-]\\)")
    genes_matched = str_remove(genes, "\\.[0-9]+$") %in% fimo_matches
    
    # Heatmap (right now only doing for HS)
    if(stress.type == "heat") {
      cells = diffexp.stress[row.names(diffexp.stress) %in% genes, grep("_logFC", names(diffexp.stress))]
      ## Add protein descriptions to the name
      cells$gene = str_split_fixed(row.names(cells), "\\.", 2)[,1]
      cells = merge(cells, unique(protein.descriptions[c("gene", "ITAG4.1_description")]), all.y=F)
      row.names(cells) = paste0(gsub(" \\(AHRD.*\\)", "", cells$ITAG4.1_description), " (", cells$gene, ")")
      cells = cells[ , -which(names(cells) %in% c("gene","ITAG4.1_description"))]
      ## Heatmap Annotations
      column_ha = HeatmapAnnotation(tissue = as.factor(samples.stress$tissue), temperature = as.numeric(samples.stress$temperature), log.duration = log(as.numeric(samples.stress$stress.duration)),
                                    genotype = samples.stress$genotype, col = list(tissue = c("anther" = "#CC79A7", "fruit" = "#D55E00", "leaf" = "#009E73", "pollen" = "#F0E442", "seed" = "#E69F00", "seedling" = "#56B4E9", "ovaries"="pink")))
      row_ha = rowAnnotation(hsf_binding = genes_matched, col = list(hsf_binding = c("TRUE" = "black", "FALSE" = "white")))
      
      pdf(paste0("output/plots/1.4_core_response_",stress.type,"_",direction,".lfs.pdf"), 15, 15)
      print(Heatmap(as.matrix(cells), col = colorRamp2(c(-10, 0, 10), c("blue", "white", "red")), bottom_annotation = column_ha, left_annotation = row_ha))
      dev.off()
      
    }
  }
}


# Make Venn Diagrams for genes and GO terms
pdf("output/plots/1.4_core_response_upsets.pdf", 9, 7)
print(upset(fromList(unlist(core.genes, recursive=F)), order.by = "freq", nsets=6))
grid.text("Overlapping Genes",x = 0.65, y=0.95, gp=gpar(fontsize=16))
print(upset(fromList(unlist(core.GOs, recursive=F)), order.by = "freq", nsets=6))
grid.text("Overlapping Gene Ontology Terms (BP)",x = 0.65, y=0.95, gp=gpar(fontsize=16))
dev.off()

# @TODO output lists of which genes & GO terms are shared/unique for each stress type
