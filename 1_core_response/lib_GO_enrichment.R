library(topGO)

geneID2GO <- readMappings("input/GOMAP_annotation.map.lfs.gz")

go_enrichment <- function(genes, seed = 845052) {
  set.seed(seed)
  
  sig_gene_list <- factor(as.integer(names(geneID2GO) %in% genes))
  names(sig_gene_list) = names(geneID2GO)
  
  GOdata_sig <- new("topGOdata",
                    description = "name",
                    ontology = "BP",
                    allGenes = sig_gene_list,
                    nodeSize = 10,
                    annot = annFUN.gene2GO,
                    gene2GO = geneID2GO)
  resultsGOfisher <- runTest(GOdata_sig, algorithm = "classic",
                             statistic = "fisher")
  tableGOresults <- GenTable(GOdata_sig,
                             classicFisher=resultsGOfisher,
                             topNodes=length(resultsGOfisher@score))
  
  tableGOresults$classicFisher <- ifelse(tableGOresults$classicFisher == "< 1e-30", 1e-30, tableGOresults$classicFisher)
  tableGOresults$classicFisher <- as.numeric(tableGOresults$classicFisher)
  
  tableGOresults$q_value <- p.adjust(tableGOresults$classicFisher, method="BH")
  tableGOresults <- tableGOresults[tableGOresults$q_value<0.05,]
  
  
  return(tableGOresults)
}
