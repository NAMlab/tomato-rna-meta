library(topGO)
set.seed(845052)

leaf.upregulated = read.csv("output/leaf_upregulated_genes.csv")
geneID2GO <- readMappings("input/GOMAP_annotation.map.lfs.gz")

sig_gene_list <- factor(as.integer(names(geneID2GO) %in% leaf.upregulated$gene))
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

tableGOresults <- tableGOresults[tableGOresults$classicFisher<0.05,]
tableGOresults$q_value <- p.adjust(tableGOresults$classicFisher, method="BH")

write.csv(tableGOresults, "output/leaf_up_GO_enrichment.csv", row.names=F)

pdf("output/plots/leaf_up_GO_enrichment.pdf", 9, 7)
showSigOfNodes(GOdata_sig, score(resultsGOfisher), firstSigNodes = 50, useInfo = 'all')
dev.off()

for(attempt.number in c(21, 31, 32, 33, 34, 35)) {
  print(paste(".. Attempt Number", attempt.number))
  important.genes = read.csv(paste0("../Attempt_Results/Attempt_",attempt.number,".csv"))
  # Remove the ones that only Stephan found but not me
  important.genes = unique(important.genes[!is.na(important.genes$medianImp),1:3])
  
  network.genes = unique(as.character(important.genes$gene))
  
  # Extract Trait Category
  important.genes$trait.category = sapply(as.character(important.genes$trait.x), function(x) { paste(strsplit(x, "\\.")[[1]][2:3], collapse=".") }, USE.NAMES = F)
  # For geometry traits, spectrum doesn't matter
  important.genes[grepl("geometry.*", important.genes$trait.category),]$trait.category = "geometry"
  
  analyze_enrichment <- function(sig_gene_list, filename) {
    print(paste("....ANALYZING... ", filename))
    GOdata_sig <- new("topGOdata",
                      description = "name",
                      ontology = b,
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
    
    tableGOresults <- tableGOresults[tableGOresults$classicFisher<0.05,]
    
    tableGOresults$q_value <- p.adjust(tableGOresults$classicFisher, method="BY")
    
    #tableGOresults <- tableGOresults[tableGOresults$q_value<0.05,]
    
    if(nrow(tableGOresults)>0){
      
      sig_list <- names(sig_gene_list)[sig_gene_list == 1]
      tableGOresults$all_ids <- NULL
      tableGOresults$sig_ids <- NA
      
      for (a in 1:length(tableGOresults$GO.ID)){
        all_ids <- genesInTerm(GOdata_sig,tableGOresults[a,"GO.ID"])
        sig_ids <- intersect(all_ids[[1]],sig_list)
        tableGOresults$sig_ids[a] <- paste(sig_ids, collapse=";")
      }
      
      if(length(levels(sig_gene_list))>1){
        # Sort rows by ascending q value
        tableGOresults = tableGOresults[order(tableGOresults$q_value),]
        
        dir.create(dirname(filename), showWarnings = F, recursive = T)
        write.csv(tableGOresults, file=paste0(filename, ".csv"), row.names=FALSE)
        pdf(paste0(filename, ".pdf"))
        showSigOfNodes(GOdata_sig, score(resultsGOfisher), firstSigNodes = 15, useInfo = 'all')
        dev.off()
      }
    }
  }
  
  upregulated.genes = c(1,2,3,4,8,9)
  downregulated.genes = c(5,6,10,7)
  allgenes = c(upregulated.genes, downregulated.genes)
  regulation.types = c("upregulated.genes", "downregulated.genes", "allgenes")
  
  gene.sets = c("all.genes", "diffexp.genes", "network.genes")
  onto <- c("BP", "MF", "CC")
  
  for(regulation.type in regulation.types) {
    if(regulation.type == "allgenes")
      xregulated.genes = diffexp.genes
    else
      xregulated.genes = unique(as.character(clusters[clusters$cluster %in% get(regulation.type),]$target.id))
    
    for(gene.set.name in gene.sets) {
      gene.set = get(gene.set.name)
      
      for(b in onto){
        
        # Differentially expressed genes (even more than "ALL GENES" --> this is the set of genes that Boruta chose from)
        if(gene.set.name == "all.genes") {
          intersect_Genes <- intersect(diffexp.genes, xregulated.genes)
          # @TODO why are we using all genes in the mothertable as the basis here?
          sig_gene_list <- factor(as.integer(gene.set %in% intersect_Genes))
          
          names(sig_gene_list) <- gene.set
          
          analyze_enrichment(sig_gene_list, paste0("results/go_enrichment/Attempt_",attempt.number,"/",regulation.type,"/against_all.genes/diffexp_genes/", b))
        }
        
        # ALL GENES IN THE NETWORK
        if((gene.set.name != "diffexp.genes" | regulation.type != "allgenes") & gene.set.name != "network.genes") { # Those two sets are the same
          intersect_Genes <- intersect(unique(as.character(important.genes$gene)), xregulated.genes)
          # @TODO why are we using all genes in the mothertable as the basis here?
          sig_gene_list <- factor(as.integer(gene.set %in% intersect_Genes))
          
          names(sig_gene_list) <- gene.set
          
          analyze_enrichment(sig_gene_list, paste0("results/go_enrichment/Attempt_",attempt.number,"/",regulation.type,"/against_",gene.set.name,"/complete_network/", b))
        }
        
        # Per subnetwork
        for(category in unique(as.character(important.genes$trait.category))) {
          intersect_Genes <- intersect(unique(as.character(important.genes[important.genes$trait.category == category,]$gene)), xregulated.genes)
          # @TODO why are we using all genes in the mothertable as the basis here?
          sig_gene_list <- factor(as.integer(gene.set %in% intersect_Genes))
          
          names(sig_gene_list) <- gene.set
          
          analyze_enrichment(sig_gene_list, paste0("results/go_enrichment/Attempt_",attempt.number,"/",regulation.type,"/against_",gene.set.name,"/per_subnet/", b, "/", category))
        }
        
        ## Do VIS subnet without rgb.blue (because that's not part of the subnet)
        intersect_Genes <- intersect(unique(as.character(important.genes[important.genes$trait.category == "intensity.vis" & important.genes$trait.x != "top.intensity.vis.rgb.blue.mean",]$gene)),
                                     xregulated.genes)
        # @TODO why are we using all genes in the mothertable as the basis here?
        sig_gene_list <- factor(as.integer(gene.set %in% intersect_Genes))
        
        names(sig_gene_list) <- gene.set
        
        analyze_enrichment(sig_gene_list, paste0("results/go_enrichment/Attempt_",attempt.number,"/",regulation.type,"/against_",gene.set.name,"/per_subnet/", b, "/intensity.vis_without_rgb_blue"))
        
        # Per trait
        for(trait in unique(as.character(important.genes$trait.x))) {
          intersect_Genes <- intersect(unique(as.character(important.genes[important.genes$trait.x == trait,]$gene)), xregulated.genes)
          # @TODO why are we using all genes in the mothertable as the basis here?
          sig_gene_list <- factor(as.integer(gene.set %in% intersect_Genes))
          
          names(sig_gene_list) <- gene.set
          
          analyze_enrichment(sig_gene_list, paste0("results/go_enrichment/Attempt_",attempt.number,"/",regulation.type,"/against_",gene.set.name,"/per_trait/", b, "/", trait))
          
          
        }
      }
    }
  }
}