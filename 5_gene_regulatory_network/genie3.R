library(GENIE3)
library(stringr)

abundance = read.csv("../data/combined_abundance.lfs.tsv.gz", sep="\t", check.names=F)
row.names(abundance) = abundance$target_id
abundance = abundance[grep("_tpm", names(abundance))]
names(abundance) = str_remove(names(abundance), "_tpm")

exprMatr = as.matrix(abundance)

set.seed(2023)
# Try to find the 5 most probable regulators for each HS core gene
core.genes = read.csv("../1_core_response/output/1.4_core_response_genes.csv")
targets = core.genes[core.genes$stress.type == "heat",]$target.id

linkList = data.frame(regulatoryGene=character(), targetGene=character(), weight=numeric())
for(t in targets) {
  weightMat <- GENIE3(exprMatr, targets = t, verbose=T, nCores=20)
  linkList <- rbind(linkList, getLinkList(weightMat, reportMax = 3))
  write.csv(linkList, "top5_regulators.csv", row.names=F)
  #system("Rscript 2_plot.R")
}



