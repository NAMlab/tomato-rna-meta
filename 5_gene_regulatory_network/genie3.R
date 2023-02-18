library(GENIE3)
library(stringr)

abundance = read.csv("../data/combined_abundance.lfs.tsv.gz", sep="\t", check.names=F)
row.names(abundance) = abundance$target_id
abundance = abundance[grep("_tpm", names(abundance))]
names(abundance) = str_remove(names(abundance), "_tpm")

exprMatr = as.matrix(abundance)

set.seed(2023)
# Step 1: Find the probable regulators of our core HS genes
core.genes = read.csv("../1_core_response/output/1.4_core_response_genes.csv")
targets = core.genes[core.genes$stress.type == "heat",]$target.id
weightMat <- GENIE3(exprMatr, targets = targets, verbose=T, nCores = 20)
linkList <- getLinkList(weightMat, reportMax = 100)

write.csv(linkList, "regulators.csv", row.names=F)


# Step 2: Find other genes probably regulated by our found regulators



