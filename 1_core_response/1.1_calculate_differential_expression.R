# Calculates fold-changes and corresponding p- and q values to determine differentially expressed
# genes in all the sample groups present in the 
library(edgeR)
library(stringr)

samples.annotation = read.csv("input/samples_annotation.csv")
samples.annotation[samples.annotation$sra_run_id == "",]$sra_run_id = samples.annotation[samples.annotation$sra_run_id == "",]$sample.name
abundance = read.csv("input/combined_abundance.tsv.gz", sep="\t", check.names=F)
row.names(abundance) = abundance$target_id
quantified.samples = unique(str_split_fixed(names(abundance)[3:ncol(abundance)], "_", 2)[,1])

samples.annotation = samples.annotation[samples.annotation$sra_run_id %in% quantified.samples,]

diffexp = data.frame(gene = abundance$target_id)

for(trt.group in unique(samples.annotation$sample.group)) {
  print(trt.group)
  ctrl.group = unique(samples.annotation[samples.annotation$sample.group == trt.group,]$respective_control)
  if(length(ctrl.group) > 1)
    stop(paste0("Multiple different control groups for sample group ", trt.group))
  if(ctrl.group == "") {
    print("This is a control group (skip)")
    next
  }
  trt.samples = samples.annotation[samples.annotation$sample.group == trt.group,]$sra_run_id
  ctrl.samples = samples.annotation[samples.annotation$sample.group == ctrl.group,]$sra_run_id

  if(length(trt.samples) < 2 | length(ctrl.samples) < 2) {
    print("No replication (skip)")
    next
  }
    
  counts = abundance[paste0(c(ctrl.samples, trt.samples),"_est_counts")]
  groups = factor(c(rep("control", length(ctrl.samples)), rep("treatment", length(trt.samples))))
  y <- DGEList(counts=counts,group=groups)
  y <- calcNormFactors(y)
  y <- estimateDisp(y)
  et <- exactTest(y)
  res <- topTags(et, n=Inf)
  res <- as.data.frame(res)[c(1,3,4)]
  names(res) = paste0(trt.group, "_", names(res))
  res$gene = row.names(res)
  
  diffexp = merge(diffexp, res)
}

write.csv(diffexp, "output/differential_expression.lfs.csv", row.names=F)
system("pigz -11 output/differential_expression.lfs.csv")
