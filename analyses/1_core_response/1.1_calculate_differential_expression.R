# Calculates fold-changes and corresponding p- and q values to determine differentially expressed
# genes in all the sample groups present in the data set.
library(edgeR)
library(stringr)

samples.annotation = read.csv("input/samples_annotation.csv")
samples.annotation[samples.annotation$sra_run_id == "",]$sra_run_id = samples.annotation[samples.annotation$sra_run_id == "",]$sample.name
abundance = read.csv("input/combined_abundance.tsv.gz", sep="\t", check.names=F)
row.names(abundance) = abundance$target_id

diffexp = data.frame(gene = abundance$target_id)

for(trt.group in unique(samples.annotation$sample.group)) {
  print(trt.group)
  ctrl.group = unique(samples.annotation[samples.annotation$sample.group == trt.group,]$respective_control)
  if(length(ctrl.group) > 1)
    stop(paste0("Multiple different control groups for sample group ", trt.group))
  if(ctrl.group == "" || ctrl.group == "X") {
    print("This is a control group (skip)")
    next
  }
  trt.samples = samples.annotation[samples.annotation$sample.group == trt.group,]$sra_run_id
  ctrl.samples = samples.annotation[samples.annotation$sample.group == ctrl.group,]$sra_run_id
  
  counts = abundance[paste0(c(ctrl.samples, trt.samples),"_est_counts")]
  groups = factor(c(rep("control", length(ctrl.samples)), rep("treatment", length(trt.samples))))
  y <- DGEList(counts=counts,group=groups)
  y <- calcNormFactors(y)
  
  if(length(trt.samples) < 2 | length(ctrl.samples) < 2) {
    # When there are no replicates, we just assume a dispersion in order to be able to calculate
    # our statistics. The dispersion has no influence on the logFC, but it does on the p-value.
    # Therefore, since our disperion might be wrong, we will set the p-value to NA in cases when there
    # are no replicates (see below) but at least we have a logFC.
    et <- exactTest(y, dispersion = 0.2^2)
  } else {
    y <- estimateDisp(y)
    et <- exactTest(y)
  }
  
  res <- topTags(et, n=Inf)
  res <- as.data.frame(res)[c(1,3,4)]
  names(res) = paste0(trt.group, "_", names(res))
  res$gene = row.names(res)
  
  # If we don't have replicates, we cannot know the p-value (see above), so we're setting it to NA
  if(length(trt.samples) < 2 | length(ctrl.samples) < 2)
    res[2:3] <- NA
  
  diffexp = merge(diffexp, res)
}

write.csv(diffexp, "output/differential_expression.lfs.csv", row.names=F)
system("pigz -11 output/differential_expression.lfs.csv") # compress
