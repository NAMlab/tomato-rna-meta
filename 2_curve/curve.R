library(stringr)
library(Rmisc)

annotated.genes = read.csv("output/annotated_genes.csv")

diffexp = read.csv("../1_core_response/output/differential_expression.lfs.csv.gz")
row.names(diffexp) = diffexp$gene
samples.annotation = read.csv("../1_core_response/input/samples_annotation.csv")
samples.annotation[samples.annotation$sra_run_id == "",]$sra_run_id = samples.annotation[samples.annotation$sra_run_id == "",]$sample.name
samples.annotation$sample.group = make.names(samples.annotation$sample.group)

# Remove contrasts without replicates
b = table(samples.annotation$sample.group)
samples.annotation = samples.annotation[!samples.annotation$sample.group %in% names(b[b < 2]),]

plt = data.frame(n.contrasts=numeric(), n.genes.upper=numeric(), n.genes.mean=numeric(), n.genes.lower=numeric(), 
                  tp.upper=numeric(), tp.mean=numeric(), tp.lower=numeric(), 
                  fp.upper=numeric(), fp.mean=numeric(), fp.lower=numeric(),
                  stress=character())
set.seed(3528)
for(stress in c("heat", "drought", "salt")) {
  samples.stress = samples.annotation[samples.annotation$stress_type == stress,]
  
  contrasts = unique(samples.stress$sample.group)
  
  res = data.frame(run = numeric(), n.contrasts = numeric(), n.genes = numeric(), true.pos = numeric(), false.pos = numeric())
  for(run in 1:50) {
    print(run)
    contrasts.queue = sample(contrasts)
    for(n.contrasts in 1:length(contrasts.queue)) {
      b = rowSums(diffexp[paste0(contrasts.queue[1:n.contrasts], "_FDR")] < 0.05 & diffexp[paste0(contrasts.queue[1:n.contrasts], "_logFC")] > 0, na.rm = T)
      found.genes = str_split_fixed(names(b[b >= (0.8 * n.contrasts)]), "\\.", 2)[,1]
      #found.genes = str_split_fixed(names(b[b >= 1]), "\\.", 2)[,1]
      true.pos = sum(found.genes %in% annotated.genes[annotated.genes$stress == stress & annotated.genes$reasonable_term == "true",]$gene)
      false.pos = sum(found.genes %in% annotated.genes[annotated.genes$stress == stress & annotated.genes$reasonable_term == "false",]$gene)
      res = rbind(res, list(run = run, n.contrasts = n.contrasts, n.genes = length(found.genes), true.pos = true.pos, false.pos=false.pos))
    }
  }
  # Number of genes in intersect set
  ci.genes = group.CI(n.genes ~ n.contrasts, res)
  
  # Number of Genes annotated with a specific term
  res$tp = (res$true.pos / (res$n.genes)) * 100 # [%]
  ci.tp = group.CI(tp ~ n.contrasts, res)
  
  res$fp = (res$false.pos / (res$n.genes)) * 100 # [%]
  ci.fp = group.CI(fp ~ n.contrasts, res)
  
  
  m = merge(ci.genes, merge(ci.tp, ci.fp))
  m$stress = stress
  plt = rbind(plt, m)
}

##### Plots #####
pdf("output/curves.pdf")
stress.cols = list(heat="red",drought="blue",salt="green")
# N Genes Plot
plot(plt$n.genes.mean ~ plt$n.contrasts, type="n", xlab="# of contrasts", las=2, ylab="# of \"core\" genes", ylim=c(min(plt$n.genes.lower), max(plt$n.genes.upper)))
for(stress in c("heat", "drought", "salt")) {
  d = plt[plt$stress == stress,]
  polygon(c(d$n.contrasts,rev(d$n.contrasts)),c(d$n.genes.upper,rev(d$n.genes.lower)),col = rgb(0.9,0.9,0.9,0.5), border = FALSE)
  lines(d$n.genes.mean, type="l", col=stress.cols[[stress]])
}


plot(plt$tp.mean ~ plt$n.contrasts, type="n", xlab="# of contrasts", las=2, ylab="Genes Annotated [%]", ylim=c(0,max(plt$tp.upper)))
for(stress in c("heat", "drought", "salt")) {
  d = plt[plt$stress == stress,]
  # True positives
  polygon(c(d$n.contrasts,rev(d$n.contrasts)),c(d$tp.upper,rev(d$tp.lower)),col = rgb(0.9,0.9,0.9,0.5), border = FALSE)
  lines(d$tp.mean, type="l", col=stress.cols[[stress]])
  
  # False positives
  polygon(c(d$n.contrasts,rev(d$n.contrasts)),c(d$fp.upper,rev(d$fp.lower)),col = rgb(0.9,0.9,0.9,0.5), border = FALSE)
  lines(d$fp.mean, type="l", col=stress.cols[[stress]], lty=2)
}

## Do a ROC curve for what is the optimal cutoff (right now we're just using 80% because we think it's right)
# After thinking about it more, this is probably nonsense because we actually do not want all the genes annotated
# with "response to heat" in our set, because also a gene which only responds under certain circumstances should be annotated
# but we actually do not want it in our set. Therefor maximising recall is actually nonsense here.
roc = data.frame(stress=character(), p = numeric(), n.genes = numeric(), tpr = numeric(), fpr = numeric())
for(stress in c("heat", "drought", "salt")) {
  contrasts = unique(samples.stress$sample.group)
  b = rowSums(diffexp[paste0(contrasts, "_FDR")] < 0.05 & diffexp[paste0(contrasts, "_logFC")] > 0, na.rm = T)

  for(p in seq(0.05, 0.95, 0.05)) {
    found.genes = str_split_fixed(names(b[b >= (p * length(contrasts))]), "\\.", 2)[,1]
    true.pos = sum(found.genes %in% annotated.genes[annotated.genes$stress == stress & annotated.genes$reasonable_term == "true",]$gene)
    fpr = (length(found.genes) - true.pos) / (nrow(diffexp) - length(annotated.genes[annotated.genes$stress == stress & annotated.genes$reasonable_term == "true",]$gene))
    roc = rbind(roc, list(stress=stress, p = p, n.genes = length(found.genes), tpr = true.pos/length(annotated.genes[annotated.genes$stress == stress & annotated.genes$reasonable_term == "true",]$gene), fpr = fpr))
  }
}
plot(roc$tpr ~ roc$fpr, xlab="False Positive Rate (1 - Specificity)", ylab="True Positive Rate", type="n")
for(stress in c("heat", "drought", "salt")) {
  d = roc[roc$stress == stress,]
  lines(d$tpr ~ d$fpr, col=stress.cols[[stress]])
  d$youden = d$tpr - d$fpr
  points(d[d$youden == max(d$youden),]$tpr ~ d[d$youden == max(d$youden),]$fpr, col="black", pch=16)
}

dev.off()
