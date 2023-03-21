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

plt$n.genes.lower = log(plt$n.genes.lower, base=10)
plt$n.genes.mean = log(plt$n.genes.mean, base=10)
plt$n.genes.upper = log(plt$n.genes.upper, base=10)

## Do the number of genes and % annotated genes again but rather than adding more contrasts here we are
# just lowering the p-value.
plt.p = data.frame(p.value=numeric(), n.genes.upper=numeric(), n.genes.mean=numeric(), n.genes.lower=numeric(), 
                 tp.upper=numeric(), tp.mean=numeric(), tp.lower=numeric(), 
                 fp.upper=numeric(), fp.mean=numeric(), fp.lower=numeric(),
                 stress=character())
for(stress in c("heat", "drought", "salt")) {
  samples.stress = samples.annotation[samples.annotation$stress_type == stress,]
  
  contrasts = unique(samples.stress$sample.group)
  
  res = data.frame(p.val = numeric(), n.genes=numeric(), true.pos=numeric(), false.pos=numeric())
  for(p.val in c(1e-01, 0.05, 0.01, 1e-05, 1e-10, 1e-15, 1e-20, 1e-25, 1e-30, 1e-40, 1e-50, 1e-60, 1e-70, 1e-80, 1e-100)) {
    print(p.val)
    #@TODO speed this up by applying it across all columns rather than looping
    for(co in contrasts) {
      found.genes = str_split_fixed(diffexp[diffexp[paste0(co, "_FDR")] < p.val & diffexp[paste0(co, "_logFC")] > 0,]$gene, "\\.", 2)[,1]
      true.pos = sum(found.genes %in% annotated.genes[annotated.genes$stress == stress & annotated.genes$reasonable_term == "true",]$gene)
      false.pos = sum(found.genes %in% annotated.genes[annotated.genes$stress == stress & annotated.genes$reasonable_term == "false",]$gene)
      res = rbind(res, list(p.val=p.val, n.genes=length(found.genes), true.pos=true.pos, false.pos=false.pos))
    }
  }

  # Number of genes in intersect set
  ci.genes = group.CI(n.genes ~ p.val, res)
  
  # Number of Genes annotated with a specific term
  res$tp = (res$true.pos / (res$n.genes)) * 100 # [%]
  ci.tp = group.CI(tp ~ p.val, res)
  
  res$fp = (res$false.pos / (res$n.genes)) * 100 # [%]
  ci.fp = group.CI(fp ~ p.val, res)
  
  
  m = merge(ci.genes, merge(ci.tp, ci.fp))
  m$stress = stress
  plt.p = rbind(plt.p, m)
}

plt.p$p.val = -log(plt.p$p.val, base=10)
plt.p$n.genes.lower = log(plt.p$n.genes.lower, base=10)
plt.p$n.genes.mean = log(plt.p$n.genes.mean, base=10)
plt.p$n.genes.upper = log(plt.p$n.genes.upper, base=10)

##### Plots #####
stress.cols = list(heat="#D55E00", drought="#0072B2", salt="#009E73")

pdf("output/curves.pdf", 6.5, 5)
par(mfrow=c(2,2), family="serif", cex=0.6)

# N Genes Plots
par(mar=c(2,3,0,0))
plot(plt$n.genes.mean ~ plt$n.contrasts, type="n", ylim=c(0, 4), axes=F)
axis(side = 1, las=1, lwd=0, lwd.ticks = 1, labels=F)
axis(side = 2, las=1, lwd=0, lwd.ticks = 1, at = 0:4, labels=c(1,10,100,"1k","10k"))
text(max(plt$n.contrasts), 3.6, "Number of selected genes \nby number of contrasts (50 simulations)", pos=2, cex=1.2)
for(stress in c("heat", "drought", "salt")) {
  d = plt[plt$stress == stress,]
  polygon(c(d$n.contrasts,rev(d$n.contrasts)),c(d$n.genes.upper,rev(d$n.genes.lower)),col = rgb(0.9,0.9,0.9,0.5), border = FALSE)
  lines(d$n.genes.mean, type="l", col=stress.cols[[stress]])
}

par(mar=c(2,1,0,0))
plot(plt.p$n.genes.mean ~ plt.p$p.val, type="n", axes=F, ylim=c(0, 4))
axis(side = 1, las=1, lwd=0, lwd.ticks = 1, labels=F)
axis(side = 2, las=1, lwd=0, lwd.ticks = 1, labels=F)
text(max(plt.p$p.val), 3.6, "Number of selected genes \nby p-value cutoff", pos=2, cex=1.2)
for(stress in c("heat", "drought", "salt")) {
  d = plt.p[plt.p$stress == stress,]
  polygon(c(d$p.val,rev(d$p.val)),c(d$n.genes.upper,rev(d$n.genes.lower)),col = rgb(0.9,0.9,0.9,0.5), border = FALSE)
  lines(d$n.genes.mean ~ d$p.val, type="l", col=stress.cols[[stress]])
}

# % Annotated plots
par(mar=c(2,3,1,0))
plot(plt$tp.mean ~ plt$n.contrasts, type="n", ylim=c(0, 50), axes=F)
axis(side = 1, las=1, lwd=0, lwd.ticks = 1)
axis(side = 2, las=1, lwd=0, lwd.ticks = 1, at=seq(0,50,10), labels=paste0(seq(0,50,10), "%"))
for(stress in c("heat", "drought", "salt")) {
  d = plt[plt$stress == stress,]
  # True positives
  polygon(c(d$n.contrasts,rev(d$n.contrasts)),c(d$tp.upper,rev(d$tp.lower)),col = rgb(0.9,0.9,0.9,0.5), border = FALSE)
  lines(d$tp.mean, type="l", col=stress.cols[[stress]])
  
  # False positives
  polygon(c(d$n.contrasts,rev(d$n.contrasts)),c(d$fp.upper,rev(d$fp.lower)),col = rgb(0.9,0.9,0.9,0.5), border = FALSE)
  lines(d$fp.mean, type="l", col=stress.cols[[stress]], lty=2)
}

par(mar=c(2,1,1,0))
plot(plt.p$tp.mean ~ plt.p$p.val, type="n", ylim=c(0, 50), axes=F)
axis(side = 1, las=1, lwd=0, lwd.ticks = 1, at = seq(0,100,20), labels=paste0("1e-", seq(0,100,20)))
axis(side = 2, las=1, lwd=0, lwd.ticks = 1, labels=F)
#axis(side = 1, las=1, lwd=0, lwd.ticks = 1, at = seq(0,100,20), labels=as.expression(sapply(seq(0,100,20), FUN=function(x) {bquote(1 ^- .(x))} )))
for(stress in c("heat", "drought", "salt")) {
  d = plt.p[plt.p$stress == stress,]
  # True positives
  polygon(c(d$p.val,rev(d$p.val)),c(d$tp.upper,rev(d$tp.lower)),col = rgb(0.9,0.9,0.9,0.5), border = FALSE)
  lines(d$tp.mean ~ d$p.val, type="l", col=stress.cols[[stress]])
  
  # False positives
  polygon(c(d$p.val,rev(d$p.val)),c(d$fp.upper,rev(d$fp.lower)),col = rgb(0.9,0.9,0.9,0.5), border = FALSE)
  lines(d$fp.mean ~ d$p.val, type="l", col=stress.cols[[stress]], lty=2)
}

dev.off()

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