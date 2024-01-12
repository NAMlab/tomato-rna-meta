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

# Now again but using the score of https://doi.org/10.1093/bioinformatics/btr671 
plt.s = data.frame(s=numeric(), n.genes.upper=numeric(), n.genes.mean=numeric(), n.genes.lower=numeric(), 
                   tp.upper=numeric(), tp.mean=numeric(), tp.lower=numeric(), 
                   fp.upper=numeric(), fp.mean=numeric(), fp.lower=numeric(),
                   stress=character())
for(stress in c("heat", "drought", "salt")) {
  samples.stress = samples.annotation[samples.annotation$stress_type == stress,]
  
  contrasts = unique(samples.stress$sample.group)
  
  res = data.frame(s = numeric(), n.genes=numeric(), true.pos=numeric(), false.pos=numeric())
  for(s in c(0, 1, 2, 3, 5, 10, 25, 100, 200, 300, 400, 500)) {
    print(s)
    for(co in contrasts) {
      found.genes = str_split_fixed(
        # They abs() their logFC but since we only care about upregulated genes we can skip that and then
        # any genes where the overall π score is negative are filtered out as well.
        diffexp[-log(diffexp[paste0(co, "_FDR")], 10) * diffexp[paste0(co, "_logFC")] > s,]$gene,"\\.", 2)[,1]
      true.pos = sum(found.genes %in% annotated.genes[annotated.genes$stress == stress & annotated.genes$reasonable_term == "true",]$gene)
      false.pos = sum(found.genes %in% annotated.genes[annotated.genes$stress == stress & annotated.genes$reasonable_term == "false",]$gene)
      res = rbind(res, list(s=s, n.genes=length(found.genes), true.pos=true.pos, false.pos=false.pos))
    }
  }
  
  # Number of genes in intersect set
  ci.genes = group.CI(n.genes ~ s, res)
  
  # Number of Genes annotated with a specific term
  res$tp = (res$true.pos / (res$n.genes)) * 100 # [%]
  ci.tp = group.CI(tp ~ s, res)
  
  res$fp = (res$false.pos / (res$n.genes)) * 100 # [%]
  ci.fp = group.CI(fp ~ s, res)
  
  
  m = merge(ci.genes, merge(ci.tp, ci.fp))
  m$stress = stress
  plt.s = rbind(plt.s, m)
}

plt.s$n.genes.lower = log(plt.s$n.genes.lower, base=10)
plt.s$n.genes.mean = log(plt.s$n.genes.mean, base=10)
plt.s$n.genes.upper = log(plt.s$n.genes.upper, base=10)

##### Plots #####
source("../config.R")
stress.cols = colors[["stress"]]

pdf("output/curves.pdf", 6.5, 5)
par(mfrow=c(2,3), family="serif", cex=0.6)

# N Genes Plots
par(mar=c(2,4,0,0))
plot(plt.p$n.genes.mean ~ plt.p$p.val, type="n", axes=F, ylim=c(0, 4), ylab="")
axis(side = 1, las=1, col=F, col.axis="gray52", at = seq(0,100,25), labels=paste0("1e-", seq(0,100,25)))
axis(side = 2, las=1, at = 0:4, labels=c(1,10,100,"1k","10k"), col=F, col.axis="gray52")
#text(max(plt.p$p.val), 3.6, "number of selected genes \nby p-value cutoff", pos=2, cex=1.1)
for(stress in c("heat", "drought", "salt")) {
  d = plt.p[plt.p$stress == stress,]
  polygon(c(d$p.val,rev(d$p.val)),c(d$n.genes.upper,rev(d$n.genes.lower)), col = rgb(t(col2rgb(stress.cols[[stress]])/255), alpha = 0.1), border = FALSE)
  lines(d$n.genes.mean ~ d$p.val, type="l", col=stress.cols[[stress]])
  mtext("|", side=2, at=min(d$n.genes.mean), col=stress.cols[[stress]], cex=0.5)
}
mtext("selected genes", side=2, line=2.6, col="gray52", cex=0.6)

par(mar=c(2,2.5,0,0))
plot(plt.s$n.genes.mean ~ plt.s$s, type="n", axes=F, ylim=c(0, 4), ylab="")
axis(side = 1, las=1, col=F, col.axis="gray52", at = seq(0,500,100), labels=seq(0,500,100))
axis(side = 2, las=1, at = 0:4, labels=c(1,10,100,"1k","10k"), col=F, col.axis="gray52")
#text(max(plt.p$p.val), 3.6, "number of selected genes \nby π score cutoff", pos=2, cex=1.1)
for(stress in c("heat", "drought", "salt")) {
  d = plt.s[plt.s$stress == stress,]
  polygon(c(d$s,rev(d$s)),c(d$n.genes.upper,rev(d$n.genes.lower)), col = rgb(t(col2rgb(stress.cols[[stress]])/255), alpha = 0.1), border = FALSE)
  lines(d$n.genes.mean ~ d$s, type="l", col=stress.cols[[stress]])
  mtext("|", side=2, at=min(d$n.genes.mean), col=stress.cols[[stress]], cex=0.5)
}

plot(plt$n.genes.mean ~ plt$n.contrasts, type="n", ylim=c(0, 4), axes=F, ylab="")
axis(side = 1, las=1, col=F, col.axis="gray52")
axis(side = 2, las=1, at = 0:4, labels=c(1,10,100,"1k","10k"), col=F, col.axis="gray52")
#text(max(plt$n.contrasts), 3.6, "number of selected genes \nby number of contrasts (50 simulations)", pos=2, cex=1.1)
for(stress in c("heat", "drought", "salt")) {
  d = plt[plt$stress == stress,]
  polygon(c(d$n.contrasts,rev(d$n.contrasts)),c(d$n.genes.upper,rev(d$n.genes.lower)), col = rgb(t(col2rgb(stress.cols[[stress]])/255), alpha = 0.1), border = FALSE)
  lines(d$n.genes.mean, type="l", col=stress.cols[[stress]])
  mtext("|", side=2, at=min(d$n.genes.mean), col=stress.cols[[stress]], cex=0.5)
}
text(12.5, log(240, 10), "salt", col=stress.cols[["salt"]])
text(20, log(20, 10), "drought", col=stress.cols[["drought"]])
text(25.3, log(100, 10), "heat", col=stress.cols[["heat"]])

# % Annotated plots
par(mar=c(3.5,4,1,0))
plot(plt.p$tp.mean ~ plt.p$p.val, type="n", ylim=c(0, 50), axes=F, xlab="", ylab="")
axis(side = 1, las=1, col=F, col.axis="gray52", at = seq(0,100,25), labels=paste0("1e-", seq(0,100,25)))
axis(side = 2, las=1, col=F, col.axis="gray52", at=seq(0,50,10), labels=paste0(seq(0,50,10), "%"))
#axis(side = 1, las=1, lwd=0, lwd.ticks = 1, at = seq(0,100,20), labels=as.expression(sapply(seq(0,100,20), FUN=function(x) {bquote(1 ^- .(x))} )))
for(stress in c("heat", "drought", "salt")) {
  d = plt.p[plt.p$stress == stress,]
  # True positives
  polygon(c(d$p.val,rev(d$p.val)),c(d$tp.upper,rev(d$tp.lower)), col = rgb(t(col2rgb(stress.cols[[stress]])/255), alpha = 0.1), border = FALSE)
  lines(d$tp.mean ~ d$p.val, type="l", col=stress.cols[[stress]])
  mtext("|", side=2, at=max(d$tp.mean), col=stress.cols[[stress]], cex=0.5)
  
  # False positives
  polygon(c(d$p.val,rev(d$p.val)),c(d$fp.upper,rev(d$fp.lower)), col = rgb(t(col2rgb(stress.cols[[stress]])/255), alpha = 0.1), border = FALSE)
  lines(d$fp.mean ~ d$p.val, type="l", col=stress.cols[[stress]], lty=2)
}
#text(0.5, 46, "proportion of correctly/incorrectly annotated genes \nby p-value cutoff", pos=4, cex=1.1)
#segments(39, 47, 55, lty=2)
mtext("p-value cutoff", side=1, line=2, col="gray52", cex=0.6)
mtext("(in-)correctly annotated genes", side=2, line=2.6, col="gray52", cex=0.6)

par(mar=c(3.5,2.5,1,0))
plot(plt.s$tp.mean ~ plt.s$s, type="n", ylim=c(0, 50), axes=F, xlab="", ylab="")
axis(side = 1, las=1, col=F, col.axis="gray52", at = seq(0,500,100), labels=seq(0,500,100))
axis(side = 2, las=1, col=F, col.axis="gray52", at=seq(0,50,10), labels=paste0(seq(0,50,10), "%"))
#axis(side = 1, las=1, lwd=0, lwd.ticks = 1, at = seq(0,100,20), labels=as.expression(sapply(seq(0,100,20), FUN=function(x) {bquote(1 ^- .(x))} )))
for(stress in c("heat", "drought", "salt")) {
  d = plt.s[plt.s$stress == stress,]
  # True positives
  polygon(c(d$s,rev(d$s)),c(d$tp.upper,rev(d$tp.lower)), col = rgb(t(col2rgb(stress.cols[[stress]])/255), alpha = 0.1), border = FALSE)
  lines(d$tp.mean ~ d$s, type="l", col=stress.cols[[stress]])
  mtext("|", side=2, at=max(d$tp.mean), col=stress.cols[[stress]], cex=0.5)
  
  # False positives
  polygon(c(d$s,rev(d$s)),c(d$fp.upper,rev(d$fp.lower)), col = rgb(t(col2rgb(stress.cols[[stress]])/255), alpha = 0.1), border = FALSE)
  lines(d$fp.mean ~ d$s, type="l", col=stress.cols[[stress]], lty=2)
}
#text(0.5, 46, "proportion of correctly/incorrectly annotated genes \nby p-value cutoff", pos=4, cex=1.1)
#segments(39, 47, 55, lty=2)
mtext(expression(pi * "-score cutoff"), side=1, line=2, col="gray52", cex=0.6)

plot(plt$tp.mean ~ plt$n.contrasts, type="n", ylim=c(0, 50), axes=F, xlab="", ylab="")
axis(side = 1, las=1, col=F, col.axis="gray52")
axis(side = 2, las=1, col=F, col.axis="gray52", at=seq(0,50,10), labels=paste0(seq(0,50,10), "%"))
#text(0.5, 46, "proportion of correctly/incorrectly annotated genes \nby number of contrasts (50 simulations)", pos=4, cex=1.1)
mtext("number of contrasts (50 simulations)", side=1, line=2, col="gray52", cex=0.6)
#segments(11.8, 47, 16.8, lty=2)
for(stress in c("heat", "drought", "salt")) {
  d = plt[plt$stress == stress,]
  # True positives
  polygon(c(d$n.contrasts,rev(d$n.contrasts)),c(d$tp.upper,rev(d$tp.lower)), col = rgb(t(col2rgb(stress.cols[[stress]])/255), alpha = 0.1), border = FALSE)
  lines(d$tp.mean, type="l", col=stress.cols[[stress]])
  mtext("|", side=2, at=max(d$tp.mean), col=stress.cols[[stress]], cex=0.5)
  
  # False positives
  polygon(c(d$n.contrasts,rev(d$n.contrasts)),c(d$fp.upper,rev(d$fp.lower)), col = rgb(t(col2rgb(stress.cols[[stress]])/255), alpha = 0.1), border = FALSE)
  lines(d$fp.mean, type="l", col=stress.cols[[stress]], lty=2)
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