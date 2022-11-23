library(stringr)
library(Rmisc)
target.genes = scan("output/hs_genes.txt", character())
pollen.genes = scan("output/pollen_genes.txt", character())

diffexp = read.csv("../1_core_response/output/differential_expression.lfs.csv.gz")
samples.annotation = read.csv("../1_core_response/input/samples_annotation.csv")
samples.annotation[samples.annotation$sra_run_id == "",]$sra_run_id = samples.annotation[samples.annotation$sra_run_id == "",]$sample.name
samples.annotation$sample.group = make.names(samples.annotation$sample.group)
# Only keep HS lines
samples.annotation = samples.annotation[samples.annotation$stress_type == "heat",]

# Remove contrasts without replicates
b = table(samples.annotation$sample.group)
samples.annotation = samples.annotation[!samples.annotation$sample.group %in% names(b[b < 2]),]

# Remove all samples from diffexp that are not part of our samples.annotation
diffexp = diffexp[grep(paste0(c("gene", unique(samples.annotation$sample.group)), collapse="|"), names(diffexp))]
row.names(diffexp) = diffexp$gene

contrasts = unique(samples.annotation$sample.group)

# For now, only one run-through:
res = data.frame(run = numeric(), n.contrasts = numeric(), n.genes = numeric(), true.pos = numeric(), pollen = numeric())
for(run in 1:200) {
  contrasts.queue = sample(contrasts)
  for(n.contrasts in 1:length(contrasts.queue)) {
    b = rowSums(diffexp[paste0(contrasts.queue[1:n.contrasts], "_FDR")] < 0.05, na.rm = T)
    found.genes = str_split_fixed(names(b[b >= (0.8 * n.contrasts)]), "\\.", 2)[,1]
    #found.genes = str_split_fixed(names(b[b >= 1]), "\\.", 2)[,1]
    true.pos = sum(found.genes %in% target.genes)
    pollen = sum(found.genes %in% pollen.genes)
    res = rbind(res, list(run = run, n.contrasts = n.contrasts, n.genes = length(found.genes), true.pos = true.pos, pollen=pollen))
  }
}
# Number of genes in intersect set
ci.genes = group.CI(n.genes ~ n.contrasts, res)

plot(ci.genes$n.genes.mean ~ ci.genes$n.contrasts, type="n", xlab="# of contrasts", las=2, ylab="# of \"core\" genes")
polygon(c(ci.genes$n.contrasts,rev(ci.genes$n.contrasts)),c(ci.genes$n.genes.upper,rev(ci.genes$n.genes.lower)),col = "grey90", border = FALSE)
lines(ci.genes$n.genes.mean, type="l")


# Number of Genes annotated with a specific term
res$prop = (res$true.pos / (res$n.genes)) * 100 # [%]
ci.prop = group.CI(prop ~ n.contrasts, res)

res$pollen.prop = (res$pollen / (res$n.genes)) * 100 # [%]
ci.pollen = group.CI(pollen ~ n.contrasts, res)

plot(ci.prop$prop.mean ~ ci.prop$n.contrasts, type="n", xlab="# of contrasts", las=2, ylab="Genes Annotated [%]", ylim=c(0,18))
polygon(c(ci.prop$n.contrasts,rev(ci.prop$n.contrasts)),c(ci.prop$prop.upper,rev(ci.prop$prop.lower)),col = "grey90", border = FALSE)
lines(ci.prop$prop.mean, type="l")
polygon(c(ci.pollen$n.contrasts,rev(ci.pollen$n.contrasts)),c(ci.pollen$pollen.upper,rev(ci.pollen$pollen.lower)),col = "grey90", border = FALSE)
lines(ci.pollen$pollen.mean, type="l", col="red")

