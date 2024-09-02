heat.genes = read.csv("../analyses/1_core_response/output/1.4_core_response_genes.csv")
heat.genes = heat.genes[heat.genes$stress.type == "heat",]
heat.genes.up = heat.genes[heat.genes$direction == "upregulated",]

family.map = read.csv("output/gene_family_annotations.csv.gz")
m = merge(heat.genes.up, family.map, by.x="target.id", by.y="target_id")

library(stringr)

diffexp = read.csv("../analyses/1_core_response/output/differential_expression.lfs.csv.gz")
row.names(diffexp) = diffexp$gene
diffexp = diffexp[-c(1)]

samples.annotation = read.csv("../analyses/1_core_response/input/samples_annotation.csv")
samples.annotation[samples.annotation$sra_run_id == "",]$sra_run_id = samples.annotation[samples.annotation$sra_run_id == "",]$sample.name
samples.annotation$sample.group = make.names(samples.annotation$sample.group)
samples.annotation = samples.annotation[samples.annotation$sample.group %in% unique(str_split_fixed(names(diffexp), "_", 2)[,1]),]
#internal.names = read.csv("../1_core_response/input/hs_core_genes_internal_names.csv")

up.genes = data.frame(sample = character(), target_id = character())

stress.type = "heat"
samples.stress = unique(samples.annotation[samples.annotation$stress_type == stress.type,c("sample.group", "tissue", "temperature", "stress.duration", "genotype.name")])
diffexp.stress = diffexp[grep(paste0(c("gene", samples.stress$sample.group), collapse="|"), names(diffexp))]

samples = unique(samples.stress$sample.group)
samples.without.p.val = is.na(diffexp.stress[1,paste0(samples, "_FDR")])
samples.without.p.val = str_remove(colnames(samples.without.p.val)[which(samples.without.p.val)], "_FDR")
samples = setdiff(samples, samples.without.p.val)

for(s in samples) {
  genes = row.names(diffexp.stress[diffexp.stress[paste0(s, "_FDR")] < 0.05  & diffexp.stress[paste0(s, "_logFC")] > 0,])
  up.genes = rbind(up.genes, data.frame(sample = s, target_id = genes))
}

up.genes = merge(up.genes, family.map)

# Find out which are the core families
core.fams = unique(up.genes[c("sample", "family_id")])
b = table(core.fams$family_id)
core.fams = names(b[b >= length(samples) * 0.9])

up.genes = up.genes[up.genes$family_id %in% core.fams,]

family.info = read.csv("output/family_info.csv.gz")               
up.genes = merge(up.genes, family.info)

# Now what does subfunctionalization mean?
