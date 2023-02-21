library(ComplexHeatmap)

evo = read.csv("output/heat_hogprof_jaccard_matrix.csv")
abundance = read.csv("../data/combined_abundance.lfs.tsv.gz", sep="\t", check.names=F)

samples.annotation = read.csv("../data/samples_annotation.csv")
samples.annotation[samples.annotation$sra_run_id == "",]$sra_run_id = samples.annotation[samples.annotation$sra_run_id == "",]$sample.name

abundance = abundance[abundance$target_id %in% evo$target,]
row.names(abundance) = abundance$target_id
abundance = abundance[grep("_tpm", names(abundance))]
# Uncomment these lines if you want to base co-expression only on heat and control samples
# samples.annotation = samples.annotation[samples.annotation$stress_type %in% c("heat", "none"),]
# abundance = abundance[grep(paste(samples.annotation$sra_run_id, collapse="|"), names(abundance))]

row.names(evo) = evo$target
evo = evo[2:ncol(evo)]

coexp = cor(t(abundance))

pdf("output/plots/4.3_coexpression_coevolution.pdf")
ht = draw(Heatmap(coexp))

draw(Heatmap(evo))

ht.evo.ordered = Heatmap(evo, row_order = row_order(ht), column_order = column_order(ht))
draw(ht.evo.ordered)
dev.off()
