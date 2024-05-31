# This script generates (for each stress type) a heatmap visualizing all genes that are differentially expressed
# in *any* sample (that is, the union of all of them). To ensure we don't just plot the whole genome,
# we are applying a somewhat stricter p-value threshold: Instead of just adjusting the p-value within each
# experiment separately, we use the raw p values and then adjust them across *all* the contrasts.
# Only genes differentially expressed in any contrast using that p value are included in the heatmap.


library(stringr)
library(ggplot2)
library(ggrepel)
source("lib_GO_enrichment.R")
diffexp = read.csv("output/differential_expression.lfs.csv.gz")
row.names(diffexp) = diffexp$gene

abundance = read.csv("input/combined_abundance.tsv.gz", sep="\t")
row.names(abundance) = abundance$target_id

samples.annotation = read.csv("input/samples_annotation.csv")
samples.annotation[samples.annotation$sra_run_id == "",]$sra_run_id = samples.annotation[samples.annotation$sra_run_id == "",]$sample.name
samples.annotation$sample.group = make.names(samples.annotation$sample.group)
samples.annotation$respective_control = make.names(samples.annotation$respective_control)
samples.annotation$sra_run_id = make.names(samples.annotation$sra_run_id)

protein.descriptions = read.csv("input/protein_descriptions.csv.gz")
id.mapping = read.csv("../data/id_mapping.csv.gz")

for(trt.group in unique(samples.annotation$sample.group)) {
        print(trt.group)
        ctrl.group = unique(samples.annotation[samples.annotation$sample.group == trt.group,]$respective_control)
        if(length(ctrl.group) > 1)
                stop(paste0("Multiple different control groups for sample group ", trt.group))
        if(ctrl.group == "" || ctrl.group == "X") {
                print("This is a control group (skip)")
                next
        }
        dir.create(paste0("output/individual_contrasts_supplement/", trt.group))
        samples = samples.annotation[samples.annotation$sample.group %in% c(ctrl.group, trt.group),]$sra_run_id
        
        #### TPMs, logFC, PValue, protein info #####
        out.csv = abundance[grep(paste0(samples, collapse="|"), names(abundance))]
        out.csv$target_id = row.names(out.csv)
        out.csv$Gramene.ID = str_remove(out.csv$target_id, "\\.[0-9]+$")
        out.csv$gene = str_remove(out.csv$target_id, "\\.[0-9]+\\.[0-9]+$")
        out.csv = merge(id.mapping, out.csv, by="Gramene.ID", all.y=T)
        out.csv = merge(protein.descriptions[c("gene", "ITAG4.1_description", "OMA_orthologues")], out.csv, by="gene")
        
        out.csv = merge(out.csv, diffexp[c(1, grep(trt.group, names(diffexp)))], by.x="target_id", by.y="gene")
        # Remove Gramene.ID and gene
        out.csv = out.csv[-c(2,5)]
        names(out.csv)[(ncol(out.csv)-2):ncol(out.csv)] = c("logFC", "PValue", "FDR")
        
        write.csv(out.csv, paste0("output/individual_contrasts_supplement/", trt.group, "/gene_expression.csv"), row.names=F, na = "")
        
        if(sum(is.na(out.csv$PValue)) > 0) {
                print("No replicates/P-Value (skip)")
                next
        }
        
        #### GO terms ####
        genes.up = unique(str_remove(out.csv[out.csv$FDR < 0.05 & out.csv$logFC > 0,]$target_id, "\\.[0-9]+\\.[0-9]+$"))
        write.csv(go_enrichment(genes.up), paste0("output/individual_contrasts_supplement/", trt.group, "/GO_upregulated.csv"), row.names=F)
        genes.down = unique(str_remove(out.csv[out.csv$FDR < 0.05 & out.csv$logFC < 0,]$target_id, "\\.[0-9]+\\.[0-9]+$"))
        write.csv(go_enrichment(genes.down), paste0("output/individual_contrasts_supplement/", trt.group, "/GO_downregulated.csv"), row.names=F)
        
        ggplot(data=out.csv, aes(x=logFC, y=-log10(FDR)))+
                geom_point(aes(color=ifelse(FDR < 0.05, ifelse(logFC < 0, "down", "up"), "ns"), alpha=0.5))+
                scale_color_manual(values=c("up"= "#D55E00",
                                            "down" = "#0072B2",
                                            "ns"="black"))+
                geom_text_repel(label=ifelse(-log10(out.csv$FDR)>10,out.csv$target_id,""))+
                theme_minimal()+
                theme(legend.position="none")
        ggsave(paste0("output/individual_contrasts_supplement/", trt.group, "/volcano_plot.jpg"), width = 30, height = 20, units = "cm")
}

# @TODO compress each folder and put into lfs?