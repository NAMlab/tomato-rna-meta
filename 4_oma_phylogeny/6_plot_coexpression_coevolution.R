library(ComplexHeatmap)

evo = read.csv("output/heat_hogprof_jaccard_matrix.csv")
abundance = read.csv("../data/combined_abundance.lfs.tsv.gz", sep="\t", check.names=F)
internal.names = read.csv("../1_core_response/input/hs_core_genes_internal_names.csv")

abundance = abundance[abundance$target_id %in% evo$target,]
abundance = merge(abundance, internal.names, by.x = "target_id", by.y="target")
row.names(abundance) = abundance$internal.name
abundance = abundance[grep("_tpm", names(abundance))]
# Uncomment these lines if you want to base co-expression only on heat and control samples
samples.annotation = read.csv("../data/samples_annotation.csv")
samples.annotation[samples.annotation$sra_run_id == "",]$sra_run_id = samples.annotation[samples.annotation$sra_run_id == "",]$sample.name
samples.annotation = samples.annotation[samples.annotation$stress_type %in% c("heat", "none"),]
abundance = abundance[grep(paste(samples.annotation$sra_run_id, collapse="|"), names(abundance))]

evo = merge(evo, internal.names, sort=F)
row.names(evo) = evo$internal.name
evo = evo[2:(ncol(evo)-1)]
names(evo) = row.names(evo)

coexp = cor(t(abundance))

# Sort both matrices alphabetically so the row and column order is the same (and thus transferrable)
coexp = coexp[order(row.names(coexp)),order(colnames(coexp))]
evo = evo[order(row.names(evo)),order(names(evo))]


pdf("output/plots/4.3_coexpression_coevolution.pdf", 15, 15)

ht.coexp = draw(Heatmap(coexp))
ht.evo = draw(Heatmap(evo))

ht.evo.ordered = draw(Heatmap(evo, row_order = row_order(ht.coexp), column_order = column_order(ht.coexp)))
ht.coexp.ordered = draw(Heatmap(coexp, row_order = row_order(ht.evo), column_order = column_order(ht.evo)))
dev.off()

# Do plotting semi-manually to achieve triangle matrix
library(circlize)
col.coexp = colorRamp2(c(-1, 0, 1), c("purple", "white", "orange"))
pdf("output/plots/4.3_coexpression_coevolution_triangle.pdf", 6, 6)
h1 <- draw(Heatmap(coexp, rect_gp = gpar(type = "none"), show_column_names = F, show_row_dend = F,
             row_names_gp = gpar(fontsize = 7), name="Co-Expression (TPM Correlation)", col=col.coexp,
             heatmap_legend_param = list(legend_direction = "horizontal", legend_width = unit(6, "cm")),
             cell_fun = function(j, i, x, y, w, h, fill) {
               if(i == j) {
                 grid.text(colnames(coexp)[i], x, y, just="right", gp = gpar(fontsize=7))
               } else if(as.numeric(x) >= 1 - as.numeric(y) + 1e-6) {
                 grid.rect(x, y, w, h, gp = gpar(fill = fill, col = fill))
               }
             }), heatmap_legend_side="top", padding = unit(c(0.1, 1.2, 0.1, 0.1), "cm"), background = "transparent")

# @TODO Wouldn't it make a lot more sense to just have one quadratic heatmap with the clustering of Coexp,
# the upper triangle with coexp values (i.e. what the clustering is based on) and the lower triangle with
# co-evo values? Diagonal black or sth to separate the two
# GAWD. That would be so much more intuitive than all the triangles I just did lol. Although the are pretty!
# Maybe I can somehow use the leftover triangle...

h2 <- draw(Heatmap(evo, rect_gp = gpar(type = "none"), show_row_names = F, show_column_dend = F,
        column_names_gp = gpar(fontsize = 7), name="Co-Evolution (Jaccard Similarity)",
        heatmap_legend_param = list(legend_direction = "horizontal", legend_width = unit(6, "cm")),
        cell_fun = function(j, i, x, y, w, h, fill) {
          if(i == j) {
            grid.text(names(evo)[i], x, y, just="left", gp = gpar(fontsize=7))
          } else if(as.numeric(x) <= 1 - as.numeric(y) + 1e-6) {
            grid.rect(x, y, w, h, gp = gpar(fill = fill, col = fill))
          }
        }), heatmap_legend_side="bottom", padding = unit(c(0.1, 0.1, 0.1, 1.2), "cm"), background = "transparent")

dev.off()

