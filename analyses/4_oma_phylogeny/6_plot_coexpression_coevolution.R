library(ComplexHeatmap)

evo = read.csv("output/heat_hogprof_jaccard_matrix.csv")
evo.all = read.csv("output/heat_hogprof_jaccard_matrix_all_taxa.csv")
internal.names = read.csv("../1_core_response/input/hs_core_genes_internal_names.csv")

evo = merge(evo, internal.names, sort=F)
row.names(evo) = evo$internal.name
evo = evo[2:(ncol(evo)-1)]
names(evo) = row.names(evo)

evo.all = merge(evo.all, internal.names, sort=F)
row.names(evo.all) = evo.all$internal.name
evo.all = evo.all[2:(ncol(evo.all)-1)]
names(evo.all) = row.names(evo.all)

# Sort both matrices alphabetically so the row and column order is the same (and thus transferrable)
evo = evo[order(row.names(evo)),order(names(evo))]
evo.all = evo.all[order(row.names(evo.all)),order(names(evo.all))]


pdf("output/plots/4.6_coexpression.pdf", 15, 15)

ht.evo = draw(Heatmap(evo))
ht.evo.all = draw(Heatmap(evo.all))

ht.evo.ordered = draw(Heatmap(evo, row_order = row_order(ht.evo.all), column_order = column_order(ht.evo.all)))
ht.evo.all.ordered = draw(Heatmap(evo.all, row_order = row_order(ht.evo), column_order = column_order(ht.evo)))
dev.off()


# Read MRCAs so we can add them to the heatmap
d = read.csv("output/combined_deepest_levels.csv")
# Prepare for plotting the genes at each level
d.h = d[d$set == "heat" & d$level < 16,]
d.h$level = d.h$level + 1
d.h = merge(d.h, read.csv("../1_core_response/input/hs_core_genes_internal_names.csv"))
# Order them for the heatmap
d.h2 <- d.h[match(names(evo), d.h$internal.name), ]

# Do plotting semi-manually to achieve triangle matrix
library(circlize)
col.evo = colorRamp2(c(0, 0.5, 1), c("blue", "white", "red"))

pdf("output/plots/4.6_coevolution_triangle.pdf", 6.5, 6.95)
column_ha = HeatmapAnnotation(MRCA = d.h2$level, annotation_legend_param = list(
  title = "MRCA", at = c(1, 15), 
  labels = c("LUCA", "Solanoideae"),
  direction = "horizontal"
))

h3 <- draw(Heatmap(evo, rect_gp = gpar(type = "none"), show_row_dend = F,
                   row_names_gp = gpar(fontsize = 7), column_names_gp = gpar(fontsize = 7), name="Co-Evolution across Streptophyta", col=col.evo,
                   heatmap_legend_param = list(legend_direction = "horizontal", legend_width = unit(6, "cm")),
                   bottom_annotation = column_ha,
                   cell_fun = function(j, i, x, y, w, h, fill) {
                     if(i == j) {
                       grid.rect(x, y, w, h, gp = gpar(fill = "gray40", col = "gray50"))
                     } else if(as.numeric(x) >= 1 - as.numeric(y) + 1e-6) {
                       grid.rect(x, y, w, h, gp = gpar(fill = fill, col = fill))
                     } else {
                       grid.rect(x, y, w, h, gp = gpar(fill = col.evo(evo.all[i,j]), col = col.evo(evo.all[i,j])))
                     }
                   }), heatmap_legend_side="top", annotation_legend_side="top", padding = unit(c(0.1, 1.2, 0.1, 0.1), "cm"), background = "transparent")
dev.off()

