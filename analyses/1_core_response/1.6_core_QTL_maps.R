library(chromoMap)

a = read.csv("input/core_genes.csv", header=F)
q = read.csv("input/qtls.csv", header=F)
# Transform MB positions to b
q$V3 = q$V3 * 1000000
q$V4 = q$V4 * 1000000
a = rbind(a, q)
write.table(a, "input/out.txt", sep="\t", row.names=F, col.names = F, quote=F)

chromoMap("input/chromosomes.txt", "input/out.txt", labels=T, segment_annotation = T,
          win.summary.display=T, n_win.factor = 10, 
          data_based_color_map = T, interactivity = T, label_font = 20, label_angle = -20,
          data_type = "categorical", export.options = T, chr_width = 25, chr_length = 2,
          ch_gap = 2, text_font_size = 15,
          data_colors = list(c("#D55E00", "#0072B2", "#F0E442")))
