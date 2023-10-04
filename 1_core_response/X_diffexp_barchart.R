library(stringr)
# library(ggplot2)

source("../config.R")

diffexp = read.csv("output/differential_expression.lfs.csv.gz")
row.names(diffexp) = diffexp$gene
diffexp = diffexp[-c(1)]

samples.annotation = read.csv("input/samples_annotation.csv")
samples.annotation[samples.annotation$sra_run_id == "",]$sra_run_id = samples.annotation[samples.annotation$sra_run_id == "",]$sample.name
samples.annotation$sample.group = make.names(samples.annotation$sample.group)
samples.annotation = samples.annotation[samples.annotation$sample.group %in% unique(str_split_fixed(names(diffexp), "_", 2)[,1]),]

heat.contrasts = unique(samples.annotation[samples.annotation$stress_type == "heat",]$sample.group)

n.diffexp.genes = as.data.frame(t(sapply(heat.contrasts, FUN = function(h) {
  c(
    up = nrow(diffexp[!is.na(diffexp[[paste0(h,"_FDR")]]) & diffexp[[paste0(h,"_FDR")]] < 0.05 & diffexp[[paste0(h,"_logFC")]] > 0,]),
    down = nrow(diffexp[!is.na(diffexp[[paste0(h,"_FDR")]]) & diffexp[[paste0(h,"_FDR")]] < 0.05 & diffexp[[paste0(h,"_logFC")]] < 0,])
  )
})))
n.diffexp.genes = n.diffexp.genes[rowSums(n.diffexp.genes) > 0,]
n.diffexp.genes$contrast = row.names(n.diffexp.genes)

# n.diffexp.genes.reshaped = reshape(n.diffexp.genes, direction="long", idvar = "contrast", v.names = "Value", varying=c("up", "down"), times=c("up", "down"))
# n.diffexp.genes.reshaped = merge(n.diffexp.genes.reshaped, unique(samples.annotation[c("sample.group", "temperature", "stress.duration")]), by.x="contrast", by.y="sample.group")
# n.diffexp.genes.reshaped = n.diffexp.genes.reshaped[order(n.diffexp.genes.reshaped$temperature, n.diffexp.genes.reshaped$stress.duration, decreasing = TRUE),]
# n.diffexp.genes.reshaped$contrast = factor(n.diffexp.genes.reshaped$contrast, levels=unique(n.diffexp.genes.reshaped$contrast))
# 
# ggplot(n.diffexp.genes.reshaped, aes(fill=time, y=Value, x=contrast)) + 
#   geom_bar(position="stack", stat="identity") +
#   theme_minimal()

n.diffexp.genes = merge(n.diffexp.genes, unique(samples.annotation[c("sample.group", "temperature", "stress.duration", "tissue")]), by.x="contrast", by.y="sample.group")
n.diffexp.genes = n.diffexp.genes[order(n.diffexp.genes$temperature, n.diffexp.genes$stress.duration, decreasing = TRUE),]
n.diffexp.genes$col = colors[["tissue"]][match(n.diffexp.genes$tissue, names(colors[["tissue"]]))]

color_segments <- function(df, x.coords, print.text = F) {
  df[is.na(df$temperature),]$temperature <- "?"
  current.tmp = df[1,]$temperature
  x.start <- 0
  for(row in 1:nrow(df)) {
    print(df[row,]$contrast)
    print(paste(df[row,]$temperature, current.tmp))
    if(df[row,]$temperature != current.tmp) {
      if(print.text) {
        mtext(paste0(current.tmp, "°C"), at = mean(c(x.start, x.coords[row])), adj=1, cex=0.3)
      }
      current.tmp = df[row,]$temperature
      abline(v = x.coords[row] - 0.6)
      x.start = x.coords[row]
    }
  }
}

pdf("output/plots/1.X_barplots.pdf", 4, 3)
par(mfrow=c(2,1), mar=c(0.5,4,2,1), cex=0.4)
bp.x.coords <- barplot(n.diffexp.genes$up, axes=F, col=n.diffexp.genes$col, border = "white")
color_segments(n.diffexp.genes, bp.x.coords, T)
axis(1, at=bp.x.coords, labels = n.diffexp.genes$contrast, las=2, lwd = 0)
axis(2, las=2)
#polygon(x=c(0,0,bp.x.coords[7] + 0.6, bp.x.coords[7] + 0.6), y=c(-7000,max(n.diffexp.genes$up),max(n.diffexp.genes$up),-7000), col="#0000FF22", border=F, xpd=T)
#polygon(x=c(bp.x.coords[11] - 0.6, bp.x.coords[11] - 0.6, bp.x.coords[12] + 0.6, bp.x.coords[12] + 0.6), y=c(-7000,max(n.diffexp.genes$up),max(n.diffexp.genes$up),-7000), col="#0000FF22", border=F, xpd=T)
par(mar=c(2,4,0,1))
barplot(-n.diffexp.genes$down, axes=F, col=n.diffexp.genes$col, border="white")
color_segments(n.diffexp.genes, bp.x.coords)
#polygon(x=c(0,0,bp.x.coords[7] + 0.6, bp.x.coords[7] + 0.6), y=c(-7000,max(n.diffexp.genes$up),max(n.diffexp.genes$up),-7000), col="#0000FF22", border=F, xpd=T)
axis(2, las=2)
dev.off()

