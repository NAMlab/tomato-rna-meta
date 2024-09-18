library(stringr)
# Visualize the number of differentially expressed genes in heat stress by temperature and duration.
# Also show the number of significantly enriched GO terms.

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

# Add info about number of significantly enriched GO terms
GO.info = as.data.frame(t(sapply(n.diffexp.genes$contrast, FUN=function(a) {
  up = read.csv(paste0("output/individual_contrasts_supplement/",a,"/GO_upregulated.csv"))
  down = read.csv(paste0("output/individual_contrasts_supplement/",a,"/GO_downregulated.csv"))
  c(GO.up = nrow(up[up$q_value <= 0.05,]), GO.down = nrow(down[down$q_value <= 0.05,]))
})))
n.diffexp.genes = cbind(n.diffexp.genes, GO.info)

color_segments <- function(df, x.coords, print.text = F) {
  df[is.na(df$temperature),]$temperature <- "?"
  current.tmp = df[1,]$temperature
  x.start <- 0
  for(row in 1:nrow(df)) {
    if(df[row,]$temperature != current.tmp) {
      if(print.text) {
        mtext(paste0(current.tmp, "°C"), at = mean(c(x.start, x.coords[row])), adj=1, cex=0.4)
      }
      current.tmp = df[row,]$temperature
      abline(v = x.coords[row] - 0.6)
      x.start = x.coords[row]
    }
  }
}
# n.diffexp.genes[is.na(n.diffexp.genes$stress.duration),]$stress.duration = "?"

used.tissues = colors[["tissue"]][names(colors[["tissue"]]) %in% n.diffexp.genes$tissue]

cairo_pdf("output/plots/1.4_barplots.pdf", 5, 6) # to get short hyphens in the axis labels
par(mfrow=c(4,1), mar=c(1,4,2,1), cex=0.4, family="serif")
bp.x.coords <- barplot(n.diffexp.genes$up, axes=F, col=n.diffexp.genes$col, border = "white", ylab = "upregulated", xlab="stress duration[h]")
color_segments(n.diffexp.genes, bp.x.coords, T)
axis(2, las=2, lwd=0, line=-1.5)
axis(1, las=2, lwd=0, labels = n.diffexp.genes$stress.duration, at=bp.x.coords, line=-0.75)
mtext("stress duration [h]", side=1, at=max(bp.x.coords), adj=1, xpd=T, cex=0.4, line=0.5)
mtext("DEGs", side=2, at = -500, xpd=T, cex=0.7, line=2)
legend("topleft", inset=.03, title="Tissue",
       names(used.tissues), border=used.tissues, fill=used.tissues, box.lty=0, cex=1.2)
par(mar=c(2,4,1,1))
barplot(-n.diffexp.genes$down, axes=F, col=n.diffexp.genes$col, border="white", ylab = "downregulated")
axis(1, at=bp.x.coords, labels=str_replace_all(n.diffexp.genes$contrast, "\\.", "-"), las=2, lwd = 0)
color_segments(n.diffexp.genes, bp.x.coords)
axis(2, las=2, lwd=0, line=-1.5, at=seq(-0,-7000,-1000), labels=seq(0, 7000, 1000))

par(mar=c(1,4,0.5,1))
barplot(n.diffexp.genes$GO.up, axes=F, col=n.diffexp.genes$col, border = "white", ylim=c(0, 1800), ylab = "upregulated")
color_segments(n.diffexp.genes, bp.x.coords)
axis(2, las=2, lwd=0, line=-1.5)
axis(1, las=2, lwd=0, labels = n.diffexp.genes$stress.duration, at=bp.x.coords, line=-0.75)
mtext("stress duration [h]", side=1, at=max(bp.x.coords), adj=1, xpd=T, cex=0.4, line=0.5)
mtext("GO Terms", side=2, at = -50, xpd=T, cex=0.7, line=2)
par(mar=c(2,4,1,1))
barplot(-n.diffexp.genes$GO.down, axes=F, col=n.diffexp.genes$col, border="white", ylim=c(-1700, 0), ylab = "downregulated")
color_segments(n.diffexp.genes, bp.x.coords)
axis(2, las=2, lwd=0, line=-1.5, at=seq(-0,-1500,-500), labels=seq(0, 1500, 500))
dev.off()

# Try a different method of visualizing the data.
pdf("output/plots/1.4_non_barplots.pdf", 4, 3)
par(cex = 0.5, mar=c(0.5,4,2,1), family="serif")
plot(n.diffexp.genes$up, col=n.diffexp.genes$col, pch=0, ylab = "differentially expressed genes", cex=2, axes=F)
points(n.diffexp.genes$down, col=n.diffexp.genes$col, pch=1, cex=2)
axis(2, las=2, lwd=0, line=0, lwd.ticks = 0.5)
color_segments(n.diffexp.genes, 1:28, T)
polygon(x=c(1,1:28, 28), y=c(0, n.diffexp.genes$stress.duration+2000, 0), col="#00000011", border=F)
#polygon(x=c(bp.x.coords[11] - 0.6, bp.x.coords[11] - 0.6, bp.x.coords[12] + 0.6, bp.x.coords[12] + 0.6), y=c(-7000,max(n.diffexp.genes$up),max(n.diffexp.genes$up),-7000), col="#0000FF22", border=F, xpd=T)
dev.off()
