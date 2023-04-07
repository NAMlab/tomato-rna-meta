d = read.csv("output/combined_deepest_levels.csv")

freqs = reshape(data.frame(table(d$level, d$set)), direction="wide", timevar="Var2", idvar="Var1")
freqs$Var1 = as.numeric(freqs$Var1)
freqs$p.heat = cumsum(freqs$Freq.heat)/sum(freqs$Freq.heat) * 100
freqs$p.drought = cumsum(freqs$Freq.drought)/sum(freqs$Freq.drought) * 100
freqs$p.salt = cumsum(freqs$Freq.salt)/sum(freqs$Freq.salt) * 100
freqs$p.random = cumsum(freqs$Freq.random)/sum(freqs$Freq.random) * 100

labels = scan("input/clades_hierarchy.txt", character())[(min(d$level)+1):(max(d$level))]
# Remove the last class (S. lycopersicum -- no orthologs within own species)
freqs = freqs[freqs$Var1 < 17,]

# Prepare for plotting the genes at each level
d.h = d[d$set == "heat" & d$level < 16,]
d.h$level = d.h$level + 1
d.h = merge(d.h, read.csv("../1_core_response/input/hs_core_genes_internal_names.csv"))

plot.cols = list(heat="#D55E00", drought="#0072B2", salt="#009E73", baseline="black")

pdf("output/plots/4.1_oldest_orthologs.pdf", 6.5, 3)
par(mar=c(2,7,1,4), family="serif", cex=0.6)
plot(freqs$p.random, -freqs$Var1, type="b", xlim=c(0,100), axes=F, xlab="", ylab="", xaxs="i")
axis(side = 2, labels=labels, at=-c(1:length(labels)), las=2, col=F, col.axis="gray52")
axis(side = 1, col=F, col.axis="gray52", at=seq(0,100,20), labels=paste0(seq(0,100,20), "%"))

text(10, -14, "proportion of core genes \nwith ortholog at ancestor node...", pos=4, cex=1.2, xpd=T)

lines(freqs$p.salt,-freqs$Var1, type="b", col=plot.cols[["salt"]])
lines(freqs$p.drought,-freqs$Var1, type="b", col=plot.cols[["drought"]])
lines(freqs$p.heat,-freqs$Var1, type="b", col=plot.cols[["heat"]])

text(35, -2.7, "salt", col=plot.cols[["salt"]])
text(34.5, -5, "drought", col=plot.cols[["drought"]], pos=2)
text(52, -2.3, "heat", col=plot.cols[["heat"]], pos=2)
text(21, -2.7, "random baseline", col="black", pos=2)
text(86, -5, "heat stress core genes\noriginating at this level", col="gray70", pos=4, xpd=T)

# Plot 
for(l in 1:16) {
  x = freqs[freqs$Var1 == l,]$p.heat
  if(l == 2) {
    # Special case (huge block of genes)
    text(x, -3, 
         paste(strwrap(paste(d.h[d.h$level == l,]$internal.name, collapse=", "), 40), collapse="\n"), 
         pos=4, col="gray70", xpd=T, cex=0.7)
  } else {
    text(x, -l, 
         paste(strwrap(paste(d.h[d.h$level == l,]$internal.name, collapse=", "), 115-x), collapse="\n"), 
         pos=4, col="gray70", xpd=T, cex=0.7)
  }
}

dev.off()
