d = read.csv("results/combined_deepest_levels.csv")

freqs = reshape(data.frame(table(d$deepest_level, d$set)), direction="wide", timevar="Var2", idvar="Var1")
freqs$Var1 = as.numeric(freqs$Var1)
freqs$p.hs_core = cumsum(freqs$Freq.hs_core)/sum(freqs$Freq.hs_core) * 100
freqs$p.random = cumsum(freqs$Freq.random)/sum(freqs$Freq.random) * 100

labels = scan("clades_hierarchy.txt", character())[(min(d$deepest_level)+1):(max(d$deepest_level)+1)]
labels[length(labels)] = "S. lycopersicum"

pdf("results/plots/4.1_oldest_orthologs.pdf")
par(mar=c(8,4,2,0), family="serif")
plot(freqs$Var1, freqs$p.random, ylim=c(0,100), type="b", axes=F, xlab="", ylab="% of Genes with an Ancestor within...")
axis(side = 2, las=2)
axis(side = 1, labels=labels, at=c(1:length(labels)), las=2)

lines(freqs$Var1, freqs$p.hs_core, type="b", col="red")
dev.off()