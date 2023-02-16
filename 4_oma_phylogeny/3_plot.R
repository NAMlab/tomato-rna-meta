d = read.csv("output/combined_deepest_levels.csv")

freqs = reshape(data.frame(table(d$deepest_level, d$set)), direction="wide", timevar="Var2", idvar="Var1")
freqs$Var1 = as.numeric(freqs$Var1)
freqs$p.heat = cumsum(freqs$Freq.heat)/sum(freqs$Freq.heat) * 100
freqs$p.drought = cumsum(freqs$Freq.drought)/sum(freqs$Freq.drought) * 100
freqs$p.salt = cumsum(freqs$Freq.salt)/sum(freqs$Freq.salt) * 100
freqs$p.random = cumsum(freqs$Freq.random)/sum(freqs$Freq.random) * 100

labels = scan("clades_hierarchy.txt", character())[(min(d$deepest_level)+1):(max(d$deepest_level)+1)]
labels[length(labels)] = "S. lycopersicum"

pdf("output/plots/4.1_oldest_orthologs.pdf")
par(mar=c(8,4,2,0), family="serif")
plot(-freqs$Var1, freqs$p.random, ylim=c(0,100), type="b", axes=F, xlab="", ylab="% of Genes with an ortholog present at...")
axis(side = 2, las=2)
axis(side = 1, labels=labels, at=-c(1:length(labels)), las=2)

lines(-freqs$Var1, freqs$p.salt, type="b", col="green")
lines(-freqs$Var1, freqs$p.drought, type="b", col="blue")
lines(-freqs$Var1, freqs$p.heat, type="b", col="red")
dev.off()
