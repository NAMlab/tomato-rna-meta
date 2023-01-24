library(ggplot2)

d = read.csv("results/combined_deepest_levels.csv")

labels = scan("clades_hierarchy.txt", character())[(min(d$deepest_level)+1):(max(d$deepest_level)+1)]
labels[length(labels)] = "S. lycopersicum"

pdf("results/plots/4.1_gene_age.pdf")
p <- ggplot(d, aes(x=deepest_level, fill=set)) +
  geom_histogram(color="#e9ecef", alpha=0.6, position = 'dodge', binwidth = 1) +
  scale_fill_manual(values=c("#69b3a2", "#404080")) +
  scale_x_continuous(labels=labels, breaks = 0:16) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
print(p)
dev.off()
