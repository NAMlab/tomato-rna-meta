# Load required libraries
library(data.table)
library(ggplot2)
library(R.utils)
library(gridExtra)

# Read circadian gene transcript IDs
circadian_genes <- fread("input/circadian_clock_genes.csv")
transcript_ids <- circadian_genes$transcript_id

# Read sample annotation and filter for DS-13 samples
samples_anno <- fread("input/samples_annotation.csv")
ds13_samples <- samples_anno[grepl("^DS-13", samples_anno$datasource_id), ]

# Extract relevant columns: sample name, stress.duration, stress_type
ds13_samples <- ds13_samples[, .(sample.name, stress.duration, stress_type)]

# Map sample names to treatment and ensure factor levels
ds13_samples[, treatment := ifelse(stress_type == "none", "control (25°C)", "heat (37°C)")]
ds13_samples[, treatment := factor(treatment, levels = c("heat (37°C)", "control (25°C)"))]

# Get sample names for abundance file columns
sample_names <- ds13_samples$sample.name

# Read abundance data (TPM values only for relevant samples and transcripts)
# Read header to get column indices
header <- fread("input/combined_abundance.tsv.gz", nrows = 0)
tpm_cols <- paste0(sample_names, "_tpm")
cols_to_read <- c("target_id", tpm_cols)

# Read only needed columns
abundance <- fread("input/combined_abundance.tsv.gz", select = cols_to_read)
abundance <- abundance[target_id %in% transcript_ids]

# Melt abundance data for plotting
abundance_long <- melt(
  abundance,
  id.vars = "target_id",
  variable.name = "sample_tpm",
  value.name = "TPM"
)

# Extract sample name from column name
abundance_long[, sample.name := sub("_tpm$", "", sample_tpm)]

# Merge with annotation to get time and treatment
plot_data <- merge(
  abundance_long,
  ds13_samples,
  by = "sample.name"
)

# Merge with gene names for legend (use 'name' instead of 'gene_name')
plot_data <- merge(
  plot_data,
  circadian_genes[, .(name, transcript_id)],
  by.x = "target_id",
  by.y = "transcript_id"
)

# Log-transform TPM values
plot_data[, logTPM := log2(TPM + 1)]

# Plot: facet by gene name, show heat and control together for each gene
p <- ggplot(plot_data, aes(x = stress.duration, y = logTPM, color = treatment, group = treatment)) +
  annotate("rect", xmin = 13, xmax = 24, ymin = -Inf, ymax = Inf, fill = "grey90", color = NA) +
  geom_line() +
  geom_point() +
  facet_wrap(~ name, ncol = 1, scales = "free_y") +
  scale_color_manual(
    values = c("heat (37°C)" = "red", "control (25°C)" = "blue"),
    breaks = c("heat (37°C)", "control (25°C)"),
    labels = c("Heat (37°C)", "Control (25°C)")
  ) +
  labs(
    title = "Log2(TPM+1) of Circadian Clock Genes Over Time (DS-13)",
    x = "Stress Duration (hours)",
    y = "log2(TPM + 1)",
    color = "Treatment"
  ) +
  theme_bw() +
  theme(
    legend.position = "bottom",
    strip.background = element_rect(fill = "grey90")
  )

# Save to PDF
dir.create("output/plots", showWarnings = FALSE, recursive = TRUE)
pdf("output/plots/1.10_circadian_clock_genes.pdf", width = 8, height = 14)
grid::grid.newpage()
grid::grid.draw(gridExtra::arrangeGrob(
  p,
  bottom = grid::textGrob(
    "Gray areas indicate dark period (13-24h).\nGenes identified by BLASTing Arabidopsis genes, see https://doi.org/10.1111/tpj.70383",
    x = 0.99, hjust = 1, gp = grid::gpar(fontsize = 10, col = "black", fontface = "italic")
  )
))
dev.off()
