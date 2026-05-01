library(edgeR)
library(limma)
library(stringr)

# Set TRUE for fast testing: all 338 samples but only 3000 random genes (gives a
# meaningful corfit value in seconds). Set FALSE for the full ~30k transcripts.
SUBSET <- FALSE

abundance <- read.csv("input/combined_abundance.tsv.gz", sep = "\t", check.names = FALSE)
row.names(abundance) <- abundance$target_id

samples.annotation <- read.csv("input/samples_annotation.csv")
samples.annotation[samples.annotation$sra_run_id == "", ]$sra_run_id <-
  samples.annotation[samples.annotation$sra_run_id == "", ]$sample.name
samples.annotation$sample.group <- make.names(samples.annotation$sample.group)

samples <- samples.annotation
samples$sample_id <- paste0(samples$sra_run_id, "_est_counts")

samples$genotype.name[samples$genotype.name == ""] <- "(unspecified)"

samples$datasource_id <- factor(samples$datasource_id)
samples$stress_type   <- factor(samples$stress_type, levels = c("none", "heat", "drought", "salt"))
# Pick reference levels with the broadest cross-dataset coverage so the named
# coefficients carry meaningful interpretations.
samples$genotype.name <- relevel(factor(samples$genotype.name), ref = "M82")
samples$tissue        <- relevel(factor(samples$tissue),        ref = "leaf")

samples$stress.duration <- suppressWarnings(as.numeric(samples$stress.duration))
samples$temperature     <- suppressWarnings(as.numeric(samples$temperature))
samples$heat_duration_known <- samples$stress_type == "heat" & !is.na(samples$stress.duration)

counts <- as.matrix(abundance[, samples$sample_id, drop = FALSE])
samples <- samples[match(colnames(counts), samples$sample_id), ]
stopifnot(all(samples$sample_id == colnames(counts)))

design <- model.matrix(~ stress_type + genotype.name + tissue, data = samples)
colnames(design) <- make.names(colnames(design))
stopifnot(nrow(design) == ncol(counts))
cat("design:", nrow(design), "x", ncol(design),
    " rank:", qr(design)$rank, "\n")

dge <- DGEList(counts = counts)

keep <- filterByExpr(dge, design = design, min.count = 10, min.total.count = 15)
dge <- dge[keep, , keep.lib.sizes = FALSE]
dge <- calcNormFactors(dge)

if (SUBSET) {
  set.seed(1)
  dge <- dge[sample(seq_len(nrow(dge)), min(3000, nrow(dge))), , keep.lib.sizes = FALSE]
  cat("[SUBSET mode] testing with", nrow(dge), "random genes\n")
}

v <- voom(dge, design, plot = TRUE)
corfit <- duplicateCorrelation(v, design, block = samples$datasource_id)
cat("Intra-datasource correlation:", round(corfit$consensus, 3), "\n")

# Re-run voom now that we have an estimate of the within-block correlation.
v <- voom(dge, design, block = samples$datasource_id,
          correlation = corfit$consensus, plot = TRUE)
corfit <- duplicateCorrelation(v, design, block = samples$datasource_id)
cat("Refined intra-datasource correlation:", round(corfit$consensus, 3), "\n")

fit <- lmFit(v, design,
             block = samples$datasource_id,
             correlation = corfit$consensus)
fit <- eBayes(fit)

cat("\nAvailable coefficients:\n"); print(colnames(design))

# =============================================================================
# Contrasts and per-contrast top tables
# =============================================================================
out_dir <- "output/voom"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# All stress contrasts: each stress vs control + each pair against each other.
# Tissue/genotype contrasts: only the levels with multi-datasource coverage,
# so the cross-dataset average actually identifies the effect (see caveats doc).
contrast_matrix <- makeContrasts(
  heat_vs_none      = stress_typeheat,
  drought_vs_none   = stress_typedrought,
  salt_vs_none      = stress_typesalt,
  heat_vs_drought   = stress_typeheat - stress_typedrought,
  heat_vs_salt      = stress_typeheat - stress_typesalt,
  drought_vs_salt   = stress_typedrought - stress_typesalt,
  anther_vs_leaf    = tissueanther,
  fruit_vs_leaf     = tissuefruit,
  seedling_vs_leaf  = tissueseedling,
  microtom_vs_m82   = genotype.nameMicro.tom,
  moneymaker_vs_m82 = genotype.nameMoneymaker,
  levels = design
)

fit2 <- contrasts.fit(fit, contrast_matrix)
fit2 <- eBayes(fit2)

write_contrast_table <- function(fit_obj, cname, dir = out_dir) {
  tt <- topTable(fit_obj, coef = cname, number = Inf, sort.by = "none")
  path <- file.path(dir, paste0("voom_", cname, ".csv"))
  write.csv(tt, path)
  n_sig <- sum(tt$adj.P.Val < 0.05, na.rm = TRUE)
  cat(sprintf("  %-22s -> %d genes, %d at FDR<0.05\n", cname, nrow(tt), n_sig), file = stdout())
  invisible(tt)
}

cat("\nWriting contrast tables to", out_dir, ":\n")
for (cname in colnames(contrast_matrix)) write_contrast_table(fit2, cname)

saveRDS(list(fit = fit, fit2 = fit2, design = design, corfit = corfit, samples = samples),
        file.path(out_dir, "voom_main_fit.rds"))


# =============================================================================
# Secondary analysis: heat duration time-course
# =============================================================================
# Restricted to heat samples with a recorded duration (n = 56 across 11
# datasources). Avoids fabricating durations for samples without recorded
# duration. The coefficient `log_heat_duration` is the log2 fold-change per
# unit of log(hours+1) — the within-heat slope on a log time axis (raw
# durations span 1 h to 1200 h, so log compresses the range).
#
# Linear-on-log-time captures monotonic responses. For non-monotonic genes
# (HSPs that peak then decline), use the spline form below; the F-test from
# topTable across the spline coefficients tells you whether the time course
# is non-flat at all, without assuming a specific shape.
# =============================================================================
cat("\n=== Heat duration sub-model ===\n")

samples_dur <- samples[samples$heat_duration_known, ]
counts_dur  <- counts[, samples$heat_duration_known]
samples_dur$log_heat_duration <- log(samples_dur$stress.duration + 1)
cat(nrow(samples_dur), "samples across",
    length(unique(samples_dur$datasource_id)), "datasources\n")

design_dur <- model.matrix(~ log_heat_duration, data = samples_dur)
# Spline alternative (df = 3 captures rise-and-fall shapes):
# library(splines)
# design_dur <- model.matrix(~ ns(log_heat_duration, df = 3), data = samples_dur)

dge_dur <- DGEList(counts = counts_dur)
keep_dur <- filterByExpr(dge_dur, design = design_dur,
                         min.count = 10, min.total.count = 15)
dge_dur <- dge_dur[keep_dur, , keep.lib.sizes = FALSE]
dge_dur <- calcNormFactors(dge_dur)
if (SUBSET) {
  set.seed(1)
  dge_dur <- dge_dur[sample(seq_len(nrow(dge_dur)), min(3000, nrow(dge_dur))), ,
                     keep.lib.sizes = FALSE]
}

v_dur <- voom(dge_dur, design_dur, plot = FALSE)
corfit_dur <- duplicateCorrelation(v_dur, design_dur,
                                   block = samples_dur$datasource_id)
cat("Intra-datasource correlation (duration model):",
    round(corfit_dur$consensus, 3), "\n")
v_dur <- voom(dge_dur, design_dur, block = samples_dur$datasource_id,
              correlation = corfit_dur$consensus, plot = FALSE)
fit_dur <- lmFit(v_dur, design_dur,
                 block = samples_dur$datasource_id,
                 correlation = corfit_dur$consensus)
fit_dur <- eBayes(fit_dur)

write_contrast_table(fit_dur, "log_heat_duration")
saveRDS(list(fit = fit_dur, design = design_dur, corfit = corfit_dur,
             samples = samples_dur),
        file.path(out_dir, "voom_heat_duration_fit.rds"))
