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

mrca = read.csv("output/combined_deepest_levels.csv")
mrca = mrca[mrca$set == "heat",]
mrca = merge(mrca, internal.names)[c("level", "internal.name")]
row.names(mrca) = mrca$internal.name
mrca.vec = mrca$level
names(mrca.vec) = mrca$internal.name
mrca.diffs = dist(mrca.vec)

# Reshape mrca.diffs into long format
library(reshape2)
mrca.diffs.long = melt(as.matrix(mrca.diffs))

# Reshape evo.all into long format
evo.all.long = melt(as.matrix(evo.all), varnames = c("Var1", "Var2"), value.name="jaccard")
merged = merge(mrca.diffs.long, evo.all.long)

plot(merged$jaccard ~ merged$value)


lm_model <- lm(jaccard ~ value, data = merged)
exp_model <- nls(jaccard ~ a * exp(-b * value), data = merged, start = list(a = 1, b = 1))

r_squared_lm <- summary(lm_model)$r.squared

# Calculate AIC and BIC for the linear model
n <- nrow(merged)
k_lm <- length(coefficients(lm_model))
aic_lm <- AIC(lm_model)
bic_lm <- BIC(lm_model)

# Calculate AIC and BIC for the exponential model
k_exp <- length(coefficients(exp_model))
aic_exp <- AIC(exp_model)
bic_exp <- BIC(exp_model)

# Calculate mean squared error for the exponential model
mse_exp <- mean(residuals(exp_model)^2)
rmse_exp <- sqrt(mse_exp)

plot(merged$jaccard ~ merged$value)
abline(lm_model, col = "red")
lines(0:14, predict(exp_model, newdata = data.frame(value = 0:14)), col = "blue")

print(paste("R-squared (Linear Model):", r_squared_lm))
print(paste("AIC (Linear Model):", aic_lm))
print(paste("BIC (Linear Model):", bic_lm))

print(paste("AIC (Exponential Model):", aic_exp))
print(paste("BIC (Exponential Model):", bic_exp))
print(paste("Mean Squared Error (Exponential Model):", mse_exp))
