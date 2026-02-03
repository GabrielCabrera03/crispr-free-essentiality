library(data.table)
IN_FILE <- "data/processed/gene_features_v1.tsv.gz"
OUT_FIG <- "results/figures/coefficients_full_glm.csv"

features <- c(
  "mean_expr",
  "var_expr",
  "loeuf",
  "loeuf_missing",
  "log_paralog_count",
  "log_ppi_degree"
)

fml <- as.formula(
  paste("Essential_label ~", paste(features, collapse = " + "))
)

train_std <- copy(train)
for (f in features) {
  train_std[[f]] <- (train_std[[f]] - mean(train_std[[f]], na.rm = TRUE)) / sd(train_std[[f]], na.rm = TRUE)
}

fit_std <- glm(
  fml, 
  data = train_std,
  family = binomial(),
  weights = train_std$w
)

coef_mat <- summary(fit_std)$coefficients

coef_dt <- data.table(
  feature = row.names(coef_mat),
  beta = coef_mat[, "Estimate"],
  se = coef_mat[, "Std. Error"],
  z = coef_mat[, "z value"],
  p = coef_mat[, "Pr(>|z|)"]
)

#remove intercept
coef_dt <- coef_dt[feature != "(Intercept"]

# add 95% CI

coef_dt[, `:=`(
  ci_low = beta - 1.96 * se,
  ci_high = beta +1.96 * se
)]

#sorting by effect size

setorder(coef_dt, beta)

print(coef_dt)

coef_dt[, odds_ratio := exp(beta)]
coef_dt[, or_ci_low := exp(ci_low)]
coef_dt[, or_ci_high := exp(ci_high)]

fwrite(coef_dt, OUT_FIG)

dotchart(coef_dt$beta, labels = coef_dt$feature,
         main = "Standardized coefficients (log-odds per 1 SD)",
         xlab = "beta")
abline(v = 0, lty = 2)

png("results/figures/coef_dotplot.png", width=900, height=600)
dotchart(coef_dt$beta, labels = coef_dt$feature,
         main = "Standardized coefficients (log-odds per 1 SD)",
         xlab = "beta")
abline(v = 0, lty = 2)
dev.off()


