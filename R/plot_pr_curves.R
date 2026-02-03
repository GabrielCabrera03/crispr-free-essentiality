
suppressPackageStartupMessages({
  library(data.table)
  library(PRROC)
  library(ggplot2)
})

# -----------------------------
# Paths
# -----------------------------
IN_FILE <- "data/processed/gene_features_v1.tsv.gz"
OUT_FIG <- "results/figures/pr_curves.png"

dir.create("results/figures", showWarnings = FALSE, recursive = TRUE)

# -----------------------------
# Load data
# -----------------------------

dt <- fread(IN_FILE)

#keep labeled rows
dt <- dt[!is.na(Essential_label)]
dt[, Essential_label := as.integer(Essential_label)]

# Impute minimal (must match modeling script)
dt[, loeuf_missing := as.integer(is.na(loeuf))]
dt[is.na(loeuf), loeuf := median(loeuf, na.rm = TRUE)]

impute_median <- function(x) {
  x[!is.finite(x)] <- NA_real_
  x[is.na(x)] <- median(x, na.rm = TRUE)
  x
}

dt[, mean_expr := impute_median(mean_expr)]
dt[, var_expr := impute_median(var_expr)]
dt[, log_paralog_count := impute_median(log_paralog_count)]
dt[, log_ppi_degree := impute_median(log_ppi_degree)]

# -----------------------------
# Train/test split (same seed!)
# -----------------------------
set.seed(42)
n <- nrow(dt)
train_idx <- sample(seq_len(n), size = floor(0.8 * n))

train <- dt[train_idx]
test  <- dt[-train_idx]

# Class weights
w_pos <- sum(train$Essential_label == 0) / sum(train$Essential_label == 1)
train[, w := ifelse(Essential_label == 1, w_pos, 1)]

# -----------------------------
# Feature sets
# -----------------------------
feature_sets <- list(
  "Expression only" = c("mean_expr", "var_expr"),
  "+ LOEUF" = c("mean_expr", "var_expr", "loeuf", "loeuf_missing"),
  "+ Paralog" = c("mean_expr", "var_expr", "loeuf", "loeuf_missing", "log_paralog_count"),
  "Full model" = c("mean_expr", "var_expr", "loeuf", "loeuf_missing",
                   "log_paralog_count", "log_ppi_degree")
)

# -----------------------------
# Compute PR curves
# -----------------------------
pr_df <- list()

for (nm in names(feature_sets)) {
  feats <- feature_sets[[nm]]
  fml <- as.formula(paste("Essential_label ~", paste(feats, collapse = " + ")))
  
  fit <- glm(
    fml,
    data = train,
    family = binomial(),
    weights = train$w
  )
  
  prob <- predict(fit, newdata = test, type = "response")
  
  # PRROC expects: class0 = positives
  pos <- prob[test$Essential_label == 1]
  neg <- prob[test$Essential_label == 0]
  
  pos <- pos[is.finite(pos)]
  neg <- neg[is.finite(neg)]
  
  pr <- PRROC::pr.curve(
    scores.class0 = pos,
    scores.class1 = neg,
    curve = TRUE
  )
  
  pr_df[[nm]] <- data.table(
    recall = pr$curve[,1],
    precision = pr$curve[,2],
    model = nm
  )
}

pr_df <- rbindlist(pr_df)

# Baseline prevalence
prevalence <- mean(test$Essential_label)

# -----------------------------
# Plot
# -----------------------------
p <- ggplot(pr_df, aes(x = recall, y = precision, color = model)) +
  geom_line(linewidth = 1.1) +
  geom_hline(
    yintercept = prevalence,
    linetype = "dashed",
    color = "black"
  ) +
  labs(
    title = "Precision–Recall curves for CRISPR-free essentiality prediction",
    subtitle = paste0("Dashed line = prevalence baseline (", round(prevalence, 3), ")"),
    x = "Recall",
    y = "Precision",
    color = "Model"
  ) +
  theme_bw(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold"),
    legend.position = "right"
  )

ggsave(OUT_FIG, p, width = 8, height = 6, dpi = 300)

cat("Saved PR curve plot to:", OUT_FIG, "\n")