#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(pROC)
})

# -----------------------------
# Paths
# -----------------------------
IN_FILE <- "data/processed/gene_features_v1.tsv.gz"
OUT_METRICS_ALL <- "results/tables/metrics_v1_multiseed.csv"
OUT_METRICS_SUMMARY <- "results/tables/metrics_v1_multiseed_summary.csv"

dir.create("results/tables", showWarnings = FALSE, recursive = TRUE)

# -----------------------------
# Settings
# -----------------------------
seeds <- c(1, 2, 3, 4, 5)
train_frac <- 0.8

# Feature sets (ablation)
feature_sets <- list(
  expr_only = c("mean_expr", "var_expr"),
  expr_loeuf = c("mean_expr", "var_expr", "loeuf", "loeuf_missing"),
  expr_loeuf_paralog = c("mean_expr", "var_expr", "loeuf", "loeuf_missing", "log_paralog_count"),
  full_v1 = c("mean_expr", "var_expr", "loeuf", "loeuf_missing",
              "log_paralog_count", "log_ppi_degree")
)

# -----------------------------
# Helper: median impute numeric vector
# -----------------------------
impute_median <- function(x) {
  x <- as.numeric(x)
  x[!is.finite(x)] <- NA_real_
  if (all(is.na(x))) return(x)
  x[is.na(x)] <- median(x, na.rm = TRUE)
  x
}

# -----------------------------
# Helper: robust AUROC + AUPRC
# -----------------------------
eval_metrics <- function(y_true, y_prob) {
  ok <- !is.na(y_true) & is.finite(y_prob)
  y_true <- as.integer(y_true[ok])
  y_prob <- as.numeric(y_prob[ok])
  
  # Need both classes present
  if (length(unique(y_true)) < 2) {
    return(list(auroc = NA_real_, auprc = NA_real_))
  }
  
  # AUROC
  roc_obj <- pROC::roc(y_true, y_prob, quiet = TRUE)
  auroc <- as.numeric(pROC::auc(roc_obj))
  
  # Manual AUPRC
  ord <- order(y_prob, decreasing = TRUE)
  y <- y_true[ord]
  
  tp <- cumsum(y == 1L)
  fp <- cumsum(y == 0L)
  
  precision <- tp / pmax(tp + fp, 1)
  recall <- tp / sum(y == 1L)
  
  # Add origin point
  precision <- c(1, precision)
  recall <- c(0, recall)
  
  # Trapezoid integration over recall
  auprc <- sum((recall[-1] - recall[-length(recall)]) *
                 (precision[-1] + precision[-length(precision)]) / 2)
  
  list(auroc = auroc, auprc = auprc)
}

# -----------------------------
# Load data
# -----------------------------
dt <- fread(IN_FILE)

# Keep labeled rows only
dt <- dt[!is.na(Essential_label)]
dt[, Essential_label := as.integer(Essential_label)]
dt <- dt[Essential_label %in% c(0L, 1L)]

# LOEUF missingness + impute
dt[, loeuf_missing := as.integer(is.na(loeuf))]
dt[, loeuf := impute_median(loeuf)]

# Defensive impute for other numeric features used
if ("mean_expr" %in% names(dt)) dt[, mean_expr := impute_median(mean_expr)]
if ("var_expr" %in% names(dt)) dt[, var_expr := impute_median(var_expr)]
if ("log_paralog_count" %in% names(dt)) dt[, log_paralog_count := impute_median(log_paralog_count)]
if ("log_ppi_degree" %in% names(dt)) dt[, log_ppi_degree := impute_median(log_ppi_degree)]

# -----------------------------
# Run multi-seed evaluation
# -----------------------------
metrics_all <- data.table(
  seed = integer(),
  model = character(),
  auroc = numeric(),
  auprc = numeric(),
  prevalence = numeric(),
  n_test = integer(),
  n_pos_test = integer()
)

for (s in seeds) {
  set.seed(s)
  
  n <- nrow(dt)
  train_idx <- sample.int(n, size = floor(train_frac * n))
  train <- dt[train_idx]
  test  <- dt[-train_idx]
  
  prevalence <- mean(test$Essential_label)
  n_test <- nrow(test)
  n_pos_test <- sum(test$Essential_label == 1L)
  
  # Class weights
  w_pos <- sum(train$Essential_label == 0L) / sum(train$Essential_label == 1L)
  train[, w := ifelse(Essential_label == 1L, w_pos, 1)]
  # Ensure no NA weights
  train[!is.finite(w) | is.na(w), w := 1]
  
  for (nm in names(feature_sets)) {
    feats <- feature_sets[[nm]]
    
    # Ensure all features exist
    missing_feats <- setdiff(feats, names(train))
    if (length(missing_feats) > 0) {
      metrics_all <- rbind(metrics_all, data.table(
        seed = s, model = nm, auroc = NA_real_, auprc = NA_real_,
        prevalence = prevalence, n_test = n_test, n_pos_test = n_pos_test
      ))
      next
    }
    
    fml <- as.formula(paste("Essential_label ~", paste(feats, collapse = " + ")))
    
    fit <- tryCatch(
      glm(
        fml,
        data = train,
        family = binomial(),
        weights = train$w,
        control = glm.control(maxit = 50)
      ),
      error = function(e) NULL
    )
    
    if (is.null(fit)) {
      metrics_all <- rbind(metrics_all, data.table(
        seed = s, model = nm, auroc = NA_real_, auprc = NA_real_,
        prevalence = prevalence, n_test = n_test, n_pos_test = n_pos_test
      ))
      next
    }
    
    pred <- tryCatch(
      predict(fit, newdata = test, type = "response"),
      error = function(e) rep(NA_real_, nrow(test))
    )
    
    m <- eval_metrics(test$Essential_label, pred)
    
    metrics_all <- rbind(metrics_all, data.table(
      seed = s, model = nm,
      auroc = m$auroc,
      auprc = m$auprc,
      prevalence = prevalence,
      n_test = n_test,
      n_pos_test = n_pos_test
    ))
  }
  
  cat("Finished seed:", s, "\n")
}

# -----------------------------
# Summarize
# -----------------------------
metrics_summary <- metrics_all[, .(
  mean_auroc = mean(auroc, na.rm = TRUE),
  sd_auroc   = sd(auroc, na.rm = TRUE),
  mean_auprc = mean(auprc, na.rm = TRUE),
  sd_auprc   = sd(auprc, na.rm = TRUE),
  mean_prev  = mean(prevalence, na.rm = TRUE),
  n_runs     = sum(!is.na(auprc))
), by = model]

setorder(metrics_summary, -mean_auprc)

print(metrics_summary)

# -----------------------------
# Save
# -----------------------------
fwrite(metrics_all, OUT_METRICS_ALL)
fwrite(metrics_summary, OUT_METRICS_SUMMARY)

cat("\nWrote per-seed metrics to:", OUT_METRICS_ALL, "\n")
cat("Wrote summary metrics to:", OUT_METRICS_SUMMARY, "\n")
