suppressPackageStartupMessages({
  library(data.table)
  library(pROC)
})

# -----------------------------
# Paths
# -----------------------------
IN_FILE <- "data/processed/gene_features_v1.tsv.gz"
OUT_DIR <- "results/tables"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
OUT_PRED <- file.path(OUT_DIR, "logreg_test_predictions.csv")
OUT_SUMMARY <- file.path(OUT_DIR, "logreg_training_summary.csv")
OUT_COEF <- file.path(OUT_DIR, "logreg_coefficients.csv")   # optional but recommended

# -----------------------------
# Settings
# -----------------------------
SEED <- 42
TRAIN_FRAC <- 0.8

# Full feature set (same as your best-performing model)
FEATURES <- c(
  "mean_expr",
  "var_expr",
  "loeuf",
  "log_paralog_count",
  "log_ppi_degree"
)

# -----------------------------
# Helpers
# -----------------------------
impute_median <- function(x) {
  x <- as.numeric(x)
  x[!is.finite(x)] <- NA_real_
  if (all(is.na(x))) return(x)
  x[is.na(x)] <- median(x, na.rm = TRUE)
  x
}

# Manual AUPRC (robust; avoids PRROC crashes)
auprc_manual <- function(y_true, y_prob) {
  ok <- !is.na(y_true) & is.finite(y_prob)
  y_true <- as.integer(y_true[ok])
  y_prob <- as.numeric(y_prob[ok])
  
  if (length(unique(y_true)) < 2) return(NA_real_)
  
  ord <- order(y_prob, decreasing = TRUE)
  y <- y_true[ord]
  
  tp <- cumsum(y == 1L)
  fp <- cumsum(y == 0L)
  
  precision <- tp / pmax(tp + fp, 1)
  recall <- tp / sum(y == 1L)
  
  precision <- c(1, precision)
  recall <- c(0, recall)
  
  sum((recall[-1] - recall[-length(recall)]) *
        (precision[-1] + precision[-length(precision)]) / 2)
}

# AUROC helper
auroc <- function(y_true, y_prob) {
  ok <- !is.na(y_true) & is.finite(y_prob)
  y_true <- as.integer(y_true[ok])
  y_prob <- as.numeric(y_prob[ok])
  
  if (length(unique(y_true)) < 2) return(NA_real_)
  
  roc_obj <- pROC::roc(y_true, y_prob, quiet = TRUE)
  as.numeric(pROC::auc(roc_obj))
}

# -----------------------------
# Load + preprocess
# -----------------------------
dt <- fread(IN_FILE)

# Keep labeled rows only
dt <- dt[!is.na(Essential_label)]
dt[, Essential_label := as.integer(Essential_label)]
dt <- dt[Essential_label %in% c(0L, 1L)]

# Create LOEUF missingness flag + impute
dt[, loeuf_missing := as.integer(is.na(loeuf))]
dt[, loeuf := impute_median(loeuf)]

# Defensive imputations for other numeric predictors
dt[, mean_expr := impute_median(mean_expr)]
dt[, var_expr := impute_median(var_expr)]
dt[, log_paralog_count := impute_median(log_paralog_count)]
dt[, log_ppi_degree := impute_median(log_ppi_degree)]

# Quick sanity checks
stopifnot(all(FEATURES %in% names(dt)))
cat("Loaded rows:", nrow(dt), "\n")
cat("Prevalence overall:", mean(dt$Essential_label), "\n")



#Train

set.seed(SEED)

train_idx <- sample(nrow(dt), 0.8 * nrow(dt))
train_data <-dt[train_idx, ]
test_data <- dt[-train_idx, ]
for (f in FEATURES) {
  train[[f]] <- impute_median(train[[f]])
  test[[f]]  <- impute_median(test[[f]])
}


#fit logistic regression model
fml <- as.formula(paste("Essential_label ~", paste(FEATURES, collapse = " + ")))

# Optional: class weights to handle imbalance
w_pos <- sum(train$Essential_label == 0L) / sum(train$Essential_label == 1L)
train[, w := ifelse(Essential_label == 1L, w_pos, 1)]
train[!is.finite(w) | is.na(w), w := 1]

fit <- tryCatch(
  glm(
    fml,
    data = train,
    family = binomial(),
    weights = w,
    control = glm.control(maxit = 50)
  ),
  error = function(e) e
)

if (inherits(fit, "error")) {
  cat("\nGLM FAILED WITH ERROR:\n")
  cat(fit$message, "\n")
  quit(status = 1)
}

#prediction

predictions_prob <- predict(fit, newdata = test, type = "response")
cat("Pred NA count:", sum(is.na(predictions_prob)), "\n")


#Evaluation 

prev <- mean(test$Essential_label)
cat("Prevalence:", prev, "\n")
cat("AUROC:", auroc(test$Essential_label, predictions_prob), "\n")
cat("AUPRC:", auprc_manual(test$Essential_label, predictions_prob), "\n")
summary(fit)


# ---- Prediction table (test set) ----
# choose a simple threshold for now (we’ll tune later)
THRESH <- 0.5

pred_dt <- data.table(
  gene_symbol = test$gene_symbol,
  y_true = test$Essential_label,
  prob_essential = predictions_prob,
  pred_label = as.integer(predictions_prob >= THRESH)
)

# add rank (1 = most likely essential)
setorder(pred_dt, -prob_essential)
pred_dt[, rank := .I]

fwrite(pred_dt, OUT_PRED)
cat("Wrote test predictions to:", OUT_PRED, "\n")


# ---- Training summary ----
summary_dt <- data.table(
  model = "logreg_full_v1",
  seed = SEED,
  train_frac = TRAIN_FRAC,
  n_train = nrow(train),
  n_test = nrow(test),
  prev_test = prev,
  auroc = auroc(test$Essential_label, predictions_prob),
  auprc = auprc_manual(test$Essential_label, predictions_prob),
  pred_na = sum(is.na(predictions_prob))
)

fwrite(summary_dt, OUT_SUMMARY)
cat("Wrote training summary to:", OUT_SUMMARY, "\n")
print(summary_dt)

# ---- Coefficients table (interpretable) ----
coef_mat <- summary(fit)$coefficients
coef_dt <- data.table(
  term = rownames(coef_mat),
  estimate = coef_mat[, "Estimate"],
  std_error = coef_mat[, "Std. Error"],
  z = coef_mat[, "z value"],
  p = coef_mat[, "Pr(>|z|)"]
)

fwrite(coef_dt, OUT_COEF)
cat("Wrote coefficients to:", OUT_COEF, "\n")
