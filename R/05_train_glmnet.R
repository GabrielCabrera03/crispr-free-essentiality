
suppressPackageStartupMessages({
  library(data.table)
  library(glmnet)
  library(pROC)
})

#Paths
args    <- commandArgs(trailingOnly = TRUE)
VERSION <- if (length(args) >= 1) args[1] else "v1"

IN_FILE     <- "data/processed/gene_features_v1.tsv.gz"
OUT_DIR     <- "results/tables"
OUT_METRICS <- file.path(OUT_DIR, paste0("glmnet_", VERSION, "_metrics.csv"))
OUT_PRED    <- file.path(OUT_DIR, paste0("glmnet_", VERSION, "_test_predictions.csv"))
OUT_COEF    <- file.path(OUT_DIR, paste0("glmnet_", VERSION, "_coefficients.csv"))
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

#Settings

SEED <- 42
TRAIN_FRAC <- 0.8

# Elastic net mixing:
# alpha=1 -> lasso, alpha=0 -> ridge, alpha=0.5 -> mix
ALPHA <- 0.5
N_FOLDS <- 5

FEATURES <- c(
  "mean_expr",
  "var_expr",
  "loeuf",
  "log_paralog_count",
  "log_ppi_degree"
)

#Helpers

impute_median <- function(x) {
  x <- as.numeric(x)
  x[!is.finite(x)] <- NA_real_
  if (all(is.na(x))) return(x)
  x[is.na(x)] <- median(x, na.rm = TRUE)
  x
}

# Standardize using TRAIN stats only (no leakage)
standardize_train_test <- function(train_vec, test_vec) {
  mu <- mean(train_vec, na.rm = TRUE)
  sdv <- sd(train_vec, na.rm = TRUE)
  if (!is.finite(sdv) || sdv == 0) sdv <- 1
  list(
    train = (train_vec - mu) / sdv,
    test  = (test_vec  - mu) / sdv
  )
}

auroc <- function(y_true, y_prob) {
  ok <- !is.na(y_true) & is.finite(y_prob)
  y_true <- as.integer(y_true[ok])
  y_prob <- as.numeric(y_prob[ok])
  if (length(unique(y_true)) < 2) return(NA_real_)
  as.numeric(pROC::auc(pROC::roc(y_true, y_prob, quiet = TRUE)))
}

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

#Load + basic preprocess

dt <- fread(IN_FILE)
dt <- dt[!is.na(Essential_label)]
dt[, Essential_label := as.integer(Essential_label)]
dt <- dt[Essential_label %in% c(0L, 1L)]

stopifnot(all(FEATURES %in% names(dt)))
cat("Loaded rows:", nrow(dt), "\n")

#Train 

set.seed(SEED)

train_idx <- sample(nrow(dt), TRAIN_FRAC * nrow(dt))
train_data <- dt[train_idx, ]
test_data <- dt[-train_idx, ]

#impute using train-only medians, apply same transform to test
for (f in FEATURES) {
  train_data[[f]] <- impute_median(train_data[[f]])
  test_data[[f]]  <- impute_median(test_data[[f]])
}

#standardize features using TRAIN stats only

for (f in FEATURES) {
  s <- standardize_train_test(train_data[[f]], test_data[[f]])
  train_data[[f]] <- s$train
  test_data[[f]]  <- s$test
}

#building x_train, y_train, x_test, y_test

y_train <- train_data$Essential_label

x_train <- as.matrix(train_data[, ..FEATURES])

y_test <- test_data$Essential_label

x_test <- as.matrix(test_data[, ..FEATURES])

#computing weights for training

w_pos <- sum(y_train == 0L) / sum(y_train == 1L)
w <- ifelse(y_train == 1L, w_pos, 1)

#fit using glmnet

cvfit <- cv.glmnet(
  x = x_train, 
  y = y_train,
  family = "binomial",
  alpha = ALPHA, 
  weights = w,
  nfolds = N_FOLDS,
  type.measure = "auc"
)

cat("lambda.min:", cvfit$lambda.min, "\n")
cat("lambda.1se:", cvfit$lambda.1se, "\n")

#prediction

prob_min <- as.numeric(
  predict(cvfit, newx = x_test, s = "lambda.min", type = "response")
)
prob_1se <- as.numeric(
  predict(cvfit, newx = x_test, s = "lambda.1se", type = "response")
)

cat("NA prob_min:", sum(is.na(prob_min)), "\n")
cat("NA prob_1se:", sum(is.na(prob_1se)), "\n")

#evaluation
prev <- mean(y_test)

cat("Prevalence:", prev, "\n\n")

cat("lambda.min\n")
cat("  AUROC:", auroc(y_test, prob_min), "\n")
cat("  AUPRC:", auprc_manual(y_test, prob_min), "\n\n")

cat("lambda.1se\n")
cat("  AUROC:", auroc(y_test, prob_1se), "\n")
cat("  AUPRC:", auprc_manual(y_test, prob_1se), "\n")

# --- Metrics ---
metrics_dt <- data.table(
  model      = c("glmnet_lambda.min", "glmnet_lambda.1se"),
  auroc      = c(auroc(y_test, prob_min), auroc(y_test, prob_1se)),
  auprc      = c(auprc_manual(y_test, prob_min), auprc_manual(y_test, prob_1se)),
  prevalence = prev
)
fwrite(metrics_dt, OUT_METRICS)
cat("\nWrote metrics to:", OUT_METRICS, "\n")
print(metrics_dt)

# --- Test predictions ---
pred_dt <- data.table(
  gene_symbol     = test_data$gene_symbol,
  y_true          = y_test,
  prob_lambda_min = prob_min,
  prob_lambda_1se = prob_1se
)
setorder(pred_dt, -prob_lambda_1se)
pred_dt[, rank := .I]
fwrite(pred_dt, OUT_PRED)
cat("Wrote test predictions to:", OUT_PRED, "\n")

# --- Coefficients at lambda.1se ---
coef_sparse <- coef(cvfit, s = "lambda.1se")
coef_mat    <- as.matrix(coef_sparse)
coef_dt <- data.table(
  term = rownames(coef_mat),
  beta = as.numeric(coef_mat[, 1])
)
coef_dt <- coef_dt[term != "(Intercept)"]
coef_dt <- coef_dt[order(-abs(beta))]
fwrite(coef_dt, OUT_COEF)
cat("Wrote coefficients to:", OUT_COEF, "\n")
