
suppressPackageStartupMessages({
  library(data.table)
  library(glmnet)
})

#Paths

IN_FILE <- "data/processed/gene_features_v1.tsv.gz"
OUT_DIR <- "results/tables"
MODEL_DIR<- "models"

OUT_COEF <- file.path(OUT_DIR, "final_glmnet_coefficients.csv")
OUT_META <- file.path(OUT_DIR, "final_glmnet_model_metadata.csv")
OUT_MODEL <- file.path(MODEL_DIR, "final_glmnet_model.rds")

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(MODEL_DIR, showWarnings = FALSE, recursive = TRUE)

#Model Settings

SEED <- 42
ALPHA <- 0.5
N_FOLDS <- 5
USE_LAMBDA <- "lambda.1se"

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

standardize_inplace <- function(dt, feature_names) {
  # Standardize using dataset-wide stats (since we are fitting final model on all data)
  scaler <- list()
  for (f in feature_names) {
    mu <- mean(dt[[f]], na.rm = TRUE)
    sdv <- sd(dt[[f]], na.rm = TRUE)
    if (!is.finite(sdv) || sdv == 0) sdv <- 1
    dt[[f]] <- (dt[[f]] - mu) / sdv
    scaler[[f]] <- list(mean = mu, sd = sdv)
  }
  scaler
}

#Load Data

set.seed(SEED)
dt <- fread(IN_FILE)

dt <- dt[!is.na(Essential_label)]
dt[, Essential_label := as.integer(Essential_label)]
dt <- dt[Essential_label %in% c(0L, 1L)]

# Impute (final)
for (f in FEATURES) dt[[f]] <- impute_median(dt[[f]])

stopifnot(all(FEATURES %in% names(dt)))
cat("Loaded labeled genes:", nrow(dt), "\n")
cat("Prevalence:", mean(dt$Essential_label), "\n")

# Standardize (final)
scaler <- standardize_inplace(dt, FEATURES)

# Matrices
x <- as.matrix(dt[, ..FEATURES])
y <- dt$Essential_label

# Class weights (final)
w_pos <- sum(y == 0L) / sum(y == 1L)
w <- ifelse(y == 1L, w_pos, 1)

cat("w_pos:", w_pos, "\n")

#Fitting CV glmnet

cvfit <- cv.glmnet(
  x = x, 
  y = y, 
  family = "binomial",
  alpha = ALPHA, 
  weights = w,
  nfolds = N_FOLDS,
  type.measure = "auc"
)

lambda_min <- cvfit$lambda.min
lambda_1se <- cvfit$lambda.1se
lambda_use <- if (USE_LAMBDA == lambda_min) lambda_min else lambda_1se

cat("lambda.min:", lambda_min, "\n")
cat("lambda.1se:", lambda_1se, "\n")
cat("Using:", USE_LAMBDA, "=", lambda_use, "\n")

#Extracting coefficients at chosen lambda

coef_sparse <- coef(cvfit, s = lambda_use)
coef_mat <- as.matrix(coef_sparse)

coef_dt <- data.table(
  term = rownames(coef_mat),
  beta = as.numeric(coef_mat[, 1])
)

# Remove zeros (glmnet shrinks many terms)
coef_dt <- coef_dt[beta != 0]
setorder(coef_dt, -abs(beta))

fwrite(coef_dt, OUT_COEF)
cat("Wrote coefficients:", OUT_COEF, "\n")

#Save Metadata

meta_dt <- data.table(
  model = "glmnet_elasticnet_logistic",
  alpha = ALPHA, 
  n_folds = N_FOLDS, 
  lambda_min = lambda_min,
  lambda_1se = lambda_1se, 
  lambda_used = lambda_use,
  lambda_choice = USE_LAMBDA,
  n_genes = length(y),
  prevalence = mean(y),
  w_pos = w_pos, 
  features = paste(FEATURES, collapse = ";")
  
)
fwrite(meta_dt, OUT_META)
cat("Wrote metadata:", OUT_META, "\n")
