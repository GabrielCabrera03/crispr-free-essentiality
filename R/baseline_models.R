
suppressPackageStartupMessages({
  library(data.table)
  library(pROC)
  library(PRROC)
})

#Paths

IN_FILE <- "data/processed/gene_features_v1.tsv.gz"
OUT_METRICS <-"results/tables/metrics_v1.csv"

dir.create("results", showWarnings = FALSE)
dir.create("results/tables", showWarnings = FALSE, recursive = TRUE)


#Load

dt <- fread(IN_FILE)

# Keep rows with labels only
dt <- dt[!is.na(Essential_label) & !is.na(Chronos_median)]

# Some features may be NA. Let's impute LOEUF simply as median impute and keep a missingness flag
dt[, loeuf_missing := as.integer(is.na(loeuf))]
dt[is.na(loeuf), loeuf := median(loeuf, na.rn = TRUE)]

#making sure response is an integer 0/1
dt[, Essential_label := as.integer(Essential_label)]
dt <- dt[Essential_label %in% c(0L, 1L)]

#TRAIN/TEST Split
# columns used anywhere in models
model_cols <- c("mean_expr","var_expr","loeuf","loeuf_missing","log_paralog_count","log_ppi_degree")

# create missingness flags and simple imputations
for (cname in model_cols) {
  if (!cname %in% names(dt)) next
  if (is.numeric(dt[[cname]])) {
    # replace +/-Inf with NA
    dt[!is.finite(get(cname)), (cname) := NA_real_]
  }
}

# LOEUF already imputed; do the same for other numeric features just in case
impute_median <- function(x) { x[is.na(x)] <- median(x, na.rm=TRUE); x }

dt[, mean_expr := impute_median(mean_expr)]
dt[, var_expr := impute_median(var_expr)]
dt[, log_paralog_count := impute_median(log_paralog_count)]
dt[, log_ppi_degree := impute_median(log_ppi_degree)]

set.seed(42)
n <- nrow(dt)
train_idx <- sample(seq_len(n), size = floor(0.8 * n))

train <- dt[train_idx]
test <- dt[-train_idx]

#  class weights to handle imbalance
w_pos <- sum(train$Essential_label == 0) / sum(train$Essential_label == 1)
train[, w := ifelse(Essential_label == 1L, w_pos, 1)]
stopifnot(!anyNA(train$w))

# 
# Evaluate AUROC + AUPRC
# 

eval_metrics <- function(y_true, y_prob) {
  # keep only complete, finite pairs
  ok <- is.finite(y_prob) & !is.na(y_true)
  y_true <- y_true[ok]
  y_prob <- y_prob[ok]
  
  roc_obj <- pROC::roc(y_true, y_prob, quiet = TRUE)
  auroc <- as.numeric(pROC::auc(roc_obj))
  
  # PRROC expects no NA and uses class0 = positives, class1 = negatives
  pos <- y_prob[y_true == 1]
  neg <- y_prob[y_true == 0]
  pos <- pos[is.finite(pos)]
  neg <- neg[is.finite(neg)]
  
  pr <- PRROC::pr.curve(scores.class0 = pos, scores.class1 = neg, curve = FALSE)
  auprc <- pr$auc.integral
  
  list(auroc = auroc, auprc = auprc)
}

# 
# Define ablation feature sets
# 

feature_sets <- list(
  expr_only = c("mean_expr", "var_expr"),
  expr_loeuf = c("mean_expr", "var_expr", "loeuf", "loeuf_missing"),
  expr_loeuf_paralog = c("mean_expr", "var_expr", "loeuf", "loeuf_missing", "log_paralog_count"),
  full_v1 = c("mean_expr", "var_expr", "loeuf", "loeuf_missing", "log_paralog_count", "log_ppi_degree")
)

# 
# Fit logisitc regression for each feature set
# 

metrics_out <- data.table(
  model = character(),
  auroc = numeric(),
  auprc = numeric()
)

for (nm in names(feature_sets)) {
  feats <- feature_sets[[nm]]
  
  # Build formula
  fml <- as.formula(paste("Essential_label ~", paste(feats, collapse = " + ")))
  
  fit <- glm(
    fml,
    data = train,
    family = binomial(),
    weights = train$w
  )
  
  pred <- predict(fit, newdata = test, type = "response")
  m <- eval_metrics(test$Essential_label, pred)
  
  metrics_out <- rbind(
    metrics_out,
    data.table(model = nm, auroc = m$auroc, auprc = m$auprc)
  )
}

# Baseline prevalence for PR comparison
prevalence <- mean(test$Essential_label)
metrics_out[, prevalence := prevalence]

fwrite(metrics_out, OUT_METRICS)
print(metrics_out)

cat("\nWrote metrics to:", OUT_METRICS, "\n")
