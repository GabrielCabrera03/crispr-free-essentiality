#!/usr/bin/env Rscript

# -----------------------------
# Build gene_features_v1
# - CCLE expression (mean/var)
# - gnomAD LOEUF
# - Ensembl paralog count
# - STRING PPI degree
# - DepMap Chronos labels (eval-only)
# -----------------------------

suppressPackageStartupMessages({
  library(data.table)
  library(biomaRt)
})

# Optional: only needed if your CCLE file is GCT and you read it via cmapR
# install.packages("cmapR") if missing
suppressPackageStartupMessages({
  have_cmapR <- requireNamespace("cmapR", quietly = TRUE)
})

# -----------------------------
# CONFIG (edit these paths)
# -----------------------------
RAW_DIR  <- "data/raw"
OUT_DIR  <- "data/processed"

# DepMap Chronos matrix (rows = cell lines, cols = gene symbols)
CHRONOS_CSV <- file.path(RAW_DIR, "CRISPR_(DepMap_Public_25Q3+Score,_Chronos)_subsetted.csv")

# CCLE expression (GCT) you used
CCLE_GCT <- file.path(RAW_DIR, "CCLE_DepMap_18q3_RNAseq_RPKM_20180718.gct")

# gnomAD LOEUF by gene
GNOMAD_LOEUF_BGZ <- file.path(RAW_DIR, "gnomad.v2.1.1.lof_metrics.by_gene.txt.bgz")

# STRING full network (human 9606)
STRING_LINKS_GZ <- file.path(RAW_DIR, "9606.protein.links.v12.0.txt.gz")
# If you used physical network instead:
# STRING_LINKS_GZ <- file.path(RAW_DIR, "9606.protein.physical.links.v12.0.txt.gz")

# Network confidence threshold (STRING scores are 0–1000)
STRING_SCORE_MIN <- 700  # common high-confidence threshold

# Output files
OUT_TSV_GZ <- file.path(OUT_DIR, "gene_features_v1.tsv.gz")
OUT_RDS    <- file.path(OUT_DIR, "gene_features_v1.rds")

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

message("=== Build gene_features_v1 ===")

# -----------------------------
# Helpers
# -----------------------------
stop_if_missing <- function(path) {
  if (!file.exists(path)) stop("Missing file: ", path)
}

as_numeric_safe <- function(x) suppressWarnings(as.numeric(x))

# -----------------------------
# 1) DepMap Chronos -> gene-level labels (eval-only)
# -----------------------------
stop_if_missing(CHRONOS_CSV)
message("[1/5] Reading DepMap Chronos matrix: ", CHRONOS_CSV)

chronos <- fread(CHRONOS_CSV)

# Expect first column is cell line id and remaining columns are gene symbols
# If your file has rownames already saved, adjust here.
# We'll assume column 1 is the cell line identifier.
cell_id_col <- names(chronos)[1]
gene_cols <- setdiff(names(chronos), cell_id_col)

# Compute per-gene median across cell lines
chronos_medians <- data.table(
  gene_symbol = gene_cols,
  Chronos_median = vapply(gene_cols, function(g) median(chronos[[g]], na.rm = TRUE), numeric(1))
)

chronos_medians[, Essential_label := ifelse(Chronos_median < -0.5, 1L, 0L)]

message("  Chronos genes: ", nrow(chronos_medians),
        " | essential: ", sum(chronos_medians$Essential_label == 1, na.rm = TRUE))

# -----------------------------
# 2) CCLE expression -> mean/var features (gene-level)
# -----------------------------
stop_if_missing(CCLE_GCT)
message("[2/5] Reading CCLE expression GCT: ", CCLE_GCT)

if (!have_cmapR) {
  stop("Package 'cmapR' not installed. Install it with install.packages('cmapR').")
}

gct <- cmapR::parse_gctx(CCLE_GCT)  # parse_gctx works for gct/gctx in cmapR
expr_mat <- gct@mat
rownames(expr_mat) <- gct@rid
colnames(expr_mat) <- gct@cid

# Strip Ensembl version suffix
ensembl_id <- sub("\\..*$", "", rownames(expr_mat))

mean_expr <- rowMeans(expr_mat, na.rm = TRUE)
var_expr  <- apply(expr_mat, 1, var, na.rm = TRUE)

expr_features <- data.table(
  ensembl_gene_id = ensembl_id,
  mean_expr = mean_expr,
  var_expr  = var_expr
)

# -----------------------------
# 3) Map Ensembl gene IDs -> HGNC symbols and merge with Chronos labels
# -----------------------------
message("[3/5] Mapping Ensembl gene IDs -> HGNC symbols (biomaRt)")

mart <- biomaRt::useMart("ensembl", dataset = "hsapiens_gene_ensembl")

ens2sym <- getBM(
  attributes = c("ensembl_gene_id", "hgnc_symbol"),
  filters = "ensembl_gene_id",
  values = unique(expr_features$ensembl_gene_id),
  mart = mart
)
ens2sym <- ens2sym[ens2sym$hgnc_symbol != "", ]
ens2sym <- ens2sym[!duplicated(ens2sym$ensembl_gene_id), ]

expr_mapped <- merge(expr_features, as.data.table(ens2sym),
                     by = "ensembl_gene_id", all.x = FALSE)

setnames(expr_mapped, "hgnc_symbol", "gene_symbol")

# Merge features with labels; drop genes lacking Chronos labels
feature_table <- merge(expr_mapped, chronos_medians, by = "gene_symbol", all.x = FALSE)

message("  After merge with Chronos: ", nrow(feature_table), " genes")

# -----------------------------
# 4) Add LOEUF (gnomAD)
# -----------------------------
stop_if_missing(GNOMAD_LOEUF_BGZ)
message("[4/5] Reading gnomAD LOEUF: ", GNOMAD_LOEUF_BGZ)

gn <- fread(GNOMAD_LOEUF_BGZ)

# Detect likely symbol + loeuf columns robustly
sym_candidates <- intersect(names(gn), c("gene", "gene_name", "hgnc", "hgnc_symbol"))
if (length(sym_candidates) == 0) stop("Could not find gene symbol column in gnomAD file.")
symbol_col <- sym_candidates[1]

loeuf_candidates <- names(gn)[grepl("loeuf|oe_lof_upper|lof\\.oe_ci\\.upper|oe_lof_upper", names(gn), ignore.case = TRUE)]
if (length(loeuf_candidates) == 0) stop("Could not find LOEUF column in gnomAD file.")
loeuf_col <- loeuf_candidates[1]

gnomad_loeuf <- gn[, .(
  gene_symbol = get(symbol_col),
  loeuf = as_numeric_safe(get(loeuf_col))
)]
gnomad_loeuf <- gnomad_loeuf[gene_symbol != ""]
gnomad_loeuf <- unique(gnomad_loeuf, by = "gene_symbol")

feature_table <- merge(feature_table, gnomad_loeuf, by = "gene_symbol", all.x = TRUE)

message("  LOEUF missing: ", sum(is.na(feature_table$loeuf)), " / ", nrow(feature_table))

# -----------------------------
# 5) Add paralog count (Ensembl)
# -----------------------------
message("[5/5] Computing paralog_count via Ensembl (biomaRt)")

# Map symbols -> Ensembl gene IDs for our current table
sym2ens <- getBM(
  attributes = c("hgnc_symbol", "ensembl_gene_id"),
  filters = "hgnc_symbol",
  values = unique(feature_table$gene_symbol),
  mart = mart
)
sym2ens <- sym2ens[sym2ens$hgnc_symbol != "" & sym2ens$ensembl_gene_id != "", ]

# --- Query paralogs in Ensembl-ID space ---
# Attribute names vary across Ensembl releases, so detect dynamically
ens2par <- getBM(
  attributes = c("ensembl_gene_id", "hsapiens_paralog_ensembl_gene"),
  filters = "ensembl_gene_id",
  values = unique(sym2ens$ensembl_gene_id),
  mart = mart
)

# Identify paralog column robustly
paralog_col_candidates <- names(ens2par)[grepl("paralog", names(ens2par), ignore.case = TRUE)]
if (length(paralog_col_candidates) == 0) {
  stop(
    "Could not find paralog column. Columns returned by biomaRt are: ",
    paste(names(ens2par), collapse = ", ")
  )
}

paralog_col <- setdiff(paralog_col_candidates, "ensembl_gene_id")[1]

# Clean empty + self-matches
ens2par <- ens2par[
  ens2par[[paralog_col]] != "" &
    ens2par$ensembl_gene_id != ens2par[[paralog_col]],
]

# Count unique paralogs per Ensembl gene
par_counts_ens <- as.data.table(ens2par)[
  , .(paralog_count = uniqueN(get(paralog_col))), by = ensembl_gene_id
]

# Map back to gene symbols
par_counts_sym <- merge(
  as.data.table(sym2ens),
  par_counts_ens,
  by = "ensembl_gene_id",
  all.x = FALSE
)

# If a symbol maps to multiple Ensembl IDs, keep the maximum count (conservative)
par_counts_sym <- par_counts_sym[
  , .(paralog_count = max(paralog_count, na.rm = TRUE)), by = hgnc_symbol
]
setnames(par_counts_sym, "hgnc_symbol", "gene_symbol")

# Merge into feature table
feature_table <- merge(
  feature_table,
  par_counts_sym,
  by = "gene_symbol",
  all.x = TRUE
)

feature_table[is.na(paralog_count), paralog_count := 0L]
feature_table[, log_paralog_count := log1p(paralog_count)]
# -----------------------------
# 6) Add STRING network degree (full network)
# -----------------------------
stop_if_missing(STRING_LINKS_GZ)
message("[6/5] (extra) Computing STRING degree from: ", STRING_LINKS_GZ)
message("  Using STRING_SCORE_MIN = ", STRING_SCORE_MIN)

edges <- fread(STRING_LINKS_GZ)
# expected columns: protein1 protein2 combined_score (names may vary slightly)
score_col <- if ("combined_score" %in% names(edges)) "combined_score" else names(edges)[3]
setnames(edges, score_col, "combined_score")

edges <- edges[combined_score >= STRING_SCORE_MIN]

# strip species prefix "9606."
edges[, protein1 := sub("^9606\\.", "", protein1)]
edges[, protein2 := sub("^9606\\.", "", protein2)]

all_ensp <- unique(c(edges$protein1, edges$protein2))

# Map Ensembl peptide ID -> gene symbol
ensp2sym <- getBM(
  attributes = c("ensembl_peptide_id", "hgnc_symbol"),
  filters = "ensembl_peptide_id",
  values = all_ensp,
  mart = mart
)
ensp2sym <- ensp2sym[ensp2sym$hgnc_symbol != "", ]

# Join endpoints to symbols
e1 <- merge(edges, as.data.table(ensp2sym), by.x = "protein1", by.y = "ensembl_peptide_id", all.x = FALSE)
setnames(e1, "hgnc_symbol", "gene1")
e2 <- merge(edges, as.data.table(ensp2sym), by.x = "protein2", by.y = "ensembl_peptide_id", all.x = FALSE)
setnames(e2, "hgnc_symbol", "gene2")

deg1 <- e1[, .N, by = gene1]; setnames(deg1, c("gene1","N"), c("gene_symbol","deg"))
deg2 <- e2[, .N, by = gene2]; setnames(deg2, c("gene2","N"), c("gene_symbol","deg"))

ppi_degree <- rbindlist(list(deg1, deg2))[, .(ppi_degree = sum(deg)), by = gene_symbol]
ppi_degree[, log_ppi_degree := log1p(ppi_degree)]

feature_table <- merge(feature_table, ppi_degree, by = "gene_symbol", all.x = TRUE)
feature_table[is.na(ppi_degree), ppi_degree := 0L]
feature_table[is.na(log_ppi_degree), log_ppi_degree := 0]

# -----------------------------
# 7) Final cleanup + save
# -----------------------------
setcolorder(feature_table, c(
  "gene_symbol",
  "mean_expr","var_expr",
  "loeuf",
  "paralog_count","log_paralog_count",
  "ppi_degree","log_ppi_degree",
  "Chronos_median","Essential_label"
))

message("Final rows: ", nrow(feature_table))
message("Essential genes: ", sum(feature_table$Essential_label == 1, na.rm = TRUE))

# Save compressed TSV + RDS
fwrite(feature_table, OUT_TSV_GZ, sep = "\t")
saveRDS(feature_table, OUT_RDS)

message("Wrote: ", OUT_TSV_GZ)
message("Wrote: ", OUT_RDS)

# Optional: quick sanity correlations
message("Sanity: Spearman(mean_expr, Chronos) = ",
        cor(feature_table$mean_expr, feature_table$Chronos_median, method="spearman", use="complete.obs"))
message("Sanity: Spearman(loeuf, Chronos) = ",
        cor(feature_table$loeuf, feature_table$Chronos_median, method="spearman", use="complete.obs"))
message("Sanity: Spearman(log_paralog_count, Chronos) = ",
        cor(feature_table$log_paralog_count, feature_table$Chronos_median, method="spearman", use="complete.obs"))
message("Sanity: Spearman(log_ppi_degree, Chronos) = ",
        cor(feature_table$log_ppi_degree, feature_table$Chronos_median, method="spearman", use="complete.obs"))
