#!/bin/bash
set -e

VERSION=${1:-v1}

echo "=== CRISPR-free Essentiality Pipeline ==="
echo "VERSION : $VERSION"
echo "Started : $(date)"
echo ""

run_step() {
  local step=$1
  local script=$2
  shift 2
  echo "--- [$step] $script $* ---"
  Rscript "R/$script" "$@"
  echo ""
}

run_step 1 01_build_gene_features.R
run_step 2 02_baseline_models.R          "$VERSION"
run_step 3 03_baseline_models_multiseed.R "$VERSION"
run_step 4 04_training_logreg.R          "$VERSION"
run_step 5 05_train_glmnet.R
run_step 6 06_final_glmnet_coeff.R
run_step 7 07_interpret_coefficients.R
run_step 8 08_plot_pr_curves.R

echo "=== Pipeline complete ==="
echo "Finished : $(date)"
echo "Results  : results/tables/ and results/figures/"
echo ""
echo "Versioned outputs written with prefix: $VERSION"
