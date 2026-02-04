# CRISPR-free Essentiality

**Predicting gene essentiality from unperturbed biological priors,
without training on CRISPR perturbation data.**

------------------------------------------------------------------------

## Overview

Large-scale CRISPR knockout screens have enabled systematic
identification of essential genes, but remain costly, context-limited,
and infeasible in many biological settings. This project explores
whether **gene essentiality can be predicted without using CRISPR data
during model training**, relying instead on unperturbed molecular and
evolutionary signals.

We construct a gene-level feature table integrating: - Baseline gene
expression - Human genetic constraint - Gene redundancy -
Protein–protein interaction network topology

Models are trained exclusively on these **non-CRISPR features** and
evaluated against independent CRISPR-based benchmarks.

------------------------------------------------------------------------

## Scientific Question

> *Is gene essentiality partially encoded in unperturbed biological
> state and long-term evolutionary constraint, prior to experimental
> perturbation?*

This project treats essentiality as a **latent biological property**
that can be inferred without exposure to perturbation labels.

------------------------------------------------------------------------

## Data Sources (raw data not redistributed)

Raw datasets are not included in this repository due to size and
licensing restrictions. Scripts are provided to reproduce all derived
features.

<table>
<thead>
<tr>
<th>Source</th>
<th>Purpose</th>
</tr>
</thead>
<tbody>
<tr>
<td>DepMap Chronos</td>
<td>Evaluation-only gene essentiality labels</td>
</tr>
<tr>
<td>CCLE bulk RNA-seq</td>
<td>Baseline expression features</td>
</tr>
<tr>
<td>gnomAD v2.1.1</td>
<td>Human genetic constraint (LOEUF)</td>
</tr>
<tr>
<td>Ensembl</td>
<td>Paralog relationships</td>
</tr>
<tr>
<td>STRING</td>
<td>Protein–protein interaction network</td>
</tr>
</tbody>
</table>

------------------------------------------------------------------------

## Feature Engineering

The current milestone focuses on building a **CRISPR-free feature
table** at the gene level.

### Feature blocks

<table>
<thead>
<tr>
<th>Biological scale</th>
<th>Features</th>
</tr>
</thead>
<tbody>
<tr>
<td>Cell-state</td>
<td>Mean and variance of CCLE expression</td>
</tr>
<tr>
<td>Organism-level</td>
<td>gnomAD LOEUF constraint</td>
</tr>
<tr>
<td>Redundancy</td>
<td>Paralog count (Ensembl)</td>
</tr>
<tr>
<td>Systems-level</td>
<td>PPI network degree (STRING)</td>
</tr>
</tbody>
</table>

### Ground truth (evaluation only)

- Binary essentiality labels derived from median Chronos scores
- Threshold: Chronos &lt; −0.5
- Labels are **never used during feature construction**

------------------------------------------------------------------------

## Baseline modeling results

Binary essentiality labels are derived from DepMap Chronos medians
(evaluation only)

**Performance (single 80/20 split, seed = 42):**

<table>
<thead>
<tr>
<th>Model</th>
<th style="text-align: right;">AUROC</th>
<th style="text-align: right;">AUPRC</th>
</tr>
</thead>
<tbody>
<tr>
<td>Expression only</td>
<td style="text-align: right;">0.847</td>
<td style="text-align: right;">0.389</td>
</tr>
<tr>
<td>+ LOEUF</td>
<td style="text-align: right;">0.833</td>
<td style="text-align: right;">0.397</td>
</tr>
<tr>
<td>+ Paralog</td>
<td style="text-align: right;">0.839</td>
<td style="text-align: right;">0.424</td>
</tr>
<tr>
<td>Full (expr + LOEUF + paralogs + network)</td>
<td style="text-align: right;"><strong>0.889</strong></td>
<td style="text-align: right;"><strong>0.512</strong></td>
</tr>
</tbody>
</table>

Prevalence baseline AUPRC: **0.092**

#### Multi-Seeding

To assess robustness, the train/test split of 80/20 was repeated across
the data set to indicate that predictive performance isn’t driven by a
particular set of genes. This shows that:

- Each feature captures a partial signal, but a no single feature can
  explain essentiality on it’s own

- The full model performs consistently best across seeds

  - Expression works reasonably well, but not sufficient enough by
    itself

  - LOEUF and Paralog count adds signals

  - **Essentiality is multifactorial**

- The model achieves an average AUPRC that is about five-fold higher
  than the prevalence baseline, indicating enrichment of true essential
  genes beyond random expectation.

Please see figures for more information-
`results/tables/metrics_v1_multiseed_summary.csv` -
`results/tables/metrics_v1_multiseed.csv`

### Interpretation (standardized coefficients)

A Standardized logistic regression coefficient test was done (per 1 SD
increase):

- Mean expression and network degree are the strongest positive
  predictors of essentiality.
- Paralog counts and expression variance are strong negative predictors
- LOEUF contributes a moderate negative effect consistent with
  evolutionary constraint

Please see figures for more information: -
`results/figures/pr_curves.png` -
`results/figures/coefficients_full_glm.csv`

### Logistic Regression Training

Trained a logistic regression model on CRISPR-derived essentiality
labels, using CRISPR-free features (expression, LOEUF, paralogs, PPI
network degree)

Artifacts written by the training script:

\- `results/tables/logreg_training_summary.csv` (AUROC/AUPRC,
prevalence, split info)

\- `results/tables/logreg_test_predictions.csv` (per-gene predicted
essentiality probabilities on the held-out test set)

\- `results/tables/logreg_coefficients.csv` (model coefficients;
interpretability)
