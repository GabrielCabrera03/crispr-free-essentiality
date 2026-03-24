<hr>
title: “README” author: “Gabriel Cabrera” date: “2026-01-30” output:
html\_document editor\_options: markdown: wrap: 72
<hr>

# CRISPR-free Essentiality

> **Predicting gene essentiality from unperturbed biological priors —
> without training on CRISPR perturbation data.**

Gene essentiality can be predicted with **AUROC 0.889** using only
evolutionary and expression features, with no CRISPR training data. This
represents a **~5× enrichment** over the prevalence baseline (AUPRC:
0.512 vs. 0.092).

<hr>

## Table of Contents

- [Overview](#overview)
- [Scientific Question](#scientific-question)
- [Repository Structure](#repository-structure)
- [Data Sources](#data-sources)
- [Feature Engineering](#feature-engineering)
- [Results](#results)
- [Interpretation](#interpretation)
- [Reproducing This Work](#reproducing-this-work)
- [Limitations & Future Directions](#limitations--future-directions)

<hr>

## Overview

Large-scale CRISPR knockout screens have enabled systematic
identification of essential genes, but remain costly, context-limited,
and in-feasible in many biological settings. This project explores
whether **gene essentiality can be predicted without using CRISPR data
during model training**, relying instead on unperturbed molecular and
evolutionary signals.

We construct a gene-level feature table integrating: | Biological Scale
| Features | Cell-state | Mean and variance of CCLE bulk RNA-seq
expression | | Organism-level | gnomAD LOEUF constraint (intolerance to
loss-of-function) | | Redundancy | Paralog count (Ensembl) | |
Systems-level | PPI network degree (STRING) |

Models are trained **exclusively** on these non-CRISPR features and
evaluated against independent CRISPR-based benchmarks (DepMap Chronos).

<hr>

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
features from the sources below.

<table>
<colgroup>
<col style="width: 33%" />
<col style="width: 33%" />
<col style="width: 33%" />
</colgroup>
<thead>
<tr>
<th>Source</th>
<th>Version</th>
<th>Purpose</th>
</tr>
</thead>
<tbody>
<tr>
<td><a href="https://depmap.org/portal/">DepMap Chronos</a></td>
<td>23Q4</td>
<td>Evaluation-only essentiality labels</td>
</tr>
<tr>
<td><a href="https://depmap.org/portal/download/">CCLE RNA-seq</a></td>
<td>—</td>
<td>Baseline expression features</td>
</tr>
<tr>
<td><a href="https://gnomad.broadinstitute.org/">gnomAD</a></td>
<td>v2.1.1</td>
<td>Human genetic constraint (LOEUF)</td>
</tr>
<tr>
<td><a href="https://ensembl.org/">Ensembl</a></td>
<td>GRCh38</td>
<td>Paralog relationships</td>
</tr>
<tr>
<td><a href="https://string-db.org/">STRING</a></td>
<td>v12.0</td>
<td>Protein–protein interaction network</td>
</tr>
</tbody>
</table>

> **Ground truth**: Binary essentiality labels are derived from median
> Chronos scores (threshold: Chronos &lt; −0.5). Labels are used **for
> evaluation only** and are never incorporated into feature
> construction.

## Feature Engineering

The current milestone focuses on building a **CRISPR-free feature
table** at the gene level.

The feature table is constructed at the gene level by integrating
signals across four biological scales:

    Gene
     ├── Cell-state         →  mean expression, expression variance (CCLE)
     ├── Organism-level     →  LOEUF score (gnomAD)
     ├── Redundancy         →  paralog count (Ensembl)
     └── Systems-level      →  PPI network degree (STRING)

### Ground truth (evaluation only)

- Binary essentiality labels derived from median Chronos scores
- Threshold: Chronos &lt; −0.5
- Labels are **never used during feature construction**

------------------------------------------------------------------------

## Baseline modeling results

### Ablation — Single Split (80/20, seed = 42)

<table>
<colgroup>
<col style="width: 33%" />
<col style="width: 33%" />
<col style="width: 33%" />
</colgroup>
<thead>
<tr>
<th>Model</th>
<th>AUROC</th>
<th>AUPRC</th>
</tr>
</thead>
<tbody>
<tr>
<td>Expression only</td>
<td>0.847</td>
<td>0.389</td>
</tr>
<tr>
<td>+ LOEUF</td>
<td>0.833</td>
<td>0.397</td>
</tr>
<tr>
<td>+ Paralog count</td>
<td>0.839</td>
<td>0.424</td>
</tr>
<tr>
<td><strong>Full model</strong> (expr + LOEUF + paralogs + network)</td>
<td><strong>0.889</strong></td>
<td><strong>0.512</strong></td>
</tr>
<tr>
<td>Prevalence baseline</td>
<td>—</td>
<td>0.092</td>
</tr>
</tbody>
</table>

The full model achieves a **~5.6× AUPRC enrichment** over the random
baseline.

### Multi-Seed Robustness

To confirm that performance is not driven by a particular train/test
split, the 80/20 split was repeated across multiple random seeds. Key
findings:

- The **full model performs best across all seeds**, confirming that the
  result is not a lucky split
- **Each feature block contributes independent signal** — no single
  feature is sufficient on its own
- **Expression alone is a strong baseline** but meaningfully improves
  with constraint, redundancy, and network features
- **Essentiality is multifactorial**: the convergence of evolutionary,
  expression, and network signals outperforms any single axis

See `results/tables/metrics_v1_multiseed_summary.csv` for full
multi-seed summary statistics.

<figure>
<img src="results/figures/pr_curves.png" alt="PR Curves" />
<figcaption aria-hidden="true">PR Curves</figcaption>
</figure>

<hr>

## Interpretation

Standardized logistic regression coefficients (per 1 SD increase in each
feature):

<table>
<colgroup>
<col style="width: 33%" />
<col style="width: 33%" />
<col style="width: 33%" />
</colgroup>
<thead>
<tr>
<th>Feature</th>
<th>Direction</th>
<th>Interpretation</th>
</tr>
</thead>
<tbody>
<tr>
<td>Mean expression</td>
<td>↑ Positive</td>
<td>Highly expressed genes are more likely essential</td>
</tr>
<tr>
<td>PPI network degree</td>
<td>↑ Positive</td>
<td>Hub genes in interaction networks tend to be essential</td>
</tr>
<tr>
<td>Expression variance</td>
<td>↓ Negative</td>
<td>Stably expressed genes are more essential than variable ones</td>
</tr>
<tr>
<td>Paralog count</td>
<td>↓ Negative</td>
<td>Genes with redundant paralogs are buffered against lethality</td>
</tr>
<tr>
<td>LOEUF</td>
<td>↓ Negative</td>
<td>Genes under strong evolutionary constraint are more likely
essential</td>
</tr>
</tbody>
</table>

Mean expression and network degree are the **strongest positive
predictors**. Paralog count and expression variance are the **strongest
negative predictors**. LOEUF contributes a moderate negative effect
consistent with evolutionary constraint theory.

See
[`results/figures/coefficients_full_glm.csv`](results/figures/coefficients_full_glm.csv)
for full standardized coefficient values.

### Output Artifacts

<table>
<colgroup>
<col style="width: 50%" />
<col style="width: 50%" />
</colgroup>
<thead>
<tr>
<th>File</th>
<th>Description</th>
</tr>
</thead>
<tbody>
<tr>
<td><code>results/tables/logreg_training_summary.csv</code></td>
<td>AUROC/AUPRC, prevalence, split info</td>
</tr>
<tr>
<td><code>results/tables/logreg_test_predictions.csv</code></td>
<td>Per-gene predicted essentiality probabilities</td>
</tr>
<tr>
<td><code>results/tables/logreg_coefficients.csv</code></td>
<td>Standardized model coefficients</td>
</tr>
<tr>
<td><code>results/tables/metrics_v1_multiseed.csv</code></td>
<td>Per-seed performance metrics</td>
</tr>
<tr>
<td><code>results/tables/metrics_v1_multiseed_summary.csv</code></td>
<td>Summary statistics across seeds</td>
</tr>
</tbody>
</table>

<hr>

## Limitations & Future Directions

### Current Limitations

- **Paralog-aware splitting**: Paralogs of training genes may appear in
  the test set, creating a potential source of information leakage.
  Future splits should be paralog-aware.
- **Cell line aggregation**: Chronos scores are aggregated to a median
  across cell lines, collapsing meaningful context-specificity.
  Per-lineage models may reveal tissue-specific essentiality patterns.
- **Model class**: Logistic regression was chosen for interpretability.
  Non-linear models (gradient boosting, GNNs over the PPI network) may
  capture interaction effects not modeled here.

### Future Directions

- ☐ Paralog-aware train/test splitting
- ☐ Per-lineage essentiality modeling (hematologic vs. solid tumors)
- ☐ Negative control validation (permuted labels, random gene sets)
- ☐ Comparison to published CRISPR-free essentiality baselines
- ☐ Model calibration assessment (Brier score, reliability diagrams)
- ☐ Extension to non-cancer contexts where CRISPR screens are infeasible

<hr>

## Citation

If you use this work, please cite this repository:

    @misc{cabrera2024crisprfree,
      author = {Cabrera, Gabriel},
      title  = {CRISPR-free Essentiality},
      year   = {2024},
      url    = {https://github.com/GabrielCabrera03/crispr-free-essentiality}
    }

<hr>
