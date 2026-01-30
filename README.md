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

## Feature Engineering (Week 2 milestone)

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

## Repository Structure
