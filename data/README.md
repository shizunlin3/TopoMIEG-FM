# Data sources and integrity

Large third-party datasets are **not redistributed** in this repository. The analysis is bound to exact public resources and SHA256 checksums in `reported_asset_hashes.tsv`.

## Final K562 baseline

The final manuscript uses **ENCODE K562 total RNA-seq experiment ENCSR885DVH**, replicates **ENCFF424CXV** and **ENCFF073NHK**. The materialized `baseline_expression.tsv` used for the frozen predictions has SHA256 `beb9fe6f...e27c63` and 15,360 unique HGNC symbols.

The two downloaded ENCODE replicate TSV hashes recorded at execution were:

- ENCFF424CXV.tsv: `8ed8c904a94eb8ca5de825bc6df406aba4823285256c99d335ec7b103ee6aac8`
- ENCFF073NHK.tsv: `937c049f16e2ee97f7ca315f80fd46a83ae2799698b0653540031eb812ad570a`

## Regulatory graph

DoRothEA human signed TF-target relations, version 1.23.0, confidence A/B only. The exact exported `dorothea_AB.tsv` used for prediction has SHA256 `20d0a16f...9753549e`.

## Held-out perturbation responses

The held-out dataset is the scPerturb-harmonized Replogle et al. K562 genome-scale CRISPRi object `ReplogleWeissman2022_K562_gwps.h5ad`. The frozen object has SHA256 `d3269d9c...8278f104`, 1,989,578 cells and 8,248 genes.

Public sources are described in the manuscript and README. The original perturbation study is Replogle et al., Cell 2022, DOI 10.1016/j.cell.2022.05.013; the harmonized resource is scPerturb (Peidli et al., Nature Methods 2024).

## Fail-closed policy

Release scripts check required columns and, when `--expected-sha256` is supplied, refuse to continue if an input hash differs. This prevents silently analyzing a different dataset revision.
