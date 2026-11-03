# TopoMIEG-FM

**Code and frozen analysis provenance for:**
*When Does a Regulatory Graph Become a Mechanism? Prospective Perturbation Testing with TopoMIEG-FM* (RECOMB 2027 submission).

TopoMIEG-FM evaluates a candidate signed regulatory graph by asking whether a perturbation-induced state transition predicted from that graph agrees with an independent perturbation experiment. The release keeps candidate-graph construction, mathematical stability, perturbation-response reproducibility, and predictive agreement as distinct quantities.

## Headline frozen results

- Controlled synthetic BRIDGE-01: generating mechanism preferred to the deliberately higher-confidence wrong alternative in **128/128** worlds; median paired signed-response difference **0.09827** (95% bootstrap interval **0.09435–0.10254**).
- Controlled synthetic BRIDGE-02: generating mechanism preferred to a structurally matched hard decoy in **128/128** worlds; median paired difference **0.01654** (95% interval **0.01624–0.01688**); context swap reversed the ranking in **128/128** worlds.
- Prospective K562 test: **39** frozen source models; median held-out Spearman correlation **−0.0455** (95% source-bootstrap interval **−0.1713 to 0.0803**). Equilibrium propagation did not outperform the one-step comparator in aggregate.
- CP10K+log1p sensitivity: median held-out Spearman **−0.0603** (95% interval **−0.1169 to 0.0490**).
- The directly measured source gene decreased in all **29/29** measurable perturbations, while downstream guide-to-guide reproducibility was weak.

These values are reproduced in `results/tables/`, exported directly from the final Supplementary Information.

## Repository layout

```text
config/                  final release constants and reported endpoints
data/                    public-source instructions, frozen source list, SHA256 registry
src/synthetic/bridge01/  exact prospectively frozen BRIDGE-01 executable and contract
src/synthetic/bridge02/  exact prospectively frozen BRIDGE-02 executable and contract
src/realdata/            manuscript-facing K562 prediction, outcome, and scoring pipeline
results/tables/           Supplementary Tables S1-S7 exported verbatim
manuscript/               final submission PDF and Supplementary Information
provenance/               formal execution logs and superseded historical runner lineage
scripts/verify_release.py release integrity check
```

## Quick verification

Python 3.11 is recommended.

```bash
python -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
python scripts/verify_release.py
python src/synthetic/bridge01/LGE_BRIDGE_01_runner.py preflight
python src/synthetic/bridge02/LGE_BRIDGE_02_runner.py preflight
```

The two synthetic runners and contracts are copied byte-for-byte from their prospectively frozen versions; `verify_release.py` checks their original SHA256 identities.

## Real-data reproduction

### 1. Obtain the public inputs

Large third-party files are not redistributed. See `data/README.md` and `data/reported_asset_hashes.tsv`.

The final K562 baseline is **ENCODE ENCSR885DVH**, replicates **ENCFF424CXV** and **ENCFF073NHK**. DoRothEA is version 1.23.0, human confidence A/B. The held-out perturbation dataset is the scPerturb-harmonized Replogle K562 genome-wide CRISPRi object.

### 2. Freeze predictions before reading outcomes

Prepare the exact `dorothea_AB.tsv` and `baseline_expression.tsv`, verify their hashes, then run:

```bash
python src/realdata/freeze_predictions.py \
  --dorothea data/work/dorothea_AB.tsv \
  --baseline data/work/baseline_expression.tsv \
  --outdir predictions
```

The release implementation contains **no Perturb-seq outcome input** and selects all baseline-expressed DoRothEA A/B TFs with at least 20 baseline-expressed A/B targets. It then freezes native equilibrium predictions, one-step predictions, and 100 signed-degree-preserving rewires per source.

### 3. Materialize held-out responses

Only after prediction freezing, materialize the perturbation responses. The final aggregation is batch-specific guide mean minus same-batch pooled explicit-control mean, then equal-weight averaging across batches within guide and guides within source.

```bash
python src/realdata/materialize_outcomes.py \
  --h5ad data/raw/ReplogleWeissman2022_K562_gwps.h5ad \
  --sources data/frozen_k562_sources.txt \
  --outdir outcomes/raw
```

The script validates the H5AD SHA256 by default and fails closed if metadata columns cannot be identified unambiguously. Column names can be supplied explicitly for a dataset revision.

For the representation sensitivity, rerun with `--representation cp10k-log1p`.

### 4. Score the frozen predictions

```bash
python src/realdata/score_realdata.py \
  --predictions predictions/prediction_vectors.tsv.gz \
  --rewired predictions/rewired_predictions.tsv.gz \
  --observed outcomes/raw/observed_response_vectors.tsv.gz \
  --outdir results/reproduced_raw
```

The primary statistic is source-level Spearman correlation between predicted and observed downstream changes; aggregation is the median across sources; uncertainty uses 2,000 source-bootstrap replicates with seed 20261108.

## Important provenance note

`provenance/historical_pre_final_baseline/` contains an earlier fail-closed runner that predates the final ENCODE baseline consolidation. It is **not** the release entry point. It is retained only so reviewers can audit development lineage. The final manuscript-facing baseline is ENCODE ENCSR885DVH as stated above and in the paper.

## Reproducibility boundary

Raw public datasets are governed by their original providers and are not vendored here. The repository records exact identifiers and SHA256 hashes for the files used in the reported analysis. Gzip container hashes can depend on serialization metadata when files are regenerated; scientific reproduction should therefore verify input identities, row-level outputs, and reported statistics, while the listed historical hashes identify the exact frozen artifacts used for the manuscript.

## Authors

Zunlin Shi, Zhi Li, Lizhou He, Fan Yang. Correspondence: Fan Yang.

## License

No open-source license has been assigned in this submission release. The repository is public for transparency and reproducibility; reuse rights remain with the authors unless a license is added later.
