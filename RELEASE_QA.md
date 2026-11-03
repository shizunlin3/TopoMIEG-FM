# Public release QA

Release date: 2026-11-03

- Final manuscript identity aligned to *When Does a Regulatory Graph Become a Mechanism? Prospective Perturbation Testing with TopoMIEG-FM*.
- Final real-data baseline: ENCODE ENCSR885DVH, replicates ENCFF424CXV and ENCFF073NHK.
- Historical E-MTAB acquisition runner isolated under `provenance/historical_pre_final_baseline/` and marked non-release.
- BRIDGE-01 frozen runner/contract hashes: PASS.
- BRIDGE-02 frozen runner/contract hashes: PASS.
- BRIDGE-01 preflight: PASS.
- BRIDGE-02 preflight: PASS.
- Final Supplementary Table S7 headline values: PASS.
- Real-data prediction and scoring release scripts: Python compile PASS; toy prediction/scoring smoke test PASS.
- `materialize_outcomes.py`: Python compile PASS. Runtime requires the environment in `environment.yml` / `requirements.txt` and the 8.8 GB public H5AD; the large H5AD is intentionally not bundled.
- No public repository URL is embedded in this release package itself; that URL is inserted into the manuscript only after remote publication.
