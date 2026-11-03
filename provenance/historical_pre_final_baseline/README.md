# Historical pre-final runner — provenance only

This directory is retained for audit provenance. **Do not use `RECOMB-SA5R9R2_runner.py` as the release entry point.**

That historical executable was created before the final manuscript baseline provenance was consolidated and contains an E-MTAB-11264 acquisition path. The final manuscript and final release use ENCODE K562 total RNA-seq **ENCSR885DVH**, replicates **ENCFF424CXV** and **ENCFF073NHK**. The release implementation is `src/realdata/freeze_predictions.py`.

The historical prediction mathematics, graph construction, blind source-universe rule, perturbation magnitude, rewiring scheme and numerical thresholds are preserved here to make the development lineage auditable.
