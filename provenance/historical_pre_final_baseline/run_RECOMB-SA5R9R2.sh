#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "$0")" && pwd)"
PY="${PYTHON:-python3}"

"$PY" "$ROOT/RECOMB-SA5R9R2_runner.py" preflight

echo "1) $PY $ROOT/RECOMB-SA5R9R2_runner.py acquire-dorothea"
echo "2) $PY $ROOT/RECOMB-SA5R9R2_runner.py inventory-baseline"
echo "3) preprocess E-MTAB-11264 K562 RNA using BASELINE_PREPROCESSING_CONTRACT.md"
echo "4) $PY $ROOT/RECOMB-SA5R9R2_runner.py register-baseline --tsv /path/to/baseline_expression.tsv"
echo "5) $PY $ROOT/RECOMB-SA5R9R2_runner.py vault-outcome"
echo "6) $PY $ROOT/RECOMB-SA5R9R2_runner.py predict"
echo "7) $PY $ROOT/RECOMB-SA5R9R2_runner.py freeze"
