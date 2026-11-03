#!/usr/bin/env bash
set -euo pipefail
HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
EXPECTED_CONTRACT_SHA="f35834415be43ed26c0300c3778455ca0e27f63dc0b0333165ee4a2c64c3b7c5"
EXPECTED_RUNNER_SHA="1af2fd91924c974722fddc2dc486cb08406b4667b77cfe625f8d95bfdea7b426"
MODE="${1:-preflight}"

actual_contract=$(sha256sum "$HERE/LGE_BRIDGE_01_CONTRACT.json" | awk '{print $1}')
actual_runner=$(sha256sum "$HERE/LGE_BRIDGE_01_runner.py" | awk '{print $1}')
[[ "$actual_contract" == "$EXPECTED_CONTRACT_SHA" ]] || { echo "FAIL contract SHA drift" >&2; exit 2; }
[[ "$actual_runner" == "$EXPECTED_RUNNER_SHA" ]] || { echo "FAIL runner SHA drift" >&2; exit 3; }

case "$MODE" in
  preflight)
    exec python "$HERE/LGE_BRIDGE_01_runner.py" preflight
    ;;
  execute)
    exec python "$HERE/LGE_BRIDGE_01_runner.py" execute --outdir "$HERE/formal_results"
    ;;
  *)
    echo "Usage: $0 {preflight|execute}" >&2
    exit 64
    ;;
esac
