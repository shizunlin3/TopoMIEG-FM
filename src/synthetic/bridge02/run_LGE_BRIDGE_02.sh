#!/usr/bin/env bash
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
CONTRACT="$HERE/LGE_BRIDGE_02_CONTRACT.json"
RUNNER="$HERE/LGE_BRIDGE_02_runner.py"
EXPECTED_CONTRACT_SHA="76ff7888256a588d2b28b3c5ac0b33812c1287bfbaac3f47d24ef6bcf2a55cef"
EXPECTED_RUNNER_SHA="c1f435dafc9b00e875919dc325988af8436f00ddc20ceecb700b5e5a6c65c07a"
MODE="${1:-preflight}"
sha() { sha256sum "$1" | awk '{print $1}'; }
[[ "$(sha "$CONTRACT")" == "$EXPECTED_CONTRACT_SHA" ]] || { echo "CONTRACT_SHA_MISMATCH" >&2; exit 41; }
[[ "$(sha "$RUNNER")" == "$EXPECTED_RUNNER_SHA" ]] || { echo "RUNNER_SHA_MISMATCH" >&2; exit 42; }
case "$MODE" in
  preflight) exec python "$RUNNER" preflight ;;
  execute) exec python "$RUNNER" execute --outdir "$HERE/formal_results" ;;
  *) echo "usage: $0 [preflight|execute]" >&2; exit 2 ;;
esac
