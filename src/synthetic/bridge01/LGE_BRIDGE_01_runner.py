#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Tuple

import numpy as np

HERE = Path(__file__).resolve().parent
CONTRACT = HERE / "LGE_BRIDGE_01_CONTRACT.json"
EXPECTED_CONTRACT_SHA256 = "f35834415be43ed26c0300c3778455ca0e27f63dc0b0333165ee4a2c64c3b7c5"
AUTH_FILE = HERE / "LGE_BRIDGE_01_EXECUTION_AUTHORIZATION.txt"
AUTH_TEXT = "AUTHORIZED_AFTER_RECOMB_SA5R1_FREEZE"


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def load_contract() -> dict:
    got = sha256(CONTRACT)
    if got != EXPECTED_CONTRACT_SHA256:
        raise RuntimeError(f"contract SHA mismatch: {got} != {EXPECTED_CONTRACT_SHA256}")
    return json.loads(CONTRACT.read_text(encoding="utf-8"))


def phi(u: np.ndarray, b: np.ndarray, W: np.ndarray) -> np.ndarray:
    return np.tanh(b + W @ u)


@dataclass(frozen=True)
class EqSolution:
    u: np.ndarray
    residual_inf: float
    rho_eff: float
    iterations: int


def solve_equilibrium(b: np.ndarray, W: np.ndarray, tol: float, max_iter: int) -> EqSolution:
    u = np.tanh(b).astype(np.float64, copy=True)
    for it in range(1, max_iter + 1):
        new = phi(u, b, W)
        if float(np.max(np.abs(new - u))) <= tol:
            u = new
            break
        u = new
    else:
        raise RuntimeError("fixed-point solver did not converge")
    residual = float(np.max(np.abs(u - phi(u, b, W))))
    J = (1.0 - u * u)[:, None] * W
    rho = float(np.max(np.abs(np.linalg.eigvals(J)))) if W.size else 0.0
    return EqSolution(u=u, residual_inf=residual, rho_eff=rho, iterations=it)


def scale_operator(W: np.ndarray, target_norm: float) -> np.ndarray:
    norm = float(np.linalg.norm(W, ord=2))
    if not np.isfinite(norm) or norm <= 0:
        raise RuntimeError("non-positive operator norm")
    return W * (target_norm / norm)


def module_labels(p: int, n_modules: int) -> np.ndarray:
    if p % n_modules != 0:
        raise ValueError("p must be divisible by n_modules")
    return np.repeat(np.arange(n_modules, dtype=int), p // n_modules)


def generate_world(rng: np.random.Generator, wid: int, cfg: dict) -> Dict[str, object]:
    g = cfg["world_generator"]
    p = int(g["p"])
    n_modules = int(g["n_modules"])
    labels = module_labels(p, n_modules)
    W_raw = np.zeros((p, p), dtype=np.float64)  # row=target, column=source

    # Deterministic connectivity backbone: one directed cycle per module.
    for m in range(n_modules):
        nodes = np.where(labels == m)[0]
        for a, b in zip(nodes, np.roll(nodes, -1)):
            mag = rng.uniform(*g["within_abs_weight_range"])
            sign = 1.0 if rng.random() < g["within_positive_probability"] else -1.0
            W_raw[b, a] = sign * mag

    # Additional directed edges.
    for source in range(p):
        for target in range(p):
            if source == target or W_raw[target, source] != 0.0:
                continue
            same = labels[source] == labels[target]
            prob = g["within_edge_probability"] if same else g["cross_edge_probability"]
            if rng.random() < prob:
                lo, hi = g["within_abs_weight_range"] if same else g["cross_abs_weight_range"]
                ppos = g["within_positive_probability"] if same else g["cross_positive_probability"]
                W_raw[target, source] = (1.0 if rng.random() < ppos else -1.0) * rng.uniform(lo, hi)

    # Source must have a within-module neighbor and an external non-neighbor for the decoy.
    candidates = []
    for source in range(p):
        within = [j for j in range(p) if j != source and labels[j] == labels[source] and W_raw[j, source] != 0]
        external_non = [j for j in range(p) if labels[j] != labels[source] and W_raw[j, source] == 0]
        if within and external_non:
            candidates.append((source, within, external_non))
    if not candidates:
        raise RuntimeError("world has no eligible source/decoy pair")

    source, within_targets, external_non = candidates[int(rng.integers(0, len(candidates)))]
    true_target = max(within_targets, key=lambda j: abs(W_raw[j, source]))

    out_degree = np.count_nonzero(W_raw, axis=0)
    in_degree = np.count_nonzero(W_raw, axis=1)
    decoy_target = min(
        external_non,
        key=lambda j: (abs(int(out_degree[j]) - int(out_degree[true_target])) +
                       abs(int(in_degree[j]) - int(in_degree[true_target])), j),
    )

    W = scale_operator(W_raw, float(g["operator_norm_target"]))
    b0 = np.clip(
        rng.normal(0.0, float(g["baseline_field_sd"]), size=p),
        -float(g["baseline_field_clip"]),
        float(g["baseline_field_clip"]),
    )
    delta_abs = float(g["source_perturbation_abs_delta"])
    delta = delta_abs if wid % 2 == 0 else -delta_abs
    edge = float(W[true_target, source])
    if edge == 0:
        raise RuntimeError("true edge vanished")
    s_true = 1.0 if delta * edge > 0 else -1.0

    perm = rng.permutation(p)
    if np.array_equal(perm, np.arange(p)):
        perm = np.roll(perm, 1)
    W_perm = W[np.ix_(perm, perm)]

    return {
        "world_id": wid,
        "W": W,
        "W_perm": W_perm,
        "b0": b0,
        "source": int(source),
        "true_target": int(true_target),
        "decoy_target": int(decoy_target),
        "s_true": float(s_true),
        "delta": float(delta),
        "operator_norm": float(np.linalg.norm(W, ord=2)),
    }


def mechanism_response(base: EqSolution, pert: EqSolution, target: int, s_claim: float, delta: float) -> float:
    return float(s_claim * (pert.u[target] - base.u[target]) / abs(delta))


def admissible(sol: EqSolution, contract: dict) -> bool:
    s = contract["solver"]
    return bool(sol.residual_inf <= float(s["formal_residual_threshold"]) and
                sol.rho_eff <= float(s["formal_effective_spectral_radius_threshold"]))


def wilson_ci(k: int, n: int, z: float = 1.959963984540054) -> Tuple[float, float]:
    if n <= 0:
        return (float("nan"), float("nan"))
    phat = k / n
    den = 1 + z*z/n
    center = (phat + z*z/(2*n)) / den
    half = z * math.sqrt(phat*(1-phat)/n + z*z/(4*n*n)) / den
    return max(0.0, center-half), min(1.0, center+half)


def bootstrap_median_ci(values: np.ndarray, reps: int, seed: int) -> Tuple[float, float]:
    rng = np.random.Generator(np.random.PCG64DXSM(seed))
    n = len(values)
    meds = np.empty(reps, dtype=np.float64)
    for i in range(reps):
        idx = rng.integers(0, n, size=n)
        meds[i] = np.median(values[idx])
    return float(np.quantile(meds, 0.025)), float(np.quantile(meds, 0.975))


def run_formal(contract: dict, outdir: Path) -> dict:
    cfg = contract["world_generator"]
    solver = contract["solver"]
    rng = np.random.Generator(np.random.PCG64DXSM(int(cfg["seed"])))
    rows = []
    for wid in range(int(cfg["n_worlds"])):
        w = generate_world(rng, wid, contract)
        W = w["W"]
        Wp = w["W_perm"]
        b0 = w["b0"]
        b1 = b0.copy(); b1[w["source"]] += w["delta"]

        base = solve_equilibrium(b0, W, float(solver["tolerance"]), int(solver["max_iterations"]))
        pert = solve_equilibrium(b1, W, float(solver["tolerance"]), int(solver["max_iterations"]))
        base_p = solve_equilibrium(b0, Wp, float(solver["tolerance"]), int(solver["max_iterations"]))
        pert_p = solve_equilibrium(b1, Wp, float(solver["tolerance"]), int(solver["max_iterations"]))

        r_true = mechanism_response(base, pert, w["true_target"], w["s_true"], w["delta"])
        r_wrong = mechanism_response(base, pert, w["decoy_target"], w["s_true"], w["delta"])
        r_wrong_sign = mechanism_response(base, pert, w["true_target"], -w["s_true"], w["delta"])
        r_perm = mechanism_response(base_p, pert_p, w["true_target"], w["s_true"], w["delta"])

        row = {
            "world_id": wid,
            "source": w["source"],
            "true_target": w["true_target"],
            "decoy_target": w["decoy_target"],
            "delta": w["delta"],
            "true_claim_sign": w["s_true"],
            "confidence_true": contract["confidence_stress_test"]["true_mechanism_confidence"],
            "confidence_high_wrong": contract["confidence_stress_test"]["high_confidence_wrong_confidence"],
            "R_eq_true": r_true,
            "R_eq_high_wrong": r_wrong,
            "R_eq_wrong_sign": r_wrong_sign,
            "R_eq_permuted_graph": r_perm,
            "paired_delta": r_true-r_wrong,
            "true_beats_high_wrong": int(r_true > r_wrong),
            "true_beats_wrong_sign": int(r_true > r_wrong_sign),
            "true_beats_permuted": int(r_true > r_perm),
            "operator_norm": w["operator_norm"],
            "baseline_residual": base.residual_inf,
            "perturbed_residual": pert.residual_inf,
            "baseline_rho_eff": base.rho_eff,
            "perturbed_rho_eff": pert.rho_eff,
            "matched_admissible": int(admissible(base, contract) and admissible(pert, contract)),
            "permuted_admissible": int(admissible(base_p, contract) and admissible(pert_p, contract)),
        }
        if not all(np.isfinite(float(v)) for k, v in row.items() if k not in {"world_id","source","true_target","decoy_target"}):
            raise RuntimeError(f"nonfinite result in world {wid}")
        if row["true_target"] == row["decoy_target"]:
            raise RuntimeError(f"duplicate true/decoy target in world {wid}")
        rows.append(row)

    outdir.mkdir(parents=True, exist_ok=True)
    tsv = outdir / contract["outputs"]["world_level_tsv"]
    fields = list(rows[0].keys())
    with tsv.open("w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields, delimiter="\t", lineterminator="\n")
        w.writeheader(); w.writerows(rows)

    delta = np.asarray([r["paired_delta"] for r in rows], dtype=np.float64)
    wins = int(sum(r["true_beats_high_wrong"] for r in rows))
    signwins = int(sum(r["true_beats_wrong_sign"] for r in rows))
    permwins = int(sum(r["true_beats_permuted"] for r in rows))
    n = len(rows)
    wl, wu = wilson_ci(wins, n)
    bl, bu = bootstrap_median_ci(
        delta,
        int(contract["primary_endpoints"]["bootstrap"]["replicates"]),
        int(contract["primary_endpoints"]["bootstrap"]["seed"]),
    )
    admiss_frac = float(np.mean([r["matched_admissible"] for r in rows]))
    win_frac = wins/n
    sign_win_frac = signwins/n

    hard = (n == int(cfg["n_worlds"]) and admiss_frac == 1.0)
    go = bool(hard and win_frac >= 0.70 and wl > 0.50 and bl > 0.0 and sign_win_frac >= 0.95)
    summary = {
        "stage": "LGE-BRIDGE-01 formal execution",
        "contract_sha256": EXPECTED_CONTRACT_SHA256,
        "n_worlds": n,
        "matched_admissibility_fraction": admiss_frac,
        "paired_win_fraction_true_over_high_confidence_wrong": win_frac,
        "paired_win_wilson95": [wl, wu],
        "median_paired_delta": float(np.median(delta)),
        "median_paired_delta_bootstrap95": [bl, bu],
        "wrong_sign_sanity_win_fraction": sign_win_frac,
        "matched_true_over_permuted_win_fraction": permwins/n,
        "mean_R_eq_true": float(np.mean([r["R_eq_true"] for r in rows])),
        "mean_R_eq_high_wrong": float(np.mean([r["R_eq_high_wrong"] for r in rows])),
        "maximum_matched_residual": float(max(max(r["baseline_residual"], r["perturbed_residual"]) for r in rows)),
        "maximum_matched_effective_spectral_radius": float(max(max(r["baseline_rho_eff"], r["perturbed_rho_eff"]) for r in rows)),
        "decision": "GO" if go else "NO_GO",
        "claim_scope": contract["claim_scope"],
    }
    (outdir / contract["outputs"]["summary_json"]).write_text(json.dumps(summary, indent=2)+"\n", encoding="utf-8")
    gate = {
        "stage": "RECOMB-SA5R2 candidate gate from LGE-BRIDGE-01",
        "status": "PASS" if go else "FAIL",
        "decision": "GO_TO_MANUSCRIPT_REBUILD" if go else "NO_GO_LOCAL_EQUILIBRIUM_PRIMARY_CLAIM",
        "summary": summary,
    }
    (outdir / contract["outputs"]["gate_json"]).write_text(json.dumps(gate, indent=2)+"\n", encoding="utf-8")
    return summary


def preflight(contract: dict) -> None:
    # 1) Analytic zero-coupling identity.
    b = np.asarray([-0.4, 0.0, 0.3], dtype=float)
    W = np.zeros((3,3), dtype=float)
    sol = solve_equilibrium(b, W, 1e-13, 100)
    if not np.allclose(sol.u, np.tanh(b), atol=1e-12, rtol=0):
        raise RuntimeError("preflight zero-coupling identity failed")
    # 2) Frozen generator must be reproducible and globally contraction-safe.
    cfg = contract["world_generator"]
    r1 = np.random.Generator(np.random.PCG64DXSM(int(cfg["seed"])))
    r2 = np.random.Generator(np.random.PCG64DXSM(int(cfg["seed"])))
    w1 = generate_world(r1, 0, contract); w2 = generate_world(r2, 0, contract)
    for key in ("W","W_perm","b0"):
        if not np.array_equal(w1[key], w2[key]):
            raise RuntimeError(f"preflight reproducibility failed for {key}")
    if abs(float(np.linalg.norm(w1["W"], ord=2)) - float(cfg["operator_norm_target"])) > 1e-12:
        raise RuntimeError("preflight operator norm target failed")
    # 3) Toy matched solve must satisfy formal admissibility.
    b0 = w1["b0"]; b1=b0.copy(); b1[w1["source"]] += w1["delta"]
    s0=solve_equilibrium(b0,w1["W"],float(contract["solver"]["tolerance"]),int(contract["solver"]["max_iterations"]))
    s1=solve_equilibrium(b1,w1["W"],float(contract["solver"]["tolerance"]),int(contract["solver"]["max_iterations"]))
    if not admissible(s0,contract) or not admissible(s1,contract):
        raise RuntimeError("preflight toy equilibrium not admissible")
    print("PREFLIGHT_PASS")
    print(f"contract_sha256={EXPECTED_CONTRACT_SHA256}")
    print("formal_execution_authorized=false")
    print("formal results were not generated")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("mode", choices=["preflight","execute"])
    ap.add_argument("--outdir", default=str(HERE / "formal_results"))
    args = ap.parse_args()
    contract=load_contract()
    if args.mode == "preflight":
        preflight(contract); return 0
    if not AUTH_FILE.exists() or AUTH_FILE.read_text(encoding="utf-8").strip() != AUTH_TEXT:
        raise RuntimeError(
            "formal execution blocked: authorization token missing or incorrect; "
            "SA5R1 freezes protocol but does not authorize results"
        )
    summary=run_formal(contract,Path(args.outdir))
    print(json.dumps(summary,indent=2))
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as e:
        print(f"ERROR: {e}", file=sys.stderr)
        raise
