#!/usr/bin/env python3
from __future__ import annotations
import argparse, gzip, hashlib, json, random, shutil, subprocess, sys, tarfile
from pathlib import Path
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent
DISC = ROOT / "assets" / "discovery"
VAULT = ROOT / "assets" / "vault"
PRED = ROOT / "predictions"
for d in (DISC, VAULT, PRED):
    d.mkdir(parents=True, exist_ok=True)

DORO_URL = "https://bioconductor.org/packages/3.23/data/experiment/src/contrib/dorothea_1.23.0.tar.gz"
FIG_API = "https://api.figshare.com/v2/articles/20029387"

def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for b in iter(lambda: f.read(1024*1024), b""):
            h.update(b)
    return h.hexdigest()

def download(url, path):
    import requests
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with requests.get(url, stream=True, timeout=120) as r:
        r.raise_for_status()
        with open(path, "wb") as f:
            for chunk in r.iter_content(1024*1024):
                if chunk:
                    f.write(chunk)

def preflight(_):
    if (ROOT / "UNSEAL_AUTHORIZATION.txt").exists():
        raise SystemExit("Unexpected unseal authorization present.")
    print("SA5R9R2_PREFLIGHT_PASS")
    print("outcome_unsealed=false")
    print("prediction_manifest_exists=" + str((PRED/"prediction_manifest.json").exists()).lower())

def acquire_dorothea(_):
    tgz = DISC / "dorothea_1.23.0.tar.gz"
    if not tgz.exists():
        download(DORO_URL, tgz)
    ext = DISC / "dorothea_pkg"
    if not ext.exists():
        ext.mkdir()
        with tarfile.open(tgz, "r:gz") as tf:
            tf.extractall(ext)

    rscript = shutil.which("Rscript")
    if not rscript:
        raise SystemExit("Rscript is required to export dorothea_hs.")

    rdas = list(ext.rglob("*dorothea_hs*.rda")) + list(ext.rglob("*dorothea_hs*.RData"))
    out = DISC / "dorothea_AB.tsv"

    if rdas:
        rcode = (
            'load("' + str(rdas[0]).replace("\\","/") + '"); '
            'x <- as.data.frame(dorothea_hs); '
            'x <- x[x$confidence %in% c("A","B"),]; '
            'write.table(x,file="' + str(out).replace("\\","/") +
            '",sep="\\t",row.names=FALSE,quote=FALSE)'
        )
        subprocess.run([rscript, "-e", rcode], check=True)
    else:
        rlib = DISC / "Rlib"
        rlib.mkdir(exist_ok=True)
        subprocess.run(["R", "CMD", "INSTALL", "-l", str(rlib), str(tgz)], check=True)
        rcode = (
            '.libPaths(c("' + str(rlib).replace("\\","/") + '",.libPaths())); '
            'library(dorothea); data(dorothea_hs,package="dorothea"); '
            'x <- as.data.frame(dorothea_hs); '
            'x <- x[x$confidence %in% c("A","B"),]; '
            'write.table(x,file="' + str(out).replace("\\","/") +
            '",sep="\\t",row.names=FALSE,quote=FALSE)'
        )
        subprocess.run([rscript, "-e", rcode], check=True)

    df = pd.read_csv(out, sep="\t")
    print("DOROTHEA_AB_FROZEN", len(df), sha256(out))

def inventory_baseline(_):
    import requests
    urls = [
        "https://www.ebi.ac.uk/biostudies/api/v1/studies/E-MTAB-11264/info",
        "https://www.ebi.ac.uk/biostudies/api/v1/studies/E-MTAB-11264"
    ]
    rec = {}
    for u in urls:
        r = requests.get(u, timeout=60)
        rec[u] = {"status": r.status_code}
        if r.ok:
            try:
                rec[u]["json"] = r.json()
            except Exception:
                rec[u]["text"] = r.text[:100000]
    p = DISC / "E-MTAB-11264_metadata.json"
    p.write_text(json.dumps(rec, indent=2))
    print("BASELINE_METADATA_FROZEN", sha256(p))

def register_baseline(args):
    src = Path(args.tsv).resolve()
    if "vault" in [x.lower() for x in src.parts]:
        raise SystemExit("Refusing baseline from vault.")
    df = pd.read_csv(src, sep="\t")
    if not {"gene_symbol","baseline_expression"}.issubset(df.columns):
        raise SystemExit("Expected columns gene_symbol and baseline_expression.")
    df = df[["gene_symbol","baseline_expression"]].dropna()
    df["gene_symbol"] = df["gene_symbol"].astype(str)
    df["baseline_expression"] = pd.to_numeric(df["baseline_expression"], errors="raise")
    df = df[np.isfinite(df["baseline_expression"]) & (df["baseline_expression"] >= 0)]
    df = df.groupby("gene_symbol", as_index=False)["baseline_expression"].mean()
    df = df.sort_values("gene_symbol")
    out = DISC / "baseline_expression.tsv"
    df.to_csv(out, sep="\t", index=False)
    prov = {
        "input_path": str(src),
        "input_sha256": sha256(src),
        "output_sha256": sha256(out),
        "n_genes": int(len(df)),
        "contract": "BASELINE_PREPROCESSING_CONTRACT.md"
    }
    (DISC/"baseline_expression.provenance.json").write_text(json.dumps(prov, indent=2))
    print("BASELINE_REGISTERED", len(df), sha256(out))

def vault_outcome(_):
    import requests
    r = requests.get(FIG_API, timeout=60)
    r.raise_for_status()
    meta = r.json()
    (VAULT/"figshare_20029387_metadata.json").write_text(json.dumps(meta, indent=2))

    preferred = ["K562_gwps_normalized_bulk_01.h5ad", "K562_gwps_raw_bulk_01.h5ad"]
    hit = None
    for name in preferred:
        for rec in meta.get("files", []):
            if rec.get("name") == name:
                hit = rec
                break
        if hit:
            break
    if hit is None:
        raise SystemExit("Expected K562 GWPS bulk file not found in Figshare metadata.")

    out = VAULT / hit["name"]
    if not out.exists():
        download(hit["download_url"], out)

    record = {
        "file": hit["name"],
        "file_id": hit.get("id"),
        "size": out.stat().st_size,
        "sha256": sha256(out),
        "supplied_md5": hit.get("supplied_md5"),
        "download_url": hit.get("download_url"),
        "PARSED": False
    }
    (VAULT/"OUTCOME_VAULT_MANIFEST.json").write_text(json.dumps(record, indent=2))
    print("OUTCOME_VAULTED_HASH_ONLY", record["file"], record["sha256"])

def rank_state(vals):
    s = pd.Series(vals)
    n = len(s)
    r = s.rank(method="average").to_numpy()
    x = 2*((r-0.5)/n)-1
    return np.clip(x, -0.95, 0.95)

def solve(W, b, x0, tol=1e-10, maxiter=10000):
    x = x0.copy()
    nit = 0
    for nit in range(1, maxiter+1):
        y = np.tanh(b + W @ x)
        if np.max(np.abs(y-x)) <= tol:
            x = y
            break
        x = y
    resid = float(np.max(np.abs(x - np.tanh(b + W @ x))))
    D = np.diag(1 - x*x)
    eig = np.linalg.eigvals(D @ W)
    rho = float(np.max(np.abs(eig))) if len(eig) else 0.0
    return x, resid, rho, nit

def _signed_edges_from_subgraph(sub):
    edges = []
    seen = set()
    for row in sub.itertuples(index=False):
        u, v = str(row.source), str(row.target)
        sign = int(np.sign(row.mor))
        key = (u, v)
        if sign == 0 or u == v:
            continue
        if key in seen:
            continue
        seen.add(key)
        edges.append((u, v, sign))
    return edges

def _degree_sign_preserving_rewire(edges, seed):
    """
    Directed double-edge swaps separately within positive and negative edge classes.
    This preserves signed out-degree and signed in-degree exactly.
    If a sparse graph admits no valid swap, the unchanged graph is retained as a
    conservative null replicate and n_swaps=0 is recorded.
    """
    rng = random.Random(int(seed))
    by_sign = {
        -1: [(u, v) for u, v, s in edges if s < 0],
         1: [(u, v) for u, v, s in edges if s > 0],
    }
    all_pairs = {(u, v) for u, v, _ in edges}
    total_swaps = 0

    for sign in (-1, 1):
        es = list(by_sign[sign])
        if len(es) < 2:
            continue
        target_swaps = max(1, 5 * len(es))
        max_attempts = max(100, 50 * len(es))
        swaps = 0
        for _ in range(max_attempts):
            if swaps >= target_swaps:
                break
            i, j = rng.sample(range(len(es)), 2)
            a, b = es[i]
            c, d = es[j]
            if a == c or b == d:
                continue
            p1, p2 = (a, d), (c, b)
            if a == d or c == b or p1 == p2:
                continue

            # Remove the two original pairs before checking proposed duplicates.
            occupied = all_pairs - {(a, b), (c, d)}
            if p1 in occupied or p2 in occupied:
                continue

            all_pairs.discard((a, b))
            all_pairs.discard((c, d))
            all_pairs.add(p1)
            all_pairs.add(p2)
            es[i], es[j] = p1, p2
            swaps += 1

        by_sign[sign] = es
        total_swaps += swaps

    rewired = [(u, v, sign) for sign in (-1, 1) for (u, v) in by_sign[sign]]
    return rewired, total_swaps

def _operator_from_edges(nodes, edges):
    idx = {g: i for i, g in enumerate(nodes)}
    W = np.zeros((len(nodes), len(nodes)), dtype=float)
    for u, v, sign in edges:
        if u in idx and v in idx:
            W[idx[v], idx[u]] = float(sign)
    raw_norm = float(np.linalg.norm(W, 2))
    if raw_norm > 0:
        W *= 0.30 / raw_norm
    return W, raw_norm, idx

def predict(_):
    doro = DISC / "dorothea_AB.tsv"
    base = DISC / "baseline_expression.tsv"
    if not doro.exists() or not base.exists():
        raise SystemExit("Discovery assets missing.")

    net = pd.read_csv(doro, sep="\t")
    if "tf" in net.columns and "source" not in net.columns:
        net = net.rename(columns={"tf": "source"})
    req = {"source", "target", "confidence", "mor"}
    if not req.issubset(net.columns):
        raise SystemExit("DoRothEA columns missing: " + str(req - set(net.columns)))

    net = net[net["confidence"].astype(str).isin(["A", "B"])].copy()
    net["mor"] = np.sign(pd.to_numeric(net["mor"], errors="coerce"))
    net = net[np.isfinite(net["mor"]) & (net["mor"] != 0)]

    bdf = pd.read_csv(base, sep="\t")
    expr = dict(zip(bdf["gene_symbol"].astype(str), bdf["baseline_expression"].astype(float)))
    genes = set(expr)
    net = net[
        net["source"].astype(str).isin(genes)
        & net["target"].astype(str).isin(genes)
    ].copy()

    # Strict blind source universe: never consult held-out perturbation target list.
    counts = net.groupby("source")["target"].nunique()
    sources = sorted(counts[counts >= 20].index.astype(str))

    pred_rows = []
    mech_rows = []
    rewire_path = PRED / "rewired_predictions.tsv.gz"

    with gzip.open(rewire_path, "wt", encoding="utf-8", newline="") as rw:
        rw.write("source\trewire_id\ttarget\tdelta_equilibrium\tn_swaps\n")

        for source in sources:
            direct_targets = sorted(
                set(net.loc[net["source"].astype(str) == source, "target"].astype(str))
            )
            nodes = [source] + [g for g in direct_targets if g != source]
            sub = net[
                net["source"].astype(str).isin(nodes)
                & net["target"].astype(str).isin(nodes)
            ].copy()

            edges = _signed_edges_from_subgraph(sub)
            W, raw_norm, idx = _operator_from_edges(nodes, edges)

            vals = [expr[g] for g in nodes]
            xobs = rank_state(vals)
            b = np.arctanh(xobs) - W @ xobs
            bp = b.copy()
            bp[idx[source]] -= 0.50

            xp, resid, rho, nit = solve(W, bp, xobs)
            delta_eq = xp - xobs
            delta_one = -0.50 * W[:, idx[source]]

            for gene in nodes:
                if gene == source:
                    continue
                i = idx[gene]
                pred_rows.append(
                    (source, gene, float(delta_eq[i]), float(delta_one[i]))
                )

            # Freeze 100 signed-degree-preserving rewired-null predictions/source.
            source_seed = (
                20261107
                + int(hashlib.sha256(source.encode("utf-8")).hexdigest()[:8], 16)
            ) % (2**32)
            swap_counts = []
            for rid in range(100):
                redges, n_swaps = _degree_sign_preserving_rewire(
                    edges, source_seed + rid
                )
                Wr, _, ridx = _operator_from_edges(nodes, redges)
                br = np.arctanh(xobs) - Wr @ xobs
                bpr = br.copy()
                bpr[ridx[source]] -= 0.50
                xpr, rresid, rrho, _ = solve(Wr, bpr, xobs)
                if rresid > 1e-10 or rrho > 0.35:
                    raise SystemExit(
                        f"Rewired numerical admissibility failure: "
                        f"{source} rewire={rid} residual={rresid} rho={rrho}"
                    )
                dpr = xpr - xobs
                swap_counts.append(n_swaps)
                for gene in nodes:
                    if gene == source:
                        continue
                    rw.write(
                        f"{source}\t{rid}\t{gene}\t"
                        f"{float(dpr[ridx[gene]]):.17g}\t{n_swaps}\n"
                    )

            mech_rows.append(
                (
                    source,
                    len(nodes),
                    len(edges),
                    raw_norm,
                    resid,
                    rho,
                    nit,
                    min(swap_counts) if swap_counts else 0,
                    float(np.median(swap_counts)) if swap_counts else 0.0,
                )
            )

    pv = pd.DataFrame(
        pred_rows,
        columns=["source", "target", "delta_equilibrium", "delta_one_step"],
    )
    mm = pd.DataFrame(
        mech_rows,
        columns=[
            "source",
            "n_nodes",
            "n_edges",
            "raw_spectral_norm",
            "residual_inf",
            "effective_spectral_radius",
            "iterations",
            "min_rewire_swaps",
            "median_rewire_swaps",
        ],
    )

    if len(mm) == 0:
        raise SystemExit("No blind mechanisms passed the >=20 target rule.")
    if not (
        (mm["residual_inf"] <= 1e-10)
        & (mm["effective_spectral_radius"] <= 0.35)
    ).all():
        raise SystemExit("Numerical admissibility failure before outcome unseal.")

    pv.to_csv(
        PRED / "prediction_vectors.tsv.gz",
        sep="\t",
        index=False,
        compression="gzip",
    )
    mm.to_csv(PRED / "mechanism_manifest.tsv", sep="\t", index=False)

    manifest = {
        "blind_source_universe": True,
        "replogle_target_list_consulted": False,
        "outcome_parsed": False,
        "n_sources": int(len(mm)),
        "n_prediction_rows": int(len(pv)),
        "rewires_per_source": 100,
        "rewire_seed_base": 20261107,
        "dorothea_sha256": sha256(doro),
        "baseline_sha256": sha256(base),
        "prediction_vectors_sha256": sha256(PRED / "prediction_vectors.tsv.gz"),
        "rewired_predictions_sha256": sha256(rewire_path),
        "mechanism_manifest_sha256": sha256(PRED / "mechanism_manifest.tsv"),
    }
    (PRED / "prediction_manifest.json").write_text(json.dumps(manifest, indent=2))
    print("PREDICTION_MANIFEST_FROZEN")
    print(json.dumps(manifest, sort_keys=True))

def freeze(_):
    if not (PRED/"prediction_manifest.json").exists():
        raise SystemExit("prediction_manifest.json missing.")
    records = []
    for d in [DISC, VAULT, PRED]:
        for p in sorted(d.rglob("*")):
            if p.is_file():
                records.append({
                    "path": str(p.relative_to(ROOT)),
                    "size": p.stat().st_size,
                    "sha256": sha256(p)
                })
    out = {"outcome_unsealed": False, "files": records}
    p = ROOT/"SHA256SUMS.json"
    p.write_text(json.dumps(out, indent=2))
    print("SA5R9R2_FREEZE_COMPLETE", sha256(p))

def main():
    ap = argparse.ArgumentParser()
    sp = ap.add_subparsers(dest="cmd", required=True)
    commands = {
        "preflight": preflight,
        "acquire-dorothea": acquire_dorothea,
        "inventory-baseline": inventory_baseline,
        "vault-outcome": vault_outcome,
        "predict": predict,
        "freeze": freeze
    }
    for name, fn in commands.items():
        p = sp.add_parser(name)
        p.set_defaults(fn=fn)
    p = sp.add_parser("register-baseline")
    p.add_argument("--tsv", required=True)
    p.set_defaults(fn=register_baseline)
    args = ap.parse_args()
    args.fn(args)

if __name__ == "__main__":
    main()
