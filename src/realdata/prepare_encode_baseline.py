#!/usr/bin/env python3
"""Validate/register a materialized ENCODE K562 baseline table for the frozen release.

The exact manuscript baseline was built from ENCSR885DVH replicates ENCFF424CXV and ENCFF073NHK,
mapped to HGNC symbols. Because ENCODE quantification schemas can change, this release does not silently
reinterpret columns. Supply a two-column materialized table (`gene_symbol`, `baseline_expression`) and this
script performs the same final registration step and verifies the frozen hash.
"""
import argparse, hashlib
from pathlib import Path
import numpy as np, pandas as pd
EXPECTED='beb9fe6fb93f24b06b70cfaac418333c3a1585033d4416d2e89dab20c9e27c63'
def sha(p):
    h=hashlib.sha256();
    with open(p,'rb') as f:
        for b in iter(lambda:f.read(1<<20),b''): h.update(b)
    return h.hexdigest()
def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--input',type=Path,required=True); ap.add_argument('--output',type=Path,default=Path('baseline_expression.tsv')); ap.add_argument('--allow-hash-mismatch',action='store_true'); a=ap.parse_args()
    df=pd.read_csv(a.input,sep='\t'); req={'gene_symbol','baseline_expression'}
    if not req.issubset(df.columns): raise SystemExit(f'Expected columns {sorted(req)}')
    df=df[['gene_symbol','baseline_expression']].dropna(); df['gene_symbol']=df.gene_symbol.astype(str); df['baseline_expression']=pd.to_numeric(df.baseline_expression,errors='raise')
    df=df[np.isfinite(df.baseline_expression)&(df.baseline_expression>=0)].groupby('gene_symbol',as_index=False).baseline_expression.mean().sort_values('gene_symbol')
    df.to_csv(a.output,sep='\t',index=False); got=sha(a.output); print(f'genes={len(df)} sha256={got}')
    if not a.allow_hash_mismatch and got!=EXPECTED: raise SystemExit('Frozen baseline SHA256 mismatch. Use the exact materialization used by the manuscript or inspect provenance before proceeding.')
if __name__=='__main__': main()
