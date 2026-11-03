#!/usr/bin/env python3
from pathlib import Path
import hashlib, json, pandas as pd, sys
ROOT=Path(__file__).resolve().parents[1]
def sha(p):
    h=hashlib.sha256();
    with p.open('rb') as f:
        for b in iter(lambda:f.read(1<<20),b''): h.update(b)
    return h.hexdigest()
checks={
 'src/synthetic/bridge01/LGE_BRIDGE_01_runner.py':'1af2fd91924c974722fddc2dc486cb08406b4667b77cfe625f8d95bfdea7b426',
 'src/synthetic/bridge01/LGE_BRIDGE_01_CONTRACT.json':'f35834415be43ed26c0300c3778455ca0e27f63dc0b0333165ee4a2c64c3b7c5',
 'src/synthetic/bridge02/LGE_BRIDGE_02_runner.py':'c1f435dafc9b00e875919dc325988af8436f00ddc20ceecb700b5e5a6c65c07a',
 'src/synthetic/bridge02/LGE_BRIDGE_02_CONTRACT.json':'76ff7888256a588d2b28b3c5ac0b33812c1287bfbaac3f47d24ef6bcf2a55cef',
}
for rel,exp in checks.items():
    got=sha(ROOT/rel)
    if got!=exp: raise SystemExit(f'HASH FAIL {rel}: {got}')
print('synthetic frozen hashes: PASS')
t=pd.read_csv(ROOT/'results/tables/Table_S7.tsv',sep='\t')
if list(t['Sources'].astype(int))!=[39,39]: raise SystemExit('Table S7 sources drift')
vals=t['Median eq. rho'].astype(float).round(4).tolist()
if vals!=[-0.0455,-0.0603]: raise SystemExit(f'Table S7 median drift: {vals}')
print('final supplement headline values: PASS')
conf=json.loads((ROOT/'config/release.json').read_text())
if conf['prediction_freeze']['n_sources']!=39 or conf['prediction_freeze']['n_prediction_rows']!=3383: raise SystemExit('release config drift')
print('release config: PASS')
print('RELEASE_VERIFY_PASS')
