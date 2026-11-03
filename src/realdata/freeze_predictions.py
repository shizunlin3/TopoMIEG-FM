#!/usr/bin/env python3
"""Freeze TopoMIEG-FM K562 predictions from the final ENCODE baseline and DoRothEA A/B graph.

This is the manuscript-facing release implementation. It deliberately has no access to
Perturb-seq outcomes or perturbation target lists.
"""
from __future__ import annotations
import argparse, gzip, hashlib, json, random
from pathlib import Path
import numpy as np
import pandas as pd

EXPECTED_DOROTHEA_SHA = "20d0a16f4dd93453bb5610e00c6c75e8f4acd047c186173077d2630e9753549e"
EXPECTED_BASELINE_SHA = "beb9fe6fb93f24b06b70cfaac418333c3a1585033d4416d2e89dab20c9e27c63"
OPERATOR_NORM = 0.30
PERTURBATION = -0.50
REWIRE_SEED = 20261107
N_REWIRES = 100
RESID_LIMIT = 1e-10
RHO_LIMIT = 0.35

def sha256(path: Path) -> str:
    h=hashlib.sha256()
    with path.open('rb') as f:
        for b in iter(lambda:f.read(1024*1024), b''): h.update(b)
    return h.hexdigest()

def rank_state(vals):
    s=pd.Series(vals); n=len(s); r=s.rank(method='average').to_numpy()
    return np.clip(2*((r-0.5)/n)-1, -0.95, 0.95)

def solve(W,b,x0,tol=1e-10,maxiter=10000):
    x=x0.copy()
    for nit in range(1,maxiter+1):
        y=np.tanh(b+W@x)
        if np.max(np.abs(y-x)) <= tol:
            x=y; break
        x=y
    else: raise RuntimeError('fixed-point solver did not converge')
    resid=float(np.max(np.abs(x-np.tanh(b+W@x))))
    eig=np.linalg.eigvals(np.diag(1-x*x)@W)
    rho=float(np.max(np.abs(eig))) if len(eig) else 0.0
    return x,resid,rho,nit

def signed_edges(sub):
    out=[]; seen=set()
    for row in sub.itertuples(index=False):
        u,v=str(row.source),str(row.target); sign=int(np.sign(row.mor))
        if sign==0 or u==v or (u,v) in seen: continue
        seen.add((u,v)); out.append((u,v,sign))
    return out

def rewire(edges,seed):
    rng=random.Random(int(seed)); by_sign={s:[(u,v) for u,v,z in edges if z==s] for s in (-1,1)}
    all_pairs={(u,v) for u,v,_ in edges}; total=0
    for sign in (-1,1):
        es=list(by_sign[sign])
        if len(es)<2: continue
        target=max(1,5*len(es)); attempts=max(100,50*len(es)); swaps=0
        for _ in range(attempts):
            if swaps>=target: break
            i,j=rng.sample(range(len(es)),2); a,b=es[i]; c,d=es[j]
            if a==c or b==d: continue
            p1,p2=(a,d),(c,b)
            if a==d or c==b or p1==p2: continue
            occupied=all_pairs-{(a,b),(c,d)}
            if p1 in occupied or p2 in occupied: continue
            all_pairs.discard((a,b)); all_pairs.discard((c,d)); all_pairs.add(p1); all_pairs.add(p2)
            es[i],es[j]=p1,p2; swaps+=1
        by_sign[sign]=es; total+=swaps
    return [(u,v,s) for s in (-1,1) for u,v in by_sign[s]], total

def operator(nodes,edges):
    idx={g:i for i,g in enumerate(nodes)}; W=np.zeros((len(nodes),len(nodes)),float)
    for u,v,s in edges:
        if u in idx and v in idx: W[idx[v],idx[u]]=float(s)
    raw=float(np.linalg.norm(W,2))
    if raw>0: W*=OPERATOR_NORM/raw
    return W,raw,idx

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--dorothea',type=Path,required=True)
    ap.add_argument('--baseline',type=Path,required=True)
    ap.add_argument('--outdir',type=Path,default=Path('predictions'))
    ap.add_argument('--allow-hash-mismatch',action='store_true',help='For exploratory use only; never for manuscript reproduction.')
    a=ap.parse_args(); a.outdir.mkdir(parents=True,exist_ok=True)
    observed={'dorothea':sha256(a.dorothea),'baseline':sha256(a.baseline)}
    if not a.allow_hash_mismatch:
        if observed['dorothea']!=EXPECTED_DOROTHEA_SHA: raise SystemExit('DoRothEA SHA256 mismatch')
        if observed['baseline']!=EXPECTED_BASELINE_SHA: raise SystemExit('ENCODE baseline SHA256 mismatch')
    net=pd.read_csv(a.dorothea,sep='\t')
    if 'tf' in net.columns and 'source' not in net.columns: net=net.rename(columns={'tf':'source'})
    req={'source','target','confidence','mor'}
    if not req.issubset(net.columns): raise SystemExit(f'DoRothEA columns missing: {sorted(req-set(net.columns))}')
    net=net[net.confidence.astype(str).isin(['A','B'])].copy(); net['mor']=np.sign(pd.to_numeric(net.mor,errors='coerce'))
    net=net[np.isfinite(net.mor)&(net.mor!=0)]
    base=pd.read_csv(a.baseline,sep='\t')
    if not {'gene_symbol','baseline_expression'}.issubset(base.columns): raise SystemExit('baseline requires gene_symbol, baseline_expression')
    base=base[['gene_symbol','baseline_expression']].dropna(); base['baseline_expression']=pd.to_numeric(base.baseline_expression,errors='raise')
    base=base[np.isfinite(base.baseline_expression)&(base.baseline_expression>=0)].groupby('gene_symbol',as_index=False).baseline_expression.mean()
    expr=dict(zip(base.gene_symbol.astype(str),base.baseline_expression.astype(float))); genes=set(expr)
    net=net[net.source.astype(str).isin(genes)&net.target.astype(str).isin(genes)].copy()
    counts=net.groupby('source').target.nunique(); sources=sorted(counts[counts>=20].index.astype(str))
    pred=[]; mech=[]; rw_path=a.outdir/'rewired_predictions.tsv.gz'
    with gzip.open(rw_path,'wt',encoding='utf-8',newline='') as rw:
        rw.write('source\trewire_id\ttarget\tdelta_equilibrium\tn_swaps\n')
        for source in sources:
            direct=sorted(set(net.loc[net.source.astype(str)==source,'target'].astype(str))); nodes=[source]+[g for g in direct if g!=source]
            sub=net[net.source.astype(str).isin(nodes)&net.target.astype(str).isin(nodes)].copy(); edges=signed_edges(sub)
            W,raw,idx=operator(nodes,edges); xobs=rank_state([expr[g] for g in nodes]); b=np.arctanh(xobs)-W@xobs
            bp=b.copy(); bp[idx[source]]+=PERTURBATION; xp,resid,rho,nit=solve(W,bp,xobs); deq=xp-xobs; done=PERTURBATION*W[:,idx[source]]
            if resid>RESID_LIMIT or rho>RHO_LIMIT: raise SystemExit(f'native numerical failure: {source}')
            for g in nodes:
                if g!=source: pred.append((source,g,float(deq[idx[g]]),float(done[idx[g]])))
            seed=(REWIRE_SEED+int(hashlib.sha256(source.encode()).hexdigest()[:8],16))%(2**32); swaps=[]
            for rid in range(N_REWIRES):
                redges,ns=rewire(edges,seed+rid); Wr,_,ridx=operator(nodes,redges); br=np.arctanh(xobs)-Wr@xobs; bpr=br.copy(); bpr[ridx[source]]+=PERTURBATION
                xpr,rr,er,_=solve(Wr,bpr,xobs)
                if rr>RESID_LIMIT or er>RHO_LIMIT: raise SystemExit(f'rewired numerical failure: {source}/{rid}')
                dpr=xpr-xobs; swaps.append(ns)
                for g in nodes:
                    if g!=source: rw.write(f'{source}\t{rid}\t{g}\t{float(dpr[ridx[g]]):.17g}\t{ns}\n')
            mech.append((source,len(nodes),len(edges),raw,resid,rho,nit,min(swaps),float(np.median(swaps))))
    pv=pd.DataFrame(pred,columns=['source','target','delta_equilibrium','delta_one_step'])
    mm=pd.DataFrame(mech,columns=['source','n_nodes','n_edges','raw_spectral_norm','residual_inf','effective_spectral_radius','iterations','min_rewire_swaps','median_rewire_swaps'])
    pv.to_csv(a.outdir/'prediction_vectors.tsv.gz',sep='\t',index=False,compression='gzip'); mm.to_csv(a.outdir/'mechanism_manifest.tsv',sep='\t',index=False)
    manifest={'blind_source_universe':True,'perturbation_outcomes_accessed':False,'n_sources':len(mm),'n_prediction_rows':len(pv),'rewires_per_source':N_REWIRES,'rewire_seed_base':REWIRE_SEED,'input_sha256':observed,'max_residual_inf':float(mm.residual_inf.max()),'max_effective_spectral_radius':float(mm.effective_spectral_radius.max())}
    (a.outdir/'prediction_manifest.release.json').write_text(json.dumps(manifest,indent=2)+'\n')
    print(json.dumps(manifest,indent=2))
if __name__=='__main__': main()
