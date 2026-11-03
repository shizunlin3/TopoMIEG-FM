#!/usr/bin/env python3
"""Materialize batch/guide-balanced Perturb-seq response vectors from the frozen K562 H5AD.

Aggregation: batch-specific perturbed-guide mean minus pooled explicit-control mean in the same batch;
equal-weight batches within guide; equal-weight guides within source.
"""
from __future__ import annotations
import argparse, hashlib, json
from pathlib import Path
import numpy as np, pandas as pd
from scipy import sparse
import anndata as ad

EXPECTED_H5AD_SHA='d3269d9c863d96555eba768b42ab330cc1b57493cd33b4ac8fe7fbab8278f104'

def sha256(p):
    h=hashlib.sha256()
    with open(p,'rb') as f:
        for b in iter(lambda:f.read(1024*1024),b''): h.update(b)
    return h.hexdigest()

def pick(obs, explicit, candidates, label):
    if explicit:
        if explicit not in obs.columns: raise SystemExit(f'{label} column {explicit!r} not present; available={list(obs.columns)}')
        return explicit
    hits=[x for x in candidates if x in obs.columns]
    if len(hits)!=1: raise SystemExit(f'Cannot unambiguously identify {label} column. Candidates found={hits}; use --{label}-key.')
    return hits[0]

def mean_rows(X, idx, cp10k=False):
    if len(idx)==0: raise ValueError('empty row set')
    Z=X[idx,:]
    if sparse.issparse(Z):
        Z=Z.tocsr().astype(np.float64)
        if cp10k:
            rs=np.asarray(Z.sum(axis=1)).ravel(); scale=np.divide(10000.0,rs,out=np.zeros_like(rs),where=rs>0)
            Z=sparse.diags(scale)@Z; Z.data=np.log1p(Z.data)
        return np.asarray(Z.mean(axis=0)).ravel()
    Z=np.asarray(Z,dtype=np.float64)
    if cp10k:
        rs=Z.sum(axis=1); Z=np.divide(Z*10000.0,rs[:,None],out=np.zeros_like(Z),where=rs[:,None]>0); Z=np.log1p(Z)
    return Z.mean(axis=0)

def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--h5ad',type=Path,required=True); ap.add_argument('--sources',type=Path,required=True); ap.add_argument('--outdir',type=Path,required=True)
    ap.add_argument('--representation',choices=['distributed','cp10k-log1p'],default='distributed')
    ap.add_argument('--gene-key'); ap.add_argument('--perturbation-key'); ap.add_argument('--nperts-key'); ap.add_argument('--batch-key'); ap.add_argument('--guide-key')
    ap.add_argument('--allow-hash-mismatch',action='store_true'); a=ap.parse_args(); a.outdir.mkdir(parents=True,exist_ok=True)
    hs=sha256(a.h5ad)
    if not a.allow_hash_mismatch and hs!=EXPECTED_H5AD_SHA: raise SystemExit(f'H5AD SHA256 mismatch: {hs}')
    A=ad.read_h5ad(a.h5ad,backed='r'); obs=A.obs
    gene=pick(obs,a.gene_key,['gene'],'gene'); pert=pick(obs,a.perturbation_key,['perturbation'],'perturbation'); nperts=pick(obs,a.nperts_key,['nperts'],'nperts')
    batch=pick(obs,a.batch_key,['batch','gem_group','gemgroup'],'batch'); guide=pick(obs,a.guide_key,['guide','sgID_AB','guide_id','guide_identity'],'guide')
    names=np.asarray(A.var_names.astype(str)); sources=[x.strip() for x in a.sources.read_text().splitlines() if x.strip()]; cp=a.representation=='cp10k-log1p'
    ctrl=((obs[gene].astype(str)=='non-targeting')&(obs[pert].astype(str)=='control')&(pd.to_numeric(obs[nperts],errors='coerce')==0)).to_numpy()
    if int(ctrl.sum())!=75328: raise SystemExit(f'Explicit-control count drift: {int(ctrl.sum())} != 75328')
    obs_gene=obs[gene].astype(str).to_numpy(); obs_pert=obs[pert].astype(str).to_numpy(); obs_np=pd.to_numeric(obs[nperts],errors='coerce').to_numpy(); obs_batch=obs[batch].astype(str).to_numpy(); obs_guide=obs[guide].astype(str).to_numpy()
    batches=sorted(set(obs_batch[ctrl])); control_means={}; control_counts=[]
    for b in batches:
        idx=np.flatnonzero(ctrl & (obs_batch==b)); control_means[b]=mean_rows(A.X,idx,cp); control_counts.append((b,len(idx)))
    source_rows=[]; guide_rows=[]; audits=[]
    for source in sources:
        sm=(obs_gene==source)&(obs_np==1); idxs=np.flatnonzero(sm)
        if len(idxs)==0: raise SystemExit(f'No cells for frozen source {source}')
        guides=sorted(set(obs_guide[idxs])); per_guide=[]; units=0; used_batches=set()
        for g in guides:
            gi=idxs[obs_guide[idxs]==g]; g_resp=[]
            for b in sorted(set(obs_batch[gi])):
                if b not in control_means: raise SystemExit(f'No explicit controls in batch {b}')
                bi=gi[obs_batch[gi]==b]; g_resp.append(mean_rows(A.X,bi,cp)-control_means[b]); units+=1; used_batches.add(b)
            gr=np.mean(np.vstack(g_resp),axis=0); per_guide.append(gr)
            guide_rows.extend((source,g,gn,float(v)) for gn,v in zip(names,gr))
        sr=np.mean(np.vstack(per_guide),axis=0); source_rows.extend((source,gn,float(v)) for gn,v in zip(names,sr))
        audits.append((source,len(idxs),len(guides),len(used_batches),units,len(names)))
    pd.DataFrame(source_rows,columns=['source','gene','delta_observed']).to_csv(a.outdir/'observed_response_vectors.tsv.gz',sep='\t',index=False,compression='gzip')
    pd.DataFrame(guide_rows,columns=['source','guide','gene','delta_observed']).to_csv(a.outdir/'guide_response_vectors.tsv.gz',sep='\t',index=False,compression='gzip')
    pd.DataFrame(audits,columns=['source','n_cells','n_guides','n_source_batches','n_guide_batch_units','n_measured_genes']).to_csv(a.outdir/'source_materialization_audit.tsv',sep='\t',index=False)
    pd.DataFrame(control_counts,columns=['batch','n_control_cells']).to_csv(a.outdir/'control_batch_counts.tsv',sep='\t',index=False)
    manifest={'h5ad_sha256':hs,'representation':a.representation,'columns':{'gene':gene,'perturbation':pert,'nperts':nperts,'batch':batch,'guide':guide},'n_controls':int(ctrl.sum()),'n_sources':len(sources),'n_genes':len(names),'aggregation':['batch-specific guide mean minus batch-specific pooled control mean','equal-weight batches within guide','equal-weight guides within source']}
    (a.outdir/'outcome_materialization_manifest.release.json').write_text(json.dumps(manifest,indent=2)+'\n')
    print(json.dumps(manifest,indent=2))
if __name__=='__main__': main()
