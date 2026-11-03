#!/usr/bin/env python3
"""Score frozen equilibrium/one-step/rewired predictions against held-out K562 responses."""
from __future__ import annotations
import argparse, json
from pathlib import Path
import numpy as np, pandas as pd
from scipy.stats import spearmanr

def rho(a,b):
    a=np.asarray(a,float); b=np.asarray(b,float)
    if len(a)<3 or np.unique(a).size<2 or np.unique(b).size<2: return np.nan
    return float(spearmanr(a,b).statistic)

def boot_median(x,n=2000,seed=20261108):
    x=np.asarray(x,float); x=x[np.isfinite(x)]; rng=np.random.default_rng(seed)
    sims=np.median(x[rng.integers(0,len(x),size=(n,len(x)))],axis=1)
    return [float(np.quantile(sims,.025)),float(np.quantile(sims,.975))]

def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--predictions',type=Path,required=True); ap.add_argument('--rewired',type=Path,required=True); ap.add_argument('--observed',type=Path,required=True); ap.add_argument('--outdir',type=Path,required=True); a=ap.parse_args(); a.outdir.mkdir(parents=True,exist_ok=True)
    p=pd.read_csv(a.predictions,sep='\t'); r=pd.read_csv(a.rewired,sep='\t'); o=pd.read_csv(a.observed,sep='\t')
    if 'target' not in o.columns and 'gene' in o.columns: o=o.rename(columns={'gene':'target'})
    rows=[]; rwrows=[]
    for s in sorted(p.source.unique()):
        ps=p[p.source==s].copy(); os=o[o.source==s][['target','delta_observed']]; m=ps.merge(os,on='target',how='inner')
        req=r[r.source==s]; eq=rho(m.delta_equilibrium,m.delta_observed); one=rho(m.delta_one_step,m.delta_observed); null=[]; exact=0; predcorr=[]; native=dict(zip(m.target,m.delta_equilibrium))
        for rid,g in req.groupby('rewire_id'):
            z=g[['target','delta_equilibrium']].merge(os,on='target',how='inner'); rr=rho(z.delta_equilibrium,z.delta_observed); null.append(rr)
            common=[x for x in z.target if x in native]; rv=np.array([float(z.loc[z.target==x,'delta_equilibrium'].iloc[0]) for x in common]); nv=np.array([native[x] for x in common])
            if len(common) and np.allclose(rv,nv,atol=1e-12,rtol=0): exact+=1
            predcorr.append(rho(rv,nv))
        null=np.asarray(null,float); med=float(np.nanmedian(null)); tied=float(np.mean(np.isclose(null,eq,atol=1e-12,rtol=0)))
        rows.append((s,len(m),eq,one,eq-one if np.isfinite(one) else np.nan,med,bool(eq>med) if np.isfinite(eq) else False))
        rwrows.append((s,eq,med,float(np.nanstd(null)),tied,exact/len(null),float(np.nanmedian(req.groupby('rewire_id').n_swaps.first())),float(np.nanmedian(predcorr))))
    df=pd.DataFrame(rows,columns=['source','n_genes','rho_equilibrium','rho_one_step','rho_difference','median_rewired','beats_rewired_median'])
    rw=pd.DataFrame(rwrows,columns=['source','native_rho','null_median','null_sd','tied_fraction','exact_vector_fraction','median_swaps','median_prediction_rho_native_rewired'])
    df.to_csv(a.outdir/'source_scores.tsv',sep='\t',index=False); rw.to_csv(a.outdir/'rewiring_diagnostics.tsv',sep='\t',index=False)
    paired=df[np.isfinite(df.rho_one_step)]; summary={'n_sources':len(df),'median_equilibrium_rho':float(np.nanmedian(df.rho_equilibrium)),'equilibrium_bootstrap95':boot_median(df.rho_equilibrium),'n_paired_sources':len(paired),'median_equilibrium_minus_one_step':float(np.nanmedian(paired.rho_difference)),'paired_bootstrap95':boot_median(paired.rho_difference),'fraction_native_gt_rewired_median':float(df.beats_rewired_median.mean()),'median_tied_rewire_fraction':float(np.nanmedian(rw.tied_fraction)),'n_all_rewires_tied':int((rw.tied_fraction==1).sum()),'n_all_rewires_exact_prediction':int((rw.exact_vector_fraction==1).sum())}
    (a.outdir/'aggregate_summary.json').write_text(json.dumps(summary,indent=2)+'\n'); print(json.dumps(summary,indent=2))
if __name__=='__main__': main()
