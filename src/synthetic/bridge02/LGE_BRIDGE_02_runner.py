#!/usr/bin/env python3
from __future__ import annotations
import argparse,csv,hashlib,json,math,sys
from dataclasses import dataclass
from pathlib import Path
from typing import Tuple
import numpy as np

HERE=Path(__file__).resolve().parent
CONTRACT=HERE/"LGE_BRIDGE_02_CONTRACT.json"
EXPECTED_CONTRACT_SHA256="76ff7888256a588d2b28b3c5ac0b33812c1287bfbaac3f47d24ef6bcf2a55cef"
AUTH_FILE=HERE/"LGE_BRIDGE_02_EXECUTION_AUTHORIZATION.txt"
AUTH_TEXT="AUTHORIZED_AFTER_RECOMB_SA5R3_FREEZE"

def sha256(path:Path)->str:
    h=hashlib.sha256()
    with path.open('rb') as f:
        for b in iter(lambda:f.read(1024*1024),b''): h.update(b)
    return h.hexdigest()

def load_contract():
    got=sha256(CONTRACT)
    if got!=EXPECTED_CONTRACT_SHA256: raise RuntimeError(f"contract SHA mismatch: {got} != {EXPECTED_CONTRACT_SHA256}")
    return json.loads(CONTRACT.read_text(encoding='utf-8'))

def phi(u,b,W): return np.tanh(b+W@u)

@dataclass(frozen=True)
class EqSolution:
    u:np.ndarray; residual_inf:float; rho_eff:float; iterations:int

def solve_equilibrium(b,W,tol,max_iter):
    u=np.tanh(b).astype(np.float64,copy=True)
    for it in range(1,max_iter+1):
        new=phi(u,b,W)
        if float(np.max(np.abs(new-u)))<=tol:
            u=new; break
        u=new
    else: raise RuntimeError('fixed-point solver did not converge')
    residual=float(np.max(np.abs(u-phi(u,b,W))))
    J=(1.0-u*u)[:,None]*W
    rho=float(np.max(np.abs(np.linalg.eigvals(J)))) if W.size else 0.0
    return EqSolution(u,residual,rho,it)

def scale_operator(W,target_norm):
    norm=float(np.linalg.norm(W,ord=2))
    if not np.isfinite(norm) or norm<=0: raise RuntimeError('non-positive operator norm')
    return W*(target_norm/norm)

def labels(p,nm): return np.repeat(np.arange(nm,dtype=int),p//nm)

def generate_world(rng,wid,contract):
    g=contract['world_generator']; p=int(g['p']); nm=int(g['n_modules']); lab=labels(p,nm)
    module=int(rng.integers(0,nm)); nodes=np.where(lab==module)[0]; roles=rng.permutation(nodes)
    source,true_t,direct_d,alt_d,m1,m2,m3,m4=[int(x) for x in roles]
    role_set={source,true_t,direct_d,alt_d,m1,m2,m3,m4}
    if len(role_set)!=8: raise RuntimeError('role identity collision')
    W=np.zeros((p,p),dtype=np.float64)

    # Background within-module structure only in the three non-motif modules.
    for m in range(nm):
        ns=np.where(lab==m)[0]
        if m==module: continue
        # directed cycle ensures nonempty background module
        for a,b in zip(ns,np.roll(ns,-1)):
            mag=float(rng.uniform(*g['background_within_abs_weight_range']))
            sign=1.0 if rng.random()<float(g['background_positive_probability']) else -1.0
            W[b,a]=sign*mag
        for s in ns:
            for t in ns:
                if s==t or W[t,s]!=0: continue
                if rng.random()<float(g['background_within_edge_probability_other_modules']):
                    mag=float(rng.uniform(*g['background_within_abs_weight_range']))
                    sign=1.0 if rng.random()<float(g['background_positive_probability']) else -1.0
                    W[t,s]=sign*mag

    # Cross-module background, excluding true/direct targets so their external context can be exactly matched.
    protected={true_t,direct_d}
    for s in range(p):
        for t in range(p):
            if s==t or lab[s]==lab[t] or s in protected or t in protected: continue
            if rng.random()<float(g['background_cross_edge_probability']):
                mag=float(rng.uniform(*g['background_cross_abs_weight_range']))
                sign=1.0 if rng.random()<float(g['background_positive_probability']) else -1.0
                W[t,s]=sign*mag

    # Exactly matched external context for true and hard-direct targets.
    for j in range(p):
        if lab[j]==module: continue
        if rng.random()<float(g['background_cross_edge_probability']):
            mag=float(rng.uniform(*g['background_cross_abs_weight_range']))
            sign=1.0 if rng.random()<float(g['background_positive_probability']) else -1.0
            W[true_t,j]=W[direct_d,j]=sign*mag
        if rng.random()<float(g['background_cross_edge_probability']):
            mag=float(rng.uniform(*g['background_cross_abs_weight_range']))
            sign=1.0 if rng.random()<float(g['background_positive_probability']) else -1.0
            W[j,true_t]=W[j,direct_d]=sign*mag

    q=1.0 if rng.random()<0.5 else -1.0
    wd=float(rng.uniform(*g['motif_direct_abs_weight_range']))
    a=float(rng.uniform(*g['motif_source_to_mediator_abs_weight_range']))
    b=float(rng.uniform(*g['motif_mediator_to_target_abs_weight_range']))

    # Equal one-step evidence for true and hard-direct decoy.
    W[true_t,source]=q*wd
    W[direct_d,source]=q*wd
    W[alt_d,source]=0.0

    # Equal absolute two-step burden. True is sign-coherent; direct decoy contains one cancelling path.
    for med in (m1,m2,m3,m4): W[med,source]=q*a
    W[true_t,m1]=+b; W[true_t,m2]=+b
    W[direct_d,m3]=+b; W[direct_d,m4]=-b
    # Alternative-path target has no direct edge but two reinforcing same-sign two-step paths.
    W[alt_d,m3]=+b; W[alt_d,m4]=+b

    # Context swap: only coherence around true/direct targets is exchanged; direct edges remain identical.
    Wswap=W.copy()
    Wswap[true_t,m2]=-b
    Wswap[direct_d,m4]=+b

    target_norm=float(g['operator_norm_target'])
    W=scale_operator(W,target_norm); Wswap=scale_operator(Wswap,target_norm)

    b0=np.clip(rng.normal(0.0,float(g['baseline_field_sd']),size=p),-float(g['baseline_field_clip']),float(g['baseline_field_clip']))
    shared=float(rng.normal(0.0,float(g['baseline_field_sd'])))
    shared=float(np.clip(shared,-float(g['baseline_field_clip']),float(g['baseline_field_clip'])))
    b0[true_t]=shared; b0[direct_d]=shared
    delta_abs=float(g['source_perturbation_abs_delta']); delta=delta_abs if wid%2==0 else -delta_abs
    s_claim=1.0 if delta*W[true_t,source]>0 else -1.0

    return dict(world_id=wid,W=W,Wswap=Wswap,b0=b0,module=module,source=source,true_target=true_t,
                direct_decoy=direct_d,alt_decoy=alt_d,m1=m1,m2=m2,m3=m3,m4=m4,
                s_claim=s_claim,delta=delta,raw_direct_abs=wd,raw_path_a=a,raw_path_b=b)

def response(base,pert,target,s_claim,delta):
    return float(s_claim*(pert.u[target]-base.u[target])/abs(delta))

def admissible(sol,c):
    s=c['solver']; return bool(sol.residual_inf<=float(s['formal_residual_threshold']) and sol.rho_eff<=float(s['formal_effective_spectral_radius_threshold']))

def wilson(k,n,z=1.959963984540054):
    ph=k/n; den=1+z*z/n; center=(ph+z*z/(2*n))/den; half=z*math.sqrt(ph*(1-ph)/n+z*z/(4*n*n))/den
    return max(0.0,center-half),min(1.0,center+half)

def boot_median(x,reps,seed):
    rng=np.random.Generator(np.random.PCG64DXSM(seed)); n=len(x); vals=np.empty(reps)
    for i in range(reps): vals[i]=np.median(x[rng.integers(0,n,size=n)])
    return float(np.quantile(vals,.025)),float(np.quantile(vals,.975))

def motif_contract(w):
    W=w['W']; s=w['source']; t=w['true_target']; d=w['direct_decoy']; a=w['alt_decoy']; m1,m2,m3,m4=w['m1'],w['m2'],w['m3'],w['m4']
    tol=1e-12
    if abs(W[t,s]-W[d,s])>tol: return False
    if abs(W[a,s])>tol: return False
    # equal absolute two-step burden in the controlled motif
    bt=abs(W[m1,s]*W[t,m1])+abs(W[m2,s]*W[t,m2])
    bd=abs(W[m3,s]*W[d,m3])+abs(W[m4,s]*W[d,m4])
    if abs(bt-bd)>tol: return False
    return True

def preflight(c):
    # analytic identity
    b=np.array([-0.4,0.0,0.3]); W=np.zeros((3,3)); sol=solve_equilibrium(b,W,1e-13,100)
    if not np.allclose(sol.u,np.tanh(b),atol=1e-12,rtol=0): raise RuntimeError('zero-coupling identity failed')
    g=c['world_generator']; r1=np.random.Generator(np.random.PCG64DXSM(int(g['seed']))); r2=np.random.Generator(np.random.PCG64DXSM(int(g['seed'])))
    w1=generate_world(r1,0,c); w2=generate_world(r2,0,c)
    for key in ('W','Wswap','b0'):
        if not np.array_equal(w1[key],w2[key]): raise RuntimeError(f'reproducibility failed: {key}')
    scalar_keys=('source','true_target','direct_decoy','alt_decoy','m1','m2','m3','m4','module')
    if any(w1[k]!=w2[k] for k in scalar_keys): raise RuntimeError('role reproducibility failed')
    if abs(float(np.linalg.norm(w1['W'],ord=2))-float(g['operator_norm_target']))>1e-12: raise RuntimeError('operator norm failed')
    if abs(float(np.linalg.norm(w1['Wswap'],ord=2))-float(g['operator_norm_target']))>1e-12: raise RuntimeError('swap operator norm failed')
    lab=labels(int(g['p']),int(g['n_modules'])); s=w1['source']; targets=(w1['true_target'],w1['direct_decoy'],w1['alt_decoy'])
    if not all(lab[x]==lab[s] for x in targets): raise RuntimeError('hard decoys not same-module')
    if not motif_contract(w1): raise RuntimeError('motif structural contract failed')
    # solver admissibility only; deliberately do not inspect primary response ranking in preflight.
    b0=w1['b0']; b1=b0.copy(); b1[s]+=w1['delta']
    for WW in (w1['W'],w1['Wswap']):
        x0=solve_equilibrium(b0,WW,float(c['solver']['tolerance']),int(c['solver']['max_iterations']))
        x1=solve_equilibrium(b1,WW,float(c['solver']['tolerance']),int(c['solver']['max_iterations']))
        if not admissible(x0,c) or not admissible(x1,c): raise RuntimeError('preflight equilibrium inadmissible')
    print('PREFLIGHT_PASS')
    print(f'contract_sha256={EXPECTED_CONTRACT_SHA256}')
    print('hard_decoy_structural_contract=PASS')
    print('formal_execution_authorized=false')
    print('primary endpoint rankings were not inspected')

def run_formal(c,outdir):
    g=c['world_generator']; sv=c['solver']; rng=np.random.Generator(np.random.PCG64DXSM(int(g['seed']))); rows=[]
    for wid in range(int(g['n_worlds'])):
        w=generate_world(rng,wid,c)
        if not motif_contract(w): raise RuntimeError(f'motif contract failed world {wid}')
        b0=w['b0']; b1=b0.copy(); b1[w['source']]+=w['delta']
        base=solve_equilibrium(b0,w['W'],float(sv['tolerance']),int(sv['max_iterations']))
        pert=solve_equilibrium(b1,w['W'],float(sv['tolerance']),int(sv['max_iterations']))
        base_s=solve_equilibrium(b0,w['Wswap'],float(sv['tolerance']),int(sv['max_iterations']))
        pert_s=solve_equilibrium(b1,w['Wswap'],float(sv['tolerance']),int(sv['max_iterations']))
        rt=response(base,pert,w['true_target'],w['s_claim'],w['delta'])
        rd=response(base,pert,w['direct_decoy'],w['s_claim'],w['delta'])
        ra=response(base,pert,w['alt_decoy'],w['s_claim'],w['delta'])
        rw=response(base,pert,w['true_target'],-w['s_claim'],w['delta'])
        rts=response(base_s,pert_s,w['true_target'],w['s_claim'],w['delta'])
        rds=response(base_s,pert_s,w['direct_decoy'],w['s_claim'],w['delta'])
        matched_ok=admissible(base,c) and admissible(pert,c); swap_ok=admissible(base_s,c) and admissible(pert_s,c)
        row={
          'world_id':wid,'source_module':w['module'],'source':w['source'],'true_target':w['true_target'],'hard_direct_target':w['direct_decoy'],'alternative_path_target':w['alt_decoy'],
          'delta':w['delta'],'claim_sign':w['s_claim'],'confidence_true':c['confidence_stress_test']['true_confidence'],'confidence_hard_direct':c['confidence_stress_test']['hard_direct_confidence'],'confidence_altpath':c['confidence_stress_test']['alternative_path_confidence'],
          'R_eq_true':rt,'R_eq_hard_direct':rd,'R_eq_alternative_path':ra,'R_eq_wrong_sign':rw,'Delta_direct':rt-rd,'true_beats_hard_direct':int(rt>rd),'true_beats_altpath':int(rt>ra),'true_beats_wrong_sign':int(rt>rw),
          'R_eq_true_context_swap':rts,'R_eq_hard_direct_context_swap':rds,'Delta_context_swap_direct_minus_true':rds-rts,'context_swap_reversal':int((rt>rd) and (rds>rts)),
          'direct_edge_true':float(w['W'][w['true_target'],w['source']]),'direct_edge_hard':float(w['W'][w['direct_decoy'],w['source']]),'alternative_direct_edge':float(w['W'][w['alt_decoy'],w['source']]),
          'operator_norm':float(np.linalg.norm(w['W'],ord=2)),'swap_operator_norm':float(np.linalg.norm(w['Wswap'],ord=2)),
          'baseline_residual':base.residual_inf,'perturbed_residual':pert.residual_inf,'baseline_rho':base.rho_eff,'perturbed_rho':pert.rho_eff,'swap_baseline_residual':base_s.residual_inf,'swap_perturbed_residual':pert_s.residual_inf,'swap_baseline_rho':base_s.rho_eff,'swap_perturbed_rho':pert_s.rho_eff,
          'matched_admissible':int(matched_ok),'swap_admissible':int(swap_ok)
        }
        if not all(np.isfinite(float(v)) for k,v in row.items()): raise RuntimeError(f'nonfinite world {wid}')
        if len({w['source'],w['true_target'],w['direct_decoy'],w['alt_decoy'],w['m1'],w['m2'],w['m3'],w['m4']})!=8: raise RuntimeError(f'role collision world {wid}')
        rows.append(row)
    outdir.mkdir(parents=True,exist_ok=True)
    fields=list(rows[0].keys()); tsv=outdir/c['outputs']['world_level_tsv']
    with tsv.open('w',encoding='utf-8',newline='') as f:
        wr=csv.DictWriter(f,fieldnames=fields,delimiter='\t',lineterminator='\n'); wr.writeheader(); wr.writerows(rows)
    n=len(rows); wins=sum(r['true_beats_hard_direct'] for r in rows); rev=sum(r['context_swap_reversal'] for r in rows); signwins=sum(r['true_beats_wrong_sign'] for r in rows); altwins=sum(r['true_beats_altpath'] for r in rows)
    wci=wilson(wins,n); rci=wilson(rev,n); aci=wilson(altwins,n)
    vals=np.asarray([r['Delta_direct'] for r in rows],float); bci=boot_median(vals,int(c['primary_endpoints']['uncertainty']['bootstrap_replicates']),int(c['primary_endpoints']['uncertainty']['bootstrap_seed']))
    hard=bool(n==int(g['n_worlds']) and all(r['matched_admissible']==1 and r['swap_admissible']==1 for r in rows))
    win=wins/n; reversal=rev/n; signf=signwins/n; altf=altwins/n
    go=bool(hard and win>=0.65 and wci[0]>0.50 and bci[0]>0 and reversal>=0.65 and rci[0]>0.50 and signf>=0.95)
    if altf>=0.65 and aci[0]>0.50: altlabel='ALTPATH_SPECIFICITY_SUPPORTED'
    elif altf<0.35: altlabel='ALTPATH_DOWNSTREAM_COMPATIBILITY_DOMINATES'
    else: altlabel='ALTPATH_BOUNDARY_AMBIGUOUS'
    summary={
      'stage':'LGE-BRIDGE-02 formal execution','contract_sha256':EXPECTED_CONTRACT_SHA256,'n_worlds':n,'matched_admissibility_fraction':sum(r['matched_admissible'] for r in rows)/n,'context_swap_admissibility_fraction':sum(r['swap_admissible'] for r in rows)/n,
      'hard_direct_win_fraction':win,'hard_direct_win_wilson95':list(wci),'median_Delta_direct':float(np.median(vals)),'median_Delta_direct_bootstrap95':list(bci),
      'context_swap_reversal_fraction':reversal,'context_swap_reversal_wilson95':list(rci),'wrong_sign_sanity_win_fraction':signf,
      'alternative_path_true_win_fraction':altf,'alternative_path_wilson95':list(aci),'alternative_path_interpretation':altlabel,
      'max_residual':float(max(max(r['baseline_residual'],r['perturbed_residual'],r['swap_baseline_residual'],r['swap_perturbed_residual']) for r in rows)),
      'max_effective_spectral_radius':float(max(max(r['baseline_rho'],r['perturbed_rho'],r['swap_baseline_rho'],r['swap_perturbed_rho']) for r in rows)),
      'decision':c['go_no_go']['go_decision'] if go else c['go_no_go']['no_go_decision'],'claim_scope':c['claim_scope']
    }
    (outdir/c['outputs']['summary_json']).write_text(json.dumps(summary,indent=2)+'\n',encoding='utf-8')
    gate={'stage':'RECOMB-SA5R4 candidate gate from LGE-BRIDGE-02','status':'PASS' if go else 'FAIL','decision':summary['decision'],'summary':summary}
    (outdir/c['outputs']['gate_json']).write_text(json.dumps(gate,indent=2)+'\n',encoding='utf-8')
    return summary

def main():
    ap=argparse.ArgumentParser(); ap.add_argument('mode',choices=['preflight','execute']); ap.add_argument('--outdir',default=str(HERE/'formal_results')); a=ap.parse_args(); c=load_contract()
    if a.mode=='preflight': preflight(c); return 0
    if not AUTH_FILE.exists() or AUTH_FILE.read_text(encoding='utf-8').strip()!=AUTH_TEXT:
        raise RuntimeError('formal execution blocked: authorization token missing or incorrect; SA5R3 freezes protocol but does not authorize results')
    print(json.dumps(run_formal(c,Path(a.outdir)),indent=2)); return 0

if __name__=='__main__':
    try: raise SystemExit(main())
    except Exception as e:
        print(f'ERROR: {e}',file=sys.stderr); raise
