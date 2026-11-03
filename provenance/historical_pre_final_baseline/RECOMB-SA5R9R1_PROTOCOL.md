# RECOMB-SA5R9R1 — Mechanistic-Closure Real-World Validation Freeze

## Scientific purpose
The validation unit is a **state transition**, not an edge label.

TopoMIEG-FM receives candidate regulatory relations from an upstream resource and asks whether those relations can be closed into a local regulatory mechanism whose perturbation-induced steady-state shift predicts an independently observed transcriptomic response.

## Primary system: K562

### Baseline state
Independent non-perturbational K562 baseline expression is taken from E-MTAB-11264.

### Candidate relation set
Use the human DoRothEA signed TF-target network, confidence levels A and B only.

TopoMIEG-FM does not infer these relations. They are treated as upstream candidate relations.

### Held-out perturbation state
Use the K562 genome-wide CRISPRi Perturb-seq study of Replogle et al. (Cell 2022; DOI 10.1016/j.cell.2022.05.013; BioProject PRJNA831566).

Perturbation-specific downstream expression responses remain hidden until the complete prediction manifest is frozen.

## Eligible mechanisms
A source perturbation enters the primary benchmark only if:
1. the perturbed gene is represented as a DoRothEA TF;
2. confidence-A/B outgoing relations exist;
3. source and targets are expressed in baseline K562;
4. at least 20 target genes are evaluable in the held-out expression matrix.

At least 30 source perturbations are required for the formal benchmark.

For source `s`, the mechanism node set is `s` plus its expressed confidence-A/B targets. The mechanism edge set is the induced confidence-A/B DoRothEA subgraph on those nodes. Thus target-target regulatory structure is retained when supported upstream.

## Steady-state representation
Let `x_obs` be a monotonic rank transform of baseline expression into [-0.95,0.95].

Let `W` be the signed candidate-relation operator, rescaled to spectral norm 0.30.

Anchor the observed state with:
`b = atanh(x_obs) - W x_obs`

so:
`x_obs = tanh(b + W x_obs)`

This baseline identity is a boundary condition, not the validation result.

## Perturbation prediction
Apply a fixed source inhibition:
`b_pert = b - 0.50 e_s`

Solve:
`x*_pert = tanh(b_pert + W x*_pert)`

Predict:
`Delta_pred = x*_pert - x_obs`

The perturbed source gene itself is excluded from primary scoring; only downstream state prediction is evaluated.

## Held-out response
For each eligible source, obtain the observed pseudobulk transcriptomic response relative to non-targeting controls from the K562 Perturb-seq resource.

Primary evaluation includes every evaluable gene in the frozen mechanism node set except the source. No outcome-derived DEG threshold is allowed for the primary metric.

## Primary metric
For source `s`:
`rho_s = Spearman(Delta_pred_s, Delta_obs_s)`

Primary aggregate:
`median_s(rho_s)`

Co-primary incremental quantity:
`median_s[rho_equilibrium_s - rho_one_step_s]`

where:
`Delta_one_step = -0.50 W e_s`

## Hard GO criteria
All must pass:
1. >=30 eligible source perturbations;
2. prediction manifest SHA256-frozen before outcome access;
3. 2,000-source-bootstrap 95% lower CI of median equilibrium rho > 0;
4. 2,000-source-bootstrap 95% lower CI of median equilibrium-minus-one-step rho > 0;
5. equilibrium rho exceeds the per-source rewired-network median in >=70% of eligible sources;
6. every primary fixed-point solve has residual <=1e-10 and effective spectral radius <=0.35.

A failed K562 primary gate cannot be rescued by the BCL11A case.

## Negative control
For every eligible source, generate 100 degree/sign-count-preserving target rewires using seed 20261107.

## Secondary biological case
BCL11A erythroid enhancer -> BCL11A -| HBG1/HBG2.

This case asks whether supplied relations close into the adult erythroid state and predict the experimentally observed transition toward lower BCL11A and higher HBG/HbF after enhancer editing.

It is explanatory and retrospective; it is not the blinded primary accuracy estimate.

## Leakage boundary
Perturbation-response vectors may not influence:
- source eligibility beyond predeclared identifier/coverage rules;
- edge inclusion or sign;
- operator weights;
- node sets;
- perturbation magnitude;
- solver settings;
- primary metrics.

No post-outcome tuning is permitted.

## Current status
This stage freezes the validation design only.
No held-out perturbation response has been scored.
