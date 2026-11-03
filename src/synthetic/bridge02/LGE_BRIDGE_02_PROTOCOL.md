# RECOMB-SA5R3 — LGE-BRIDGE-02 Hard-Decoy Reliability Stress-Test Freeze

## Gate identity

- Stage: `RECOMB-SA5R3`
- Experiment: `LGE-BRIDGE-02`
- Status: `PROSPECTIVE_FROZEN_NOT_EXECUTED`
- Formal execution: **NOT AUTHORIZED IN THIS STAGE**
- Real-data access: **PROHIBITED**
- locked-3600 access: **PROHIBITED**
- Post-result threshold or world-definition changes: **PROHIBITED**

## Why BRIDGE-02 is required

LGE-BRIDGE-01 passed its prospectively frozen gate, but its deliberately high-confidence wrong target was a cross-module non-neighbor. That experiment established basic confidence–equilibrium decoupling, but the 128/128 primary win rate leaves a reviewer-facing concern that the decoy was structurally too easy.

BRIDGE-02 therefore leaves the equilibrium operator unchanged and changes only the difficulty of the candidate mechanisms. The primary decoy now has the **same source, same module, same direct-edge sign, identical direct-edge magnitude, matched motif degree, and identical absolute two-step path burden** as the true target. The discriminating feature is whether the surrounding signed paths are locally coherent or cancelling.

## Scientific hypothesis

A candidate mechanism can have strong one-step structural support and higher upstream confidence yet remain less locally self-consistent when its surrounding signed regulatory context is incoherent.

The primary hypothesis is

```text
R_eq(M_true) > R_eq(M_hard-direct)
```

when direct-edge evidence and lower-order path burden are matched.

## Equilibrium operator

The operator remains exactly the BRIDGE-01 reliability operator:

```text
Phi(u;b,W) = tanh(b + W u)
u* = Phi(u*;b,W)
J(u*) = diag(1-u*^2) W
```

Every matched and context-swap operator is independently rescaled to

```text
||W||_2 = 0.30,
```

so the map remains globally contraction-safe. The reliability quantity is unchanged:

```text
R_eq(M) = s_claim * (u1[target]-u0[target]) / |delta|.
```

`R_eq` is not a posterior probability.

## Frozen hard-decoy motif

Each world contains 32 nodes arranged into four modules of eight. One module is selected before outcomes exist and its nodes are randomly assigned to:

```text
source, true target, hard-direct target, alternative-path target, m1, m2, m3, m4.
```

Let `q` be a randomly drawn sign.

### True mechanism

- direct source→true edge: `q*w_direct`;
- two source→mediator→true paths;
- both two-step paths reinforce the claimed direction.

### Primary hard-direct decoy

- same source and same module;
- direct source→decoy edge is **exactly equal** to the true direct edge after scaling;
- same number of controlled incoming motif edges;
- same absolute two-step path burden as the true target;
- one two-step path reinforces and the other cancels;
- claimed direction is the same as the true mechanism;
- upstream confidence is higher (`0.97` vs `0.80`).

Thus a method relying only on direct-edge presence, direct-edge sign, direct-edge magnitude, module identity, or absolute path burden cannot distinguish the pair.

### Alternative-path boundary control

The alternative-path target is in the same module and has **no direct source edge**, but receives two reinforcing two-step paths. It is intentionally treated as a specificity-boundary control rather than a primary failure criterion because a strong indirect response can represent legitimate downstream compatibility even when it does not identify a direct mechanism.

### Context-swap control

A second operator is created prospectively by swapping only local sign coherence around the true and hard-direct targets while leaving their direct source edges unchanged. In the swap:

- true context becomes one reinforcing + one cancelling path;
- hard-direct context becomes two reinforcing paths.

A reversal of the ranking is therefore evidence that `R_eq` follows local regulatory context rather than target identity.

## Frozen synthetic generator

- RNG: NumPy `PCG64DXSM`
- world seed: `20261104`
- worlds: `128`
- nodes: `32`
- modules: `4 × 8`
- operator norm: `0.30`
- baseline field: Normal(0,0.08), clipped to ±0.25
- true and hard-direct target baseline fields: exactly matched
- perturbation magnitude: `|delta|=0.50`
- perturbation sign: positive for even world IDs, negative for odd world IDs
- bootstrap: 2,000 replicates, seed `20261105`

The selected motif module is constructed explicitly. Other modules and cross-module edges provide random background context. External incoming/outgoing context for the true and hard-direct targets is generated identically.

## Confidence stress test

These are deliberately adversarial synthetic confidence values, not calibrated biological probabilities:

```text
true              S = 0.80
hard direct decoy S = 0.97
alternative path  S = 0.95
wrong sign        S = 0.60
```

Confidence alone therefore ranks both hard decoys above the true mechanism.

## Primary endpoints

For every world:

```text
Delta_direct = R_eq(true) - R_eq(hard_direct).
```

Primary summaries:

1. hard-direct paired win fraction: `P(Delta_direct>0)`;
2. median `Delta_direct`;
3. context-swap reversal fraction: matched graph ranks `true>hard_direct` **and** context-swap graph ranks `hard_direct>true`.

Uncertainty:

- Wilson 95% CI for win fractions;
- 2,000-replicate world bootstrap percentile 95% CI for median `Delta_direct`.

## Secondary endpoints

- true vs alternative-path target;
- true vs wrong-sign control;
- context-swap response delta;
- fixed-point residuals and effective spectral radii.

Alternative-path interpretation is frozen before execution:

- `ALTPATH_SPECIFICITY_SUPPORTED`: win fraction ≥0.65 and Wilson lower >0.50;
- `ALTPATH_DOWNSTREAM_COMPATIBILITY_DOMINATES`: win fraction <0.35;
- otherwise `ALTPATH_BOUNDARY_AMBIGUOUS`.

This label cannot rescue or overturn the primary hard-direct gate.

## GO / NO-GO

### Hard preconditions

All 128 worlds must be present without exclusion. All role identities must be valid. Structural matching of the hard-direct pair must pass exactly. All matched and context-swap baseline/perturbed solves must satisfy residual ≤1e-10 and effective spectral radius ≤0.35. No nonfinite output is permitted. Frozen hashes must match.

### GO requires every criterion

```text
hard-direct win fraction >= 0.65
hard-direct Wilson 95% lower bound > 0.50
bootstrap 95% lower bound for median Delta_direct > 0
context-swap reversal fraction >= 0.65
context-swap Wilson 95% lower bound > 0.50
wrong-sign sanity win fraction >= 0.95
```

Failure of any criterion yields:

```text
NO_GO_HARD_DECOY_LOCAL_EQUILIBRIUM_PRIMARY_CLAIM
```

No threshold may change after authorization.

## Allowed claim if GO

> Under a prospectively frozen synthetic local-regulatory generator, local-equilibrium response distinguished a context-coherent direct mechanism from a deliberately higher-confidence direct-edge decoy matched on source, module, direct-edge sign and magnitude, motif degree, and absolute two-step path burden; a context-swap control further tested whether the ranking followed local signed regulatory context rather than target identity.

## Claims prohibited even if GO

BRIDGE-02 does not show that real genomic loci are biological steady states; does not produce a calibrated biological probability; does not establish real direct regulatory causality; does not map GWAS Z statistics into an equilibrium field; and does not authorize equilibrium-based filtering or reweighting of real-data PIPs.

## Execution governance

Formal execution requires an intentionally absent file:

```text
LGE_BRIDGE_02_EXECUTION_AUTHORIZATION.txt
```

with exact contents:

```text
AUTHORIZED_AFTER_RECOMB_SA5R3_FREEZE
```

Therefore SA5R3 can run preflight checks but cannot generate formal result files.
