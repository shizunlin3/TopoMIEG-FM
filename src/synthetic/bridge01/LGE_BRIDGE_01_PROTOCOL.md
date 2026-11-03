# RECOMB-SA5R1 — LGE-BRIDGE-01 Prospective Protocol & Executable Freeze

## Gate identity

- Stage: `RECOMB-SA5R1`
- Experiment: `LGE-BRIDGE-01`
- Status at freeze: `PROSPECTIVE_FROZEN_NOT_EXECUTED`
- Formal execution in this stage: **NOT AUTHORIZED**
- Real-data access: **PROHIBITED**
- locked-3600 access: **PROHIBITED**
- Adaptive tuning after endpoint inspection: **PROHIBITED**

## Scientific question

Can a local-equilibrium response criterion distinguish a structurally supported synthetic genetic mechanism from a deliberately higher-confidence but structurally incorrect mechanism?

This experiment is designed to test the statement:

> Statistical confidence and local mechanistic self-consistency are distinct quantities.

The experiment does **not** test whether real genomic loci are biological steady states, and it does not modify or reweight any fine-mapping posterior.

## Why a new bridge operator is frozen here

Historical TopoMIEG-FM development contains validated synthetic equilibrium prototypes, including the P5C1B1X/Y/Z series. Those stages establish that synthetic fixed-point objects can be numerically and probabilistically qualified, but their frozen Gates explicitly do not define a real-GWAS-to-equilibrium mapping. The historical `full_equilibrium` label is therefore not revived.

For LGE-BRIDGE-01, local equilibrium is treated as an **orthogonal reliability operator**, not as a posterior distribution. The operator is defined prospectively and self-contained so its semantics are auditable.

## Local-equilibrium operator

For a local signed regulatory system with state `u in (-1,1)^p`, baseline field `b`, and directed interaction operator `W`, define

```text
Phi(u; b, W) = tanh(b + W u)
```

and the equilibrium state

```text
u* = Phi(u*; b, W).
```

The Jacobian at the fixed point is

```text
J(u*) = diag(1 - u*^2) W.
```

The synthetic generator rescales every `W` so that

```text
||W||_2 = 0.30.
```

Because `|d tanh(x)/dx| <= 1`, this gives a global contraction constant no larger than 0.30. Hence a unique fixed point exists for every frozen world under the designed operator. Formal admissibility nevertheless remains an explicit recorded endpoint rather than being assumed silently.

## Candidate mechanism object

A synthetic candidate mechanism is

```text
M = (source node, target node, claimed response direction, upstream confidence).
```

The mechanism does not enter a fine-mapping posterior. Instead, a fixed perturbation is applied to the source node and the equilibrium response at the target is measured.

Baseline equilibrium:

```text
u0* = Phi(u0*; b0, W)
```

Perturbed equilibrium:

```text
u1* = Phi(u1*; b0 + delta e_source, W).
```

The primary local-equilibrium reliability quantity is the signed response

```text
R_eq(M) = s_claim * (u1*[target] - u0*[target]) / |delta|.
```

Larger `R_eq` means that the equilibrium response agrees more strongly with the mechanism's claimed direction. For notation where a loss is needed,

```text
D_eq(M) = -R_eq(M),
```

so lower is more compatible. `D_eq` is not a posterior loss and has no biological-probability interpretation.

## Frozen synthetic world generator

- RNG: NumPy `PCG64DXSM`
- seed: `20261102`
- worlds: `128`
- nodes per world: `24`
- modules: `4`
- nodes per module: `6`
- directed within-module edge probability: `0.55`
- directed cross-module edge probability: `0.08`
- within-module absolute weight range: `[0.40, 1.00]`
- cross-module absolute weight range: `[0.10, 0.35]`
- within-module positive-edge probability: `0.80`
- cross-module positive-edge probability: `0.50`
- operator norm after rescaling: `0.30`
- baseline field: Normal(0, 0.08), clipped to `[-0.25, 0.25]`
- source perturbation magnitude: `|delta|=0.50`
- perturbation sign: positive in even-numbered worlds, negative in odd-numbered worlds

Each module receives a directed cycle before random edges are added, so every module has a nonempty regulatory backbone without outcome-dependent regeneration.

## Frozen candidate construction

### True mechanism `M+`

1. Select a source node that has at least one within-module outgoing edge and at least one external non-neighbor.
2. The true target is the largest-absolute-weight within-module outgoing edge from that source.
3. Claimed direction is `sign(delta * W[target, source])`.

### High-confidence wrong mechanism `M_HC-`

1. Use the same source.
2. Choose a target in a different module with no direct edge from the source.
3. Among eligible targets, minimize total in/out-degree mismatch to the true target; ties are broken by node index.
4. Use the same claimed direction as the true mechanism.
5. Assign higher upstream confidence than the true mechanism.

### Wrong-sign sanity control

Use the true source and true target but reverse the claimed response direction.

### Permuted-graph structural control

Apply a nonidentity simultaneous row/column permutation `P W P^T`, preserving operator spectrum while breaking the original source/target structural identity. Source and target labels remain fixed.

## Frozen confidence stress test

The confidence scores are intentionally adversarial and are not claimed to be calibrated real-world probabilities:

- true mechanism: `S=0.80`
- high-confidence wrong mechanism: `S=0.95`
- wrong-sign control: `S=0.60`

This construction forces a setting in which upstream confidence alone ranks the wrong mechanism above the true mechanism. The experiment asks whether equilibrium reliability supplies independent information under that stress condition.

## Solver and admissibility

- initialization: `u^(0)=tanh(b)`
- solver: Picard iteration
- solver tolerance: `1e-12`
- maximum iterations: `10000`
- formal residual threshold: `1e-10`
- formal effective spectral-radius threshold: `0.35`

A mechanism evaluation is admissible only when both the baseline and perturbed equilibria satisfy both thresholds.

## Primary endpoints

For each world define

```text
Delta_R = R_eq(M+) - R_eq(M_HC-).
```

Co-primary summaries:

1. **Paired win fraction**: fraction of worlds with `Delta_R > 0`.
2. **Median paired delta**: median of `Delta_R` across 128 worlds.

Uncertainty is frozen as:

- Wilson 95% CI for the paired win fraction.
- 2,000 world-level bootstrap replicates for the median paired delta.
- bootstrap RNG seed: `20261103`.
- percentile 95% CI.

## Secondary endpoints

- true mechanism versus wrong-sign control;
- true mechanism under matched graph versus true labels under a permuted graph;
- fixed-point residual distribution;
- effective spectral-radius distribution;
- descriptive confidence × equilibrium-reliability quadrant counts.

The permuted-graph result is secondary and is not allowed to replace a failed primary endpoint.

## GO / NO-GO rule

### Hard preconditions

All must hold:

- exactly 128 worlds generated;
- no post-hoc world exclusion;
- no duplicated true and wrong targets within a world;
- all matched baseline and perturbed equilibrium solves are admissible;
- no nonfinite outputs;
- frozen contract and runner hashes match the freeze manifest.

### GO requires all of the following

- paired win fraction `>= 0.70`;
- Wilson 95% lower bound for paired win fraction `> 0.50`;
- bootstrap 95% lower bound for median `Delta_R` `> 0`;
- true mechanism beats the wrong-sign sanity control in `>= 0.95` of worlds.

### NO-GO

Failure of any hard precondition or any GO criterion gives:

```text
NO_GO_LOCAL_EQUILIBRIUM_PRIMARY_CLAIM
```

No threshold may be changed after formal execution.

## Allowed manuscript claim if GO

> In controlled synthetic local systems, a graph-constrained equilibrium response distinguished structurally supported mechanisms from deliberately higher-confidence structural decoys, supporting local-equilibrium compatibility as an orthogonal reliability criterion rather than another posterior probability.

## Claims prohibited even if GO

- real genomic loci are at biological steady state;
- the method infers biological time constants;
- equilibrium reliability is a calibrated biological probability;
- equilibrium reliability improves or recalibrates real-data PIPs;
- real-GWAS summary statistics have been mapped into the equilibrium field;
- the experiment establishes disease mechanism or causal regulation.

## Execution governance

`SA5R1` freezes the protocol but **does not authorize formal execution**.

The formal runner requires a file named:

```text
LGE_BRIDGE_01_EXECUTION_AUTHORIZATION.txt
```

whose exact contents are:

```text
AUTHORIZED_AFTER_RECOMB_SA5R1_FREEZE
```

That file is intentionally absent from the frozen package. Therefore an accidental `execute` call fails closed.

Formal execution should occur only in the next authorized stage, after confirming the freeze manifest and repository location.
