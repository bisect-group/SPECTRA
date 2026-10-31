# SPECTRA_CC vs FASTERCC runtime comparison

Extends the existing `RuntimeComparisons/ForSPECTRA_CCResults.m` (SPECTRA_CC vs
FASTCC vs SWIFTCC) to FASTERCC, on a common Gurobi backend.

Reference:
Pires Pacheco, Gonzalez & Sauter (2026), *FASTERCC: Accelerating Flux Consistency
Testing and Context-Specific Reconstruction for Large-Scale Metabolic Network
Models*, bioRxiv `2026.03.19.712885`.

## What was run

| Experiment | Design | Arms | Status |
|---|---|---|---|
| 1a paper sizes | the 10 input sizes plotted in the SPECTRA paper (20–10,580 rxns), 3 draws × 3 timings | spectracc, fastercc_r0, fastcc, swiftcc | complete |
| 1b removal sweep | 0–10,000 removed in steps of 500, 3 × 3 | same four | complete |
| 2 expanded Recon3D | k = 5, 10, 20 (47,820–183,420 rxns), 4 removal levels, 3 × 3 | spectracc, fastercc_r0 | k≤20 complete; k=30, 50 running |
| 3 AGORA2 communities | 50 and 100 species under random media, 3 × 3 | spectracc, fastercc_r0 | complete at two media floors and two epsilons |
| 3t tolerance sweep | feasTol 1e-9/1e-8/1e-7/1e-6, 2 sizes, 2 media floors, 3 × 3 | spectracc | complete (144 tasks) |
| Harvey fixpoint | iterate each method on its own consistent set until stable | spectracc, spectracc_relaxed, fastercc_r0 | running |

Numbers are in `results_csv/`. One row per instance per arm: median runtime over
the timing repeats, min, max, LP count, consistent-set size.

**The per-instance `.mat` files are not committed.** They carry the full
consistent-set membership mask for every arm and run to ~39 GB. Every instance
is a pure function of its seed, so `bench_make_instance` / `bench_make_community`
regenerate any of them exactly; the CSVs carry everything the figures use.

## Why Gurobi and not CPLEX

The paper runs FASTERCC on CPLEX, and `fasterCC` hard-wires it: the LP layer
calls `cplexlp` directly (`LP7_rFASTCORMICS.m:32`, `LP3_rFASTCORMICS.m:25`),
bypassing COBRA. The IBM CPLEX MATLAB API ships only up to CPLEX 12.10, which is
why the paper pins MATLAB R2019b.

That route was tried first and abandoned. The CPLEX 12.10 MEX
(`cplexlink12100.mexa64`) **segfaults intermittently** under MATLAB R2023b, the
oldest MATLAB available on fir — three crashes in about six runs, while a
75-solve stress test passed clean. Non-deterministic failure is worse across a
long array run than an outright one.

So FASTERCC's LP layer is routed through `solveCobraLP` by
`fastercc_lp_backend/`, and every arm runs on Gurobi.

**Reporting consequence:** these are timings of *FASTERCC's algorithm on a common
Gurobi backend*, not of *FASTERCC as published*. The change cuts both ways — it
also removes the per-LP COBRA wrapper overhead that SPECTRA_CC paid and FASTERCC
did not — but it is a different claim and must be described as one.

## Arms

| Arm | What it is |
|---|---|
| `spectracc` | stock `spectraCC` on the instance as generated |
| `spectracc_relaxed` | stock `spectraCC` on the same instance with bounds relaxed so the zero flux vector is feasible (`lb>0→0`, `ub<0→0`) — exactly what `fasterCC` does to its own input at `fasterCC.m:40-46` |
| `fastercc_r0` | `fasterCC(..., revopt=0)` |
| `fastercc_r1` | `fasterCC(..., revopt=1)` |

`revopt` is FASTERCC's own documented 4th argument, not a variant of ours. The
paper never says which it timed. `revopt=1` runs only *after* the consistent set
`A` is final (`fasterCC.m:322` vs `:330`), so it cannot change `A` — it only
tightens bounds via a one-reaction-at-a-time LP loop. So `revopt=0` is the
like-for-like arm for a consistency check, and `revopt=1` is reported separately.

`spectracc_relaxed` exists because FASTERCC solves a *relaxed* problem. Without
that arm, any measured gap conflates algorithm with problem difficulty.

No independent third method runs in experiments 2 and 3, so arm disagreements
there cannot be adjudicated from the benchmark alone — they are stored for a
separate FVA pass rather than resolved by picking a winner.

`fastcc` and `swiftcc` are arms of experiment 1 only — it is the one regime
small enough (≤10,600 reactions) for them to finish. `swiftcc` takes `(S, rev)`
rather than a model, with `rev = model.lb < 0` to match
`ForSPECTRA_CCResults.m:5`, and is credited 1 LP as that script records it; its
reversible part is settled by a QR factorisation rather than by LPs, so its LP
count is not comparable with the others.

`spectracc_relaxed` is implemented in `harvey_fixpoint.m`, not in
`bench_variants.m` — it is a no-op on every model in experiments 1–3, all of
which already admit the zero flux vector. Harvey does not: 167 of its 81,094
reactions carry bounds excluding zero.

### A caveat that survives every arm

Both methods now reach Gurobi through `solveCobraLP`, which rebuilds and
validates the problem struct on every call, so the per-LP wrapper overhead is
charged to both. Report these as *end-to-end runtimes on a common backend*
rather than as a claim about LP strategy; the per-arm LP counts make the per-LP
cost visible.

## Test set

The paper's prose is self-inconsistent: "20 expanded models composed of 1–10
clusters" cannot reach its stated 183,920-reaction maximum, and 20 × 5 removal
levels = 100, not the 1140 it reports.

Reading the expansion code resolves it. `compute_expanded_input_model.m:177`
builds

```
rxns = [internal_rxns × k ; intercluster_transport × k ; trans ; Ex]
```

where `trans = Trans_<exmet>` and `Ex = Ex_<exmet>` are each `numel(exmets)` and
do **not** scale with `k`. So size is affine:

```
n(k) = k·(nBase − |exRxns|) + 2·|exmets|
```

Since `|exRxns| ≈ |exmets|`, `k=2` collapses to exactly `2·10600 = 21,200` — the
paper's stated minimum — and `k=19` gives ≈183.9k, its maximum. And
**19 × 6 × 10 = 1140**, exactly its stated total. The "1–10 clusters" is a typo
for 1–19.

Removal levels are the paper's five (`10, 50, 100, 200, 500`) plus `0`, the
fully consistent control. Both ends matter: the paper's speedup is 17.5–31.5× at
500 blocked but only 2.3–5.2× at 0.

## Stages

`bench_config` still carries the original three-stage ladder (`smoke`, `pilot`,
`full`), but the experiments above do not use it: every run passes
`bench_config('smoke')` and overrides the test set through `opts` (`ks`,
`levels`, `repOnly`, `reps`). The stage now selects little more than the
results directory, and `bench_config('full')` deliberately errors because it
was never sized.

## Reproducing it

The scheduler scripts used on our cluster are **not** part of this repository —
they encode one site's module names, partitions and absolute paths, which would
not transfer. Everything below is the portable form: each call is one MATLAB
invocation that can be run directly or wrapped in whatever scheduler you use.

```bash
# 0. one-time: patch fasterCC to report its LP count
cd $FASTERCC_SRC && patch -p0 < <harness>/fastercc_lpcount.patch

# Then export the four paths bench_config reads:
#   SPECTRA_SRC  FASTERCC_SRC  SCFASTCORMICS_SRC  COBRA_SRC
# plus AGORA2_DIR for experiment 3. bench_addpaths asserts each one resolves
# and fails loudly if a name is shadowed.

# 1. smoke test first — it checks the environment too
matlab -batch "cfg=bench_config('smoke'); bench_addpaths(cfg); bench_smoke"

# 2. experiment 1: Recon3D
matlab -batch "cfg=bench_config('smoke'); bench_addpaths(cfg); bench_prepare_base(cfg)"
matlab -batch "cfg=bench_config('smoke'); bench_addpaths(cfg); \
   bench_variants(cfg,'papersizes','spectra_first',struct('fva',false))"
matlab -batch "cfg=bench_config('smoke'); bench_addpaths(cfg); \
   bench_variants(cfg,'recon','spectra_first',struct('fva',false))"

# 3. experiment 2: expanded Recon3D — cache the expanded model per k first
matlab -batch "cfg=bench_config('smoke'); cfg.clusterCounts=[5 10 20]; \
   bench_addpaths(cfg); bench_generate_models(cfg)"
matlab -batch "cfg=bench_config('smoke'); bench_addpaths(cfg); \
   bench_variants(cfg,'expanded','spectra_first',struct('fva',false,'ks',20))"

# 4. experiment 3: AGORA2 communities
matlab -batch "cfg=bench_config('smoke'); bench_addpaths(cfg); bench_freeze_agora_pool(cfg)"
matlab -batch "cfg=bench_config('smoke'); bench_addpaths(cfg); \
   bench_variants(cfg,'community','spectra_first', \
      struct('fva',false,'ks',50,'epsilon',1e-8))"

# 5. Harvey, iterated to a fixpoint
matlab -batch "cfg=bench_config('smoke'); bench_addpaths(cfg); \
   harvey_fixpoint('spectracc',cfg)"
```

`opts` fields that matter when splitting work across jobs: `ks`, `levels`,
`repOnly` (pin one draw), `reps` (timings per process), `arms`, `epsilon`,
`feasTol`, and `tagSuffix` — which must differ between concurrent runs, since
two runs sharing a tag share a results directory and resume from each other's
checkpoint.

**Two operational warnings, if you run these at scale.**

Throttle concurrent jobs. Gurobi is commonly licensed through a shared token
server, and 84 of ours starting in the same second exhausted it: 79 died with
`Gurobi error 10009: Failed to connect to token server`. `bench_variants`
catches a failure per arm and continues, so such a task still exits 0 and
writes an **empty** result. Any collection step must count rows actually
produced rather than trust the exit status.

`bench_freeze_agora_pool` writes `agora2/pool.txt`, the frozen list of the 7,302
AGORA2 models community instances sample from. It refuses to overwrite: the
pool defines which community every seed builds, so changing it would silently
redefine every instance already computed.

## LP counting

Both runtime and LP count come from **one** run of each arm; nothing executes
twice.

- `spectraCC` returns its LP count as its 2nd output (stock behaviour).
- `fasterCC` returns it as a 3rd output, `nLP`, added by
  `fastercc_lpcount.patch` — 12 added lines, all counter increments, plus the
  signature. `fasterCC.m.orig` is kept beside it so the diff is auditable.
  `bench_check_env` verifies the patch is applied via `nargout('fasterCC')`.

## Files

| File | Role |
|---|---|
| `bench_config.m` | all settings; stage selects the test set |
| `bench_check_env.m` | MATLAB/CPLEX/COBRA/patch verification — run first |
| `bench_prepare_base.m` | Recon3D → flux-consistent core (~10,600 rxns) |
| `bench_generate_models.m` | expands per `k`, fits the size law, writes manifest |
| `bench_make_instance.m` | rebuilds one instance deterministically from its seed |
| `bench_seed.m` | `seed = f(k, nRemove, rep)` |
| `bench_run_one.m` | runs every arm once; runtime + LP count + consistent set |
| `bench_aggregate.m` | results → CSV, speedups, set agreement |
| `bench_smoke.m` | end-to-end pipeline check |
| `bench_task.m` | array-job entry point |
| `bench_addpaths.m` | path setup; asserts every dependency resolves and is not shadowed |
| `fastercc_lpcount.patch` | the FASTERCC LP-counter diff |
| `bench_variants.m` | the driver used for every experiment here: regimes, arms, per-arm checkpointing |
| `bench_make_community.m` | one AGORA2 community instance under a random medium |
| `bench_freeze_agora_pool.m` | freezes the 7,302-model sampling universe, once |
| `harvey_fixpoint.m` | iterates a checker on its own consistent set until stable |
| `prep_fastercc.m` | FASTERCC's structural preprocessing, applied to SPECTRA's input |
| `fastercc_lp_backend/` | routes FASTERCC's `cplexlp` calls through `solveCobraLP` |
| `probe_fastercc_mem.m` | RSS probe used to localise the FASTERCC OOM |
| `collect_all_results.m` | pools per-instance results into the CSVs |
| `results_csv/` | every result table; the inputs to the figures |

Instances are **not** stored. One top-end instance is ~184k reactions, so the
full set would run to terabytes; only the expanded model per `k` is cached and
each instance is replayed from its seed.

## Open items

- Scheduler scripts are not tracked here; sizing a run for a new site means
  re-deriving wall-clock and memory from the timings in `results_csv/`.
- Peak memory is not recorded. The paper does not report it either, despite
  scalability being its stated motivation.

## Open issue: fasterCC builds an empty LP on fully-blocked models

**Status: separate work item, not blocking the benchmarks.**

`fasterCC(model, 1e-4, 0, 0)` fails on any model whose consistent set is
empty:

```
solveCobraLP: gurobi reporting OPTIMAL but no solution
```

Reproduce with `diag_fastercc_crash.m`; `logs/diag_501158.out` holds the run.

### What is actually happening

The LP struct reaching the solver is empty, not ill-conditioned:

```
A: []   rhs: [0x1 double]   lb: [0x1 double]   ub: [0x1 double]   obj: [0x1 double]
```

Gurobi reports OPTIMAL - an empty problem trivially is - but returns no
solution vector, and `solveCobraLP` raises. So fasterCC calls the solver
without first checking that its candidate set is non-empty.

### Not a tolerance artefact

Ruled out, because we had introduced `feasTol = 1e-9` in the same round and it
was the obvious suspect. It fails identically at **1e-9, 1e-7 and 1e-6**.

### Not a general degeneracy problem

`spectraCC` on the same model returns `|cons| = 0` in 2 LPs, in 1.12 s.

### When it bites

Only at the very end of the removal sweep. Measured collapse of the
consistent set on Recon3D:

| removed | 0 | 4000 | 7000 | 9000 | 9500 | 10000 |
|---|---|---|---|---|---|---|
| surviving rxns | 10600 | 6600 | 3600 | 1600 | 1100 | 600 |
| \|cons\| | 10600 | 2953 | 502 | 139 | 80 | **0** |

So one level out of 21. The sweep keeps that level: `bench_variants` records
the arm's error and continues, and "fasterCC cannot process a fully blocked
model" is itself a result worth reporting.

### Likely fix, if pursued

Guard the solver call in `fasterCC.m` on a non-empty candidate set and return
an empty consistent set directly. Not attempted - any edit to fasterCC breaks
the "run it as published" rule this comparison depends on.
