# PM_PTA_309_AUDIT.md

## Purpose

Retrospective read-only audit of the historical FDA Pmetrics 3.0.9 PTA workflow. The audit
tests whether `PM_pta` assigns profile-level PDI/success results to the correct regimen-weight
row for each of the three FDA maintenance regimens (II, EI, CI). The audit was conducted after
a separate non-FDA workflow was found to potentially mis-associate rows when multiple simulation
objects are combined before passing to `PM_pta`.

## FDA runtime

- R 4.4.2
- Pmetrics 3.0.9
- Historical development model: Run 38 final object loaded directly from stored
  `Pmetrics/Runs/38/outputs/PMout.Rdata` via `base::load()` into an isolated environment.
  `core$final` (NPAG, PM_final, R6; 33 support points × 10 parameters) used as `poppar`.
- `PM_result$new()` / `PM_load()` were not used; they fail internally for the stored runs due
  to a Pmetrics 3.0.9 incompatibility in `alt_mod_lib_names()$alt`.

## Existing simulation architecture

Each FDA maintenance regimen is simulated independently:

- `PM_sim$new(... split=TRUE ...)` called separately for II, EI, and CI.
- Each simulation object contains three weight templates (60, 80, 120 kg) as separate subject IDs.
- `PM_pta$new(simdata=sim, simlabels=lbls, ...)` is called on each single-regimen simulation
  object in isolation; it is never called on a combined multi-regimen object.

This architecture means each `PM_pta` object contains exactly 3 rows (one per weight class)
and the `simlabels` vector is constructed from the same template used for simulation.

## Audit method

For each regimen (II, EI, CI):

1. Reproduce `PM_sim` and `PM_pta` under identical settings as the historical workflow.
2. Extract `sim$data$obs` (profile-level ELF concentrations, outeq = 5) from the same
   simulation object.
3. Compute a direct PDI vector independently for each template ID (weight class):
   `mean(free_fraction * C_ELF >= MIC)` over all simulated time points, per profile.
4. For each `PM_pta` row, compare its internal 1000-element `pdi` list-column against all
   three direct weight-specific PDI vectors using MAD and Pearson correlation.
5. Identify the best-matching direct weight by minimum MAD and confirm it equals the
   intended weight.

Diagnostic settings:

| Parameter | Value |
|-----------|-------|
| ELF compartment | outeq = 5 |
| MIC | 8 mg/L |
| Success threshold | 68% fT>MIC |
| Free fraction | 0.8 |
| nsim | 1000 |
| seed | 12345 |
| PM_pta start | 0 |
| PM_pta end | Inf |

## Results

### II regimen (2 g q12 intermittent infusion + loading dose)

| PM row | PM label | Intended wt | Best-match wt | Intended MAD | Best MAD | Correlation | Exact? | Next-best MAD |
|--------|----------|-------------|---------------|-------------|----------|-------------|--------|--------------|
| 1 | 2 g q12 II + LD; 60kg | 60 kg | 60 kg | 0.001822 | 0.001822 | 0.999973 | No | 0.068569 |
| 2 | 2 g q12 II + LD; 80kg | 80 kg | 80 kg | 0.002369 | 0.002369 | 0.999970 | No | 0.068569 |
| 3 | 2 g q12 II + LD; 120kg | 120 kg | 120 kg | 0.002976 | 0.002976 | 0.999975 | No | 0.123166 |

`prop_success == mean(success)`: TRUE for all rows.

Intended match clearly better than all alternatives (30–100× lower MAD).

### EI regimen (2 g q12 extended infusion over 4 h + loading dose)

| PM row | PM label | Intended wt | Best-match wt | Intended MAD | Best MAD | Correlation | Exact? | Next-best MAD | Sep ratio |
|--------|----------|-------------|---------------|-------------|----------|-------------|--------|--------------|-----------|
| 1 | 2 g q12 EI + LD; 60kg | 60 kg | 60 kg | 0.001328 | 0.001328 | 0.999977 | No | 0.060434 | 46× |
| 2 | 2 g q12 EI + LD; 80kg | 80 kg | 80 kg | 0.001892 | 0.001892 | 0.999974 | No | 0.060434 | 48× |
| 3 | 2 g q12 EI + LD; 120kg | 120 kg | 120 kg | 0.002352 | 0.002352 | 0.999977 | No | 0.118908 | 51× |

`prop_success == mean(success)`: TRUE for all rows.

Intended match 46–51× better than nearest alternative by MAD.

### CI regimen (4 g q24 continuous infusion + loading dose)

| PM row | PM label | Intended wt | Best-match wt | Intended MAD | Best MAD | Correlation | Exact? | Next-best MAD | Sep ratio |
|--------|----------|-------------|---------------|-------------|----------|-------------|--------|--------------|-----------|
| 1 | 4 g q24 CI + LD; 60kg | 60 kg | 60 kg | 0.000029 | 0.000029 | 1.000000 | No | 0.047624 | 1648× |
| 2 | 4 g q24 CI + LD; 80kg | 80 kg | 80 kg | 0.000010 | 0.000010 | 1.000000 | No | 0.047642 | 4995× |
| 3 | 4 g q24 CI + LD; 120kg | 120 kg | 120 kg | 0.000032 | 0.000032 | 1.000000 | No | 0.111480 | 3481× |

`prop_success == mean(success)`: TRUE for all rows.

Intended match 1648–4995× better than nearest alternative by MAD. Correlation = 1.000 to
three decimal places for all rows.

## Operator sensitivity

A confirmatory sensitivity check repeated the direct PDI fingerprint using strict `> MIC`
in place of `>= MIC`. All nine rows retained the same intended 60/80/120 kg best-match
mapping. MAD and correlation values were numerically identical to those from the primary
`>=` audit. The mapping conclusion is therefore robust to the choice of `>=` versus `>` at
the audited MIC of 8 mg/L.

## Interpretation

The targeted PDI-vector fingerprint audit found **no deterministic weight/result
mis-association** in the historical FDA Pmetrics 3.0.9 workflow for II, EI, or CI.

For all nine rows, the PM_pta PDI vector overwhelmingly best matched the direct PDI vector
for its intended weight, with correlations ≥0.99997 and intended-weight MAD far smaller than
alternative-weight MAD. Separation ratios ranged from approximately 30× (II) to approximately
5000× (CI) over the nearest alternative.

The FDA workflow differs architecturally from the non-FDA failing workflow: II, EI, and CI are
evaluated in separate `PM_sim` / `PM_pta` objects, each containing three weight-specific IDs.
No deterministic weight/result mis-association was observed under this architecture at the
audited diagnostic condition.

The FDA workflow does not expose the specific multi-digit lexicographic ID-ordering permutation
identified in the separate non-FDA audit because each `PM_pta` object contains only IDs 1–3,
for which numeric and lexicographic order are identical. This is consistent with the observed
correct mapping for all nine audited FDA rows.

The small nonzero MAD values observed in II and EI rows (~0.001–0.003) did not affect mapping
identification and their source was not investigated. No claim is made about their mechanism.

This audit does not assert that Pmetrics 3.0.9 PM_pta is universally valid across all
architectures. It does not rule out the possibility that row mis-association could occur under
other ID structures or if multiple simulation objects were combined before passing to `PM_pta`.
It applies only to the specific FDA single-regimen-per-PM_pta architecture with IDs 1–3
described above.

## Required remediation

No remediation is required for the tested historical FDA PTA mapping workflow. No deterministic weight/result mis-association was identified under the audited architecture and diagnostic condition.
