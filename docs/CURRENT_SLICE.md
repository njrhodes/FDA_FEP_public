# CURRENT_SLICE.md

## Goal

Audit historical Pmetrics 3.0.9 PTA regimen/weight-result association.

## Status

**COMPLETE.** Retrospective PTA mapping audit finished, operator sensitivity confirmed, and
ID-ordering risk interpretation documented. See `docs/PM_PTA_309_AUDIT.md` for full results.
No deterministic weight-row mis-association was found in the FDA II, EI, or CI workflow.
Mapping conclusion is robust to `>=` versus `>` at MIC 8 mg/L. The lexicographic ID-ordering
permutation identified in the non-FDA audit is not exposed by this FDA architecture (IDs 1–3
only). No further PTA mapping work planned.

## Scope

- FDA simulation/PTA path only.
- II, EI, CI.
- Existing 60/80/120 kg simulation templates.
- Compare PM_pta profile-level PDI vectors with direct PDI vectors from the actual simulation IDs.
- Establish whether returned PM_pta labels/results correspond to the intended weight template.
- Read-only scientific audit until a defect is demonstrated.

## Explicitly out of scope

- model fitting
- reviewer analyses
- figure redesign
- table changes
- manuscript prose
- package source archaeology
- fixing PM_pta
- changing simulation design
