# CURRENT_SLICE.md

## Goal

Audit historical Pmetrics 3.0.9 PTA regimen/weight-result association.

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
