# FDA-FEP Public Analysis

Cefepime population pharmacokinetic modeling code for renal replacement therapy, plasma, post-filter plasma, effluent, and epithelial lining fluid analyses.

- FDA BAA 75F40122C00134 | PI: Jim Rhodes, PharmD | Midwestern University

---

## Repository structure

```text
FDA_FEP_public/
├── FDA_FEP_public.Rproj
├── Pmetrics/
│   ├── Rscript/
│   │   ├── Analysis.R             # canonical manuscript workflow
│   │   └── Reviewer_Analyses.R    # isolated reviewer-response analyses
│   ├── Sim/
│   │   ├── sim5.csv               # continuous infusion template
│   │   ├── sim6.csv               # extended infusion template
│   │   └── sim7.csv               # intermittent infusion template
│   ├── Runs/                      # only HTML reports may be tracked
│   ├── src/                       # local model-ready data; ignored by Git
│   └── private/                   # local preprocessing code/data; ignored by Git
├── .gitignore
└── README.md
```

The simulation templates contain no patient-level observations and are intentionally public. Clinical source data, model-ready patient data, private preprocessing code, identifiers, and non-HTML Pmetrics run outputs are excluded from the public repository.

## Working rules

- Open the project from `FDA_FEP_public.Rproj`.
- Run commands from the repository root.
- `Pmetrics/Rscript/Analysis.R` is the canonical public analysis script.
- Use repository-root-relative paths; the scripts do not call `setwd()`.
- Do not commit anything under `Pmetrics/private/` or `Pmetrics/src/`.
- Commit `Pmetrics/Sim/sim5.csv`, `sim6.csv`, and `sim7.csv`.
- Only HTML artifacts under `Pmetrics/Runs/` may be committed.

## Private preprocessing boundary

The private preprocessing script is not part of the public repository. Locally, it should be stored at:

```text
Pmetrics/private/Prepare_Local_Data.R
```

Exact local layout after preprocessing:

```text
Pmetrics/
├── private/
│   ├── Prepare_Local_Data.R
│   ├── raw/
│   │   └── authorized_source.xlsx
│   └── derived/
│       ├── development.csv
│       └── validation.csv
└── src/
    ├── development.csv
    └── validation.csv
```

The private script writes the retained copies under `Pmetrics/private/derived/`, then copies byte-identical model-ready inputs into `Pmetrics/src/`, which is the location consumed by `Analysis.R`.

All four files are ignored by Git. The `check` action verifies that the retained private copies and the `Pmetrics/src/` copies are byte-identical.

Run preprocessing:

```bash
Rscript --vanilla Pmetrics/private/Prepare_Local_Data.R doctor
Rscript --vanilla Pmetrics/private/Prepare_Local_Data.R all --overwrite
```

An explicit authorized source path may be supplied:

```bash
Rscript --vanilla Pmetrics/private/Prepare_Local_Data.R all \
  --source="Pmetrics/private/raw/authorized_source.xlsx" \
  --overwrite
```

A CSV export of the source sheet is also accepted and avoids the `readxl` dependency.

## Public simulation templates

`Analysis.R` reads Adrian's existing templates directly from `Pmetrics/Sim/`.

| File | Regimen | Weight templates | CRRT flow assumption |
|---|---|---|---|
| `sim5.csv` | 2 g loading dose, then 4 g over 24 h continuous infusion | 60, 80, 120 kg | 30 mL/kg/h: 1.8, 2.4, 3.6 L/h |
| `sim6.csv` | 2 g loading dose, then 2 g q12 h over 4 h | 60, 80, 120 kg | 30 mL/kg/h: 1.8, 2.4, 3.6 L/h |
| `sim7.csv` | 2 g loading dose, then 2 g q12 h over 0.5 h | 60, 80, 120 kg | 30 mL/kg/h: 1.8, 2.4, 3.6 L/h |

## Public run map

Public run numbers follow the manuscript rather than the historical development sequence.

| Public run | Manuscript role | Model structure | Historical provenance |
|---:|---|---|---:|
| 1 | Table S3 candidate 1; selected final model | Piecewise central volume by CRRT status; native systemic clearance | 38 |
| 2 | Table S3 candidate 2 | Native clearance scaled by `WT/70` | 48 |
| 3 | Table S3 candidate 3 | Non-CRRT `Voff` scaled by `WT/70` | 46 |
| 4 | Table S3 candidate 4 | Native clearance scaled by `CrCl/120` | 36 |
| 5 | Held-out MAP validation | Public run 1 used as the fixed informative population prior | 39 |

Historical run numbers are provenance only. The public workflow creates and consumes numeric Pmetrics run folders `Pmetrics/Runs/1/` through `Pmetrics/Runs/5/`.

## Minimum linear analysis path

1. Run the private preprocessing script.
2. Confirm `Pmetrics/src/development.csv` and `Pmetrics/src/validation.csv` exist.
3. Fit public runs 1-4 on the same development dataset.
4. Compare all four candidates in manuscript order.
5. Apply public run 1 to the held-out validation cohort as public run 5.
6. Generate parameter and PTA summaries from public run 1 only.

Reviewer-response analyses are isolated from this primary path and use run numbers beginning at 101.

## Pmetrics 3.2.1 model requirements

The canonical script enforces Pmetrics 3.2.1. All five assay/model-error declarations explicitly map to `outeq = 1` through `outeq = 5`; no output relies on list position or an implicit default.

The `check` command compiles every primary and reviewer model definition before any fit is attempted.

## Commands

Use the R installation under which Pmetrics 3.2.1 and its compiled dependencies were installed. From Git Bash, the executable may be set explicitly:

```bash
RSCRIPT="/c/Program Files/R/R-4.6.1/bin/x64/Rscript.exe"
```

Validate data, simulation templates, model definitions, the run registry, and the Git allowlist:

```bash
"$RSCRIPT" --vanilla Pmetrics/Rscript/Analysis.R check
```

Run each stage:

```bash
"$RSCRIPT" --vanilla Pmetrics/Rscript/Analysis.R fit
"$RSCRIPT" --vanilla Pmetrics/Rscript/Analysis.R compare
"$RSCRIPT" --vanilla Pmetrics/Rscript/Analysis.R validate
"$RSCRIPT" --vanilla Pmetrics/Rscript/Analysis.R simulate
```

Run the complete manuscript path:

```bash
"$RSCRIPT" --vanilla Pmetrics/Rscript/Analysis.R all
```

Alternative model-ready data paths remain available for controlled local use:

```bash
"$RSCRIPT" --vanilla Pmetrics/Rscript/Analysis.R all \
  --development="D:/authorized/development.csv" \
  --validation="D:/authorized/validation.csv"
```

Existing run folders are protected by default. Add `--overwrite` only for an intentional rerun.

## Reviewer-response analyses

`Reviewer_Analyses.R` is loaded only through the canonical script and does not alter public runs 1-5.

```bash
"$RSCRIPT" --vanilla Pmetrics/Rscript/Analysis.R reviewer audit
"$RSCRIPT" --vanilla Pmetrics/Rscript/Analysis.R reviewer fit
"$RSCRIPT" --vanilla Pmetrics/Rscript/Analysis.R reviewer report
"$RSCRIPT" --vanilla Pmetrics/Rscript/Analysis.R reviewer all
```

Reviewer-only runs are:

| Reviewer run | Analysis | Purpose |
|---:|---|---|
| 101 | Single-volume development fit | Replaces `Von/Voff` with one central-volume parameter while retaining K12/K21 and the remaining model structure |
| 102 | Single-volume held-out MAP validation | Predictive comparison with public run 5 |
| 103 | Exclude HD-only subjects | Tests robustness of the selected model and `CL1` to HD-only inclusion |
| 104 | Fixed `CL_HD = 3.6 L/h` | Lower intermittent-HD clearance sensitivity |
| 105 | Fixed `CL_HD = 10.8 L/h` | Upper intermittent-HD clearance sensitivity |
| 106 | Remove the two highest development ELF observations | Reviewer-requested influence analysis |
| 107 | Conventional single-volume covariate development fit | One central volume scaled by `WT/70` and native clearance scaled by `CrCl/120` |
| 108 | Conventional single-volume held-out MAP validation | Predictive comparison with public run 5 and reviewer run 102 |

The reviewer reports address code-resolvable questions concerning cohort and RRT composition, CRRT exposure, HD timing, simulation flow assumptions, structural sensitivity, HD-only sensitivity, HD-clearance sensitivity, ELF influence, and expanded goodness-of-fit diagnostics. Laboratory dilution integrity, BAL procedural review, breakpoint standards, and infusion-stability evidence require source documentation outside the model code.

## Scientific mapping

The selected model estimates 10 parameters: `Von`, `Voff`, `V2` (reported as `V_ELF`), `K12`, `K21`, `K15`, `K51`, `CL1`, `S_eff`, and `S_post`. It predicts pre-filter plasma, post-filter plasma, effluent concentration, cumulative effluent amount, and ELF concentration.

Held-out validation uses the selected development population distribution as the informative prior without re-estimating population support-point locations.

Monte Carlo simulation uses 1,000 semiparametric profiles per regimen and weight template, predictions from 23.9 to 48 hours every 0.1 hour, a fixed free fraction of 0.8, MICs from 0.25 to 32 mg/L, and 50%, 68%, and 100% `fT>MIC` targets.

## Public repository scope

This repository contains public analysis code, documentation, the three non-patient simulation templates, and selected HTML artifacts. Clinical data, derived patient-level datasets, private preprocessing logic, and non-HTML analysis outputs are intentionally excluded.
