# FDA-FEP Public Analysis

Cefepime population pharmacokinetic modeling code for CRRT, plasma, post-filter, effluent, and epithelial lining fluid analyses.

- FDA BAA 75F40122C00134 | PI: Jim Rhodes, PharmD | Midwestern University

---

## Repository structure

```text
FDA_FEP_public/
├── FDA_FEP_public.Rproj
├── Pmetrics/
│   ├── Rscript/
│   │   └── Analysis.R
│   ├── Runs/       # local model runs; only HTML reports are tracked
│   ├── Sim/        # local simulation inputs; not tracked
│   └── src/        # local source data; not tracked
├── .gitignore
└── README.md
```

## Working rules

- Open the project from `FDA_FEP_public.Rproj`.
- `Pmetrics/Rscript/Analysis.R` is the canonical analysis script.
- Use repository-root-relative paths.
- Do not commit source data, derived data, simulation inputs, or model-run files.
- HTML reports generated inside `Pmetrics/Runs/` may be committed.
- Do not commit patient-level or otherwise sensitive information.

## Local files

The analysis requires data and supporting files that are not distributed with this repository. Place authorized local files under `Pmetrics/src/`, `Pmetrics/Sim/`, or the appropriate run directory.

## Reproducibility

Run the analysis from the repository root. Model settings, run numbering, data preparation, validation, and report-generation code should be documented in `Analysis.R`.

## Public repository scope

This repository contains analysis code, project documentation, and selected HTML artifacts only. Data and non-HTML analysis outputs are intentionally excluded.
