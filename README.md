# FDA-FEP Public Analysis

Public R/Pmetrics code supporting the cefepime population pharmacokinetic analysis in critically ill adults receiving renal replacement therapy.

- FDA BAA 75F40122C00134 | PI: Jim Rhodes, PharmD | Midwestern University

## Repository structure

```text
FDA_FEP_public/
├── FDA_FEP_public.Rproj
├── Pmetrics/
│   ├── Rscript/
│   │   ├── Analysis.R
│   │   └── Reviewer_Analyses.R
│   └── Sim/
│       ├── sim5.csv
│       ├── sim6.csv
│       └── sim7.csv
├── .gitignore
└── README.md
```

The repository is intentionally source-focused. Generated Pmetrics run folders, HTML reports, figures, tables, rendered documents, model-ready clinical data, and local preprocessing materials are not tracked.

## Runtime

The reproducible workflow is pinned to **R 4.4.2** and **Pmetrics 3.0.9**. Model definitions use the Pmetrics 3.0.9 default solver and output-specific error declarations by list order.

Each fitted development or reviewer run uses explicit run-specific settings of `cycles = 1000`, `points = 300`, and `seed = 12345`.

From Git Bash:

```bash
RSCRIPT="/c/Program Files/R/R-4.4.2/bin/x64/Rscript.exe"
"$RSCRIPT" --vanilla -e 'cat(R.version.string, "\n"); cat("Pmetrics ", as.character(packageVersion("Pmetrics")), "\n", sep = "")'
```

## Public simulation templates

`Analysis.R` reads the three non-patient simulation templates in `Pmetrics/Sim/`.

| File | Regimen | Weight templates | CRRT flow assumption |
|---|---|---|---|
| `sim5.csv` | 2 g loading dose, then 4 g over 24 h continuous infusion | 60, 80, 120 kg | 30 mL/kg/h: 1.8, 2.4, 3.6 L/h |
| `sim6.csv` | 2 g loading dose, then 2 g q12 h over 4 h | 60, 80, 120 kg | 30 mL/kg/h: 1.8, 2.4, 3.6 L/h |
| `sim7.csv` | 2 g loading dose, then 2 g q12 h over 0.5 h | 60, 80, 120 kg | 30 mL/kg/h: 1.8, 2.4, 3.6 L/h |

## Primary run map

| Public run | Role | Model structure | Source-run provenance |
|---:|---|---|---:|
| 1 | Selected development model | Piecewise central volume by CRRT status; native systemic clearance | 38 |
| 2 | Development candidate | Native clearance scaled by `WT/70` | 48 |
| 3 | Development candidate | Non-CRRT `Voff` scaled by `WT/70` | 46 |
| 4 | Development candidate | Native clearance scaled by `CrCl/120` | 36 |
| 5 | Held-out MAP validation | Run 1 population distribution used as the fixed informative prior | 39 |

Source-run numbers preserve provenance. Public runs 1-5 are deterministic reproduction runs keyed to the same structures.

## Primary workflow

Validate the runtime, required local inputs, simulation templates, model definitions, run registry, and Git allowlist:

```bash
"$RSCRIPT" --vanilla Pmetrics/Rscript/Analysis.R version
"$RSCRIPT" --vanilla Pmetrics/Rscript/Analysis.R check
```

Run individual stages:

```bash
"$RSCRIPT" --vanilla Pmetrics/Rscript/Analysis.R fit
"$RSCRIPT" --vanilla Pmetrics/Rscript/Analysis.R compare
"$RSCRIPT" --vanilla Pmetrics/Rscript/Analysis.R validate
"$RSCRIPT" --vanilla Pmetrics/Rscript/Analysis.R parameters
"$RSCRIPT" --vanilla Pmetrics/Rscript/Analysis.R simulate
```

Run the complete primary workflow:

```bash
"$RSCRIPT" --vanilla Pmetrics/Rscript/Analysis.R all
```

Existing run folders are protected by default. Use `--overwrite` only for an intentional rerun.

## Reviewer-response analyses

`Reviewer_Analyses.R` is loaded through `Analysis.R` and is isolated from public runs 1-5.

```bash
"$RSCRIPT" --vanilla Pmetrics/Rscript/Analysis.R reviewer audit
"$RSCRIPT" --vanilla Pmetrics/Rscript/Analysis.R reviewer fit
"$RSCRIPT" --vanilla Pmetrics/Rscript/Analysis.R reviewer report
"$RSCRIPT" --vanilla Pmetrics/Rscript/Analysis.R reviewer all
```

Reviewer-only runs are:

| Run | Analysis |
|---:|---|
| 101 | Single-volume development fit |
| 102 | Single-volume held-out MAP validation |
| 103 | Development fit excluding HD-only subjects |
| 104 | Fixed `CL_HD = 3.6 L/h` sensitivity |
| 105 | Fixed `CL_HD = 10.8 L/h` sensitivity |
| 106 | Development fit excluding the two highest ELF observations |
| 107 | Conventional single-volume covariate development fit |
| 108 | Conventional single-volume covariate held-out MAP validation |
| 109 | Development fit with OUTEQ 1-3 values >100 mg/L right-censored at 100 mg/L |
| 110 | Held-out MAP validation using run 109 as prior with the same right-censoring rule |

Runs 109-110 address uncertainty above the stated assay range without asserting dilution integrity. Reviewer-only runs are sensitivity analyses and do not replace or alter the primary run sequence.

## Repository policy

Track source code, project/configuration files, this documentation, and the three public non-patient simulation templates. Do not track generated HTML reports or any other rendered analysis output. Local run directories and local clinical-data inputs remain outside version control.
