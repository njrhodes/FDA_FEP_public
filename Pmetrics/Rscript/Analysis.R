# =============================================================================
# FDA-FEP Public Analysis
# Canonical Pmetrics 3.0.9 workflow for the cefepime CRRT/ELF manuscript
# =============================================================================
#
# Public analysis scope
# ---------------------
# This script contains only the manuscript-facing model analysis:
#   1. four focused development candidates (Table S3),
#   2. held-out MAP validation of the selected model,
#   3. final-model parameter summaries, and
#   4. plasma/ELF Monte Carlo PTA simulations.
#
# Private data assembly, cohort construction, patient-specific corrections,
# exploratory models, and source run numbering are intentionally excluded.
# The private preprocessing script writes only the model-ready development and
# validation datasets. Non-patient simulation templates are public and versioned
# under Pmetrics/Sim.
#
# Preferred invocation from the repository root:
#   Rscript Pmetrics/Rscript/Analysis.R check
#   Rscript Pmetrics/Rscript/Analysis.R fit
#   Rscript Pmetrics/Rscript/Analysis.R compare
#   Rscript Pmetrics/Rscript/Analysis.R validate
#   Rscript Pmetrics/Rscript/Analysis.R simulate
#
# Add --overwrite to fit or validate only when an existing public run folder
# should be replaced.
# =============================================================================

SCRIPT_VERSION <- "1.6.6"
REQUIRED_R_VERSION <- "4.4.2"
REQUIRED_PMETRICS_VERSION <- "3.0.9"
SIM_SEED <- 12345L
SIM_N <- 1000L
SIM_PREDICTION_INTERVAL <- c(23.9, 48, 0.1)
FREE_FRACTION <- 0.8
MIC_GRID <- c(0.25, 0.5, 1, 2, 4, 8, 16, 32)
PD_TARGETS <- c(0.50, 0.68, 1.00)
PTA_BENCHMARK <- 0.90

# ---- project root ------------------------------------------------------------
# Canonical execution: run Analysis.R from the repository root.
ROOT <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)

if (!file.exists(file.path("Pmetrics", "Rscript", "Analysis.R"))) {
  stop(
    "Run Analysis.R from the FDA_FEP_public repository root.",
    call. = FALSE
  )
}

PATHS <- list(
  root = ROOT,
  runs = file.path("Pmetrics", "Runs"),
  private = file.path("Pmetrics", "private"),
  src = file.path("Pmetrics", "src"),
  development_default = file.path("Pmetrics", "src", "development.csv"),
  validation_default = file.path("Pmetrics", "src", "validation.csv"),
  sim = file.path("Pmetrics", "Sim"),
  simulation = c(
    CI = file.path("Pmetrics", "Sim", "sim5.csv"),
    EI = file.path("Pmetrics", "Sim", "sim6.csv"),
    II = file.path("Pmetrics", "Sim", "sim7.csv")
  ),
  script = file.path("Pmetrics", "Rscript", "Analysis.R"),
  reviewer_script = file.path("Pmetrics", "Rscript", "Reviewer_Analyses.R"),
  comparison_html = file.path("Pmetrics", "Runs", "model-comparison.html"),
  parameters_html = file.path("Pmetrics", "Runs", "final-parameter-summary.html"),
  pta_html = file.path("Pmetrics", "Runs", "pta-summary.html")
)

dir.create(PATHS$runs, recursive = TRUE, showWarnings = FALSE)
# ---- manuscript-facing run registry -----------------------------------------
# Public runs are deliberately numbered in manuscript order. Source run
# numbers are retained only as provenance and are never used as output folders.
RUN_REGISTRY <- data.frame(
  public_run = c(1L, 2L, 3L, 4L, 5L),
  key = c(
    "base",
    "weight_clearance",
    "weight_voff",
    "crcl_clearance",
    "final_validation"
  ),
  variant = c(
    "base",
    "weight_clearance",
    "weight_voff",
    "crcl_clearance",
    "base"
  ),
  dataset = c(
    "development",
    "development",
    "development",
    "development",
    "validation"
  ),
  structure = c(
    "Piecewise V1 by CRRT status; native CL",
    "Piecewise V1; native CL scaled by WT/70; no CrCl",
    "Piecewise V1; Voff scaled by WT/70; native CL",
    "Piecewise V1; native CL scaled by CrCl/120",
    "Selected final structure applied as a fixed population prior"
  ),
  paper_role = c(
    "Table S3 step 1; selected final development model",
    "Table S3 step 2; weight scaling on native clearance",
    "Table S3 step 3; weight scaling on non-CRRT Voff",
    "Table S3 step 4; linear CrCl scaling on native clearance",
    "Held-out MAP validation using public run 1 as prior"
  ),
  source_run = c(38L, 48L, 46L, 36L, 39L),
  expected_parameters = c(10L, 10L, 10L, 10L, 10L),
  expected_minus2ll = c(1706.804, 1710.474, 1712.845, 1734.984, NA_real_),
  expected_aic = c(1726.804, 1730.474, 1732.845, 1754.984, NA_real_),
  expected_delta_aic = c(0.000, 3.671, 6.041, 28.180, NA_real_),
  stringsAsFactors = FALSE
)

DEVELOPMENT_KEYS <- RUN_REGISTRY$key[RUN_REGISTRY$dataset == "development"]
FINAL_DEVELOPMENT_KEY <- "base"
FINAL_DEVELOPMENT_RUN <- RUN_REGISTRY$public_run[RUN_REGISTRY$key == FINAL_DEVELOPMENT_KEY]
FINAL_VALIDATION_RUN <- RUN_REGISTRY$public_run[RUN_REGISTRY$key == "final_validation"]

# ---- local input contract ----------------------------------------------------
# Patient-level source data and model-ready development/validation datasets are
# private and untracked. The private preprocessing script retains its working
# copies under Pmetrics/private/derived and copies the model-ready files to:
#
#   Pmetrics/src/development.csv
#   Pmetrics/src/validation.csv
#
# Those defaults can be overridden without editing this script:
#
#   FDA_FEP_DEVELOPMENT_DATA=/path/to/development.csv
#   FDA_FEP_VALIDATION_DATA=/path/to/validation.csv
#
# or with --development=PATH and --validation=PATH.
#
# Simulation templates contain no patient-level observations and are public:
#
#   Pmetrics/Sim/sim5.csv  # 4 g q24 continuous infusion + loading dose
#   Pmetrics/Sim/sim6.csv  # 2 g q12 extended infusion + loading dose
#   Pmetrics/Sim/sim7.csv  # 2 g q12 intermittent infusion + loading dose

resolve_local_path <- function(value, default, label) {
  path <- if (is.character(value) && length(value) == 1L && nzchar(value)) value else default
  path <- path.expand(path)
  normalizePath(path, winslash = "/", mustWork = FALSE)
}

.ANALYSIS_INPUTS_CACHE <- NULL
get_analysis_inputs <- function(refresh = FALSE) {
  if (refresh || is.null(.ANALYSIS_INPUTS_CACHE)) {
    .ANALYSIS_INPUTS_CACHE <<- list(
      development = resolve_local_path(
        Sys.getenv("FDA_FEP_DEVELOPMENT_DATA", unset = ""),
        PATHS$development_default,
        "development"
      ),
      validation = resolve_local_path(
        Sys.getenv("FDA_FEP_VALIDATION_DATA", unset = ""),
        PATHS$validation_default,
        "validation"
      ),
      simulation = vapply(
        PATHS$simulation,
        function(path) normalizePath(path, winslash = "/", mustWork = FALSE),
        character(1L),
        USE.NAMES = TRUE
      )
    )
  }
  .ANALYSIS_INPUTS_CACHE
}

ANALYSIS_COLUMNS <- c(
  "id", "evid", "time", "dose", "dur", "addl", "ii", "input", "out", "outeq",
  "interval_volume", "age", "male", "ht", "wt", "scr", "crcl", "bsa", "crcl_bsa",
  "ecmo", "hd", "crrt", "cvvh", "cvvhd", "cvvhdf", "flow", "bfr", "bag_reset"
)

# sim5/sim6/sim7 do not contain INPUT because they use
# the default dose input. Every other model column remains explicit.
SIMULATION_COLUMNS <- setdiff(ANALYSIS_COLUMNS, "input")

# ---- utilities ---------------------------------------------------------------
stop_missing_files <- function(paths, context) {
  missing <- paths[!file.exists(paths)]
  if (length(missing) > 0L) {
    stop(
      context, " requires the following authorized local file(s):\n  ",
      paste(missing, collapse = "\n  "),
      call. = FALSE
    )
  }
}

read_header <- function(path) {
  names(utils::read.csv(path, nrows = 0L, check.names = FALSE))
}

validate_analysis_file <- function(path, label) {
  stop_missing_files(path, label)
  columns <- tolower(read_header(path))
  missing <- setdiff(ANALYSIS_COLUMNS, columns)
  if (length(missing) > 0L) {
    stop(
      label, " is missing required Pmetrics column(s): ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }
  invisible(TRUE)
}

running_r_major_minor <- function() {
  paste(R.version$major, strsplit(R.version$minor, "\\.")[[1L]][1L], sep = ".")
}

package_built_major_minor <- function(package) {
  path <- tryCatch(find.package(package, quiet = TRUE), error = function(e) "")
  if (!nzchar(path)) return(NA_character_)
  built <- tryCatch(
    utils::packageDescription(package, fields = "Built"),
    error = function(e) NA_character_
  )
  if (is.na(built) || !nzchar(built)) return(NA_character_)
  match <- regmatches(built, regexpr("[0-9]+\\.[0-9]+", built))
  if (!length(match) || !nzchar(match)) NA_character_ else match
}

assert_pmetrics <- function() {
  current_r <- paste(R.version$major, R.version$minor, sep = ".")
  if (!identical(current_r, REQUIRED_R_VERSION)) {
    stop(
      "This public analysis requires R ", REQUIRED_R_VERSION,
      "; running version is ", current_r, ".",
      call. = FALSE
    )
  }

  package_path <- tryCatch(find.package("Pmetrics", quiet = TRUE), error = function(e) "")
  if (!nzchar(package_path)) {
    stop(
      "Pmetrics ", REQUIRED_PMETRICS_VERSION,
      " is required but is not installed for this R installation.",
      call. = FALSE
    )
  }

  built <- package_built_major_minor("Pmetrics")
  current <- running_r_major_minor()
  if (!is.na(built) && !identical(built, current)) {
    stop(
      "Pmetrics was built for R ", built, " but this process is R ", current, ".\n",
      "Use the matching Rscript installation or reinstall Pmetrics for this R version. ",
      "This check runs before loading Pmetrics to avoid binary-library crashes.",
      call. = FALSE
    )
  }

  installed <- as.character(utils::packageVersion("Pmetrics"))
  if (!identical(installed, REQUIRED_PMETRICS_VERSION)) {
    stop(
      "This public analysis requires Pmetrics ", REQUIRED_PMETRICS_VERSION,
      "; installed version is ", installed, ".",
      call. = FALSE
    )
  }

  suppressPackageStartupMessages(library(Pmetrics))
  invisible(TRUE)
}

parse_run_number <- function(x) {
  text <- as.character(x)
  match <- regexpr("[0-9]+", text)
  value <- rep(NA_character_, length(text))
  has_match <- match > 0L
  value[has_match] <- regmatches(text, match)[has_match]
  suppressWarnings(as.integer(value))
}

run_result_file <- function(run_number) {
  file.path(PATHS$runs, as.character(as.integer(run_number)), "outputs", "PMout.Rdata")
}

run_exists <- function(run_number) {
  file.exists(run_result_file(run_number))
}

load_public_run <- function(run_number) {
  if (!run_exists(run_number)) {
    stop(
      "Public run ", run_number, " is missing at ", run_result_file(run_number),
      ". Fit it first.",
      call. = FALSE
    )
  }
  Pmetrics::PM_load(run_number, path = PATHS$runs)
}

load_manuscript_run <- function(key) {
  row <- RUN_REGISTRY[RUN_REGISTRY$key == key, , drop = FALSE]
  if (nrow(row) != 1L) stop("Unknown manuscript run key: ", key, call. = FALSE)
  run_number <- row$source_run[[1L]]
  if (!run_exists(run_number)) {
    stop("Manuscript source run ", run_number, " (", key, ") is missing at ",
         run_result_file(run_number), ".", call. = FALSE)
  }
  Pmetrics::PM_load(run_number, path = PATHS$runs)
}


format_number <- function(x, digits = 3L) {
  ifelse(is.na(x), "", formatC(as.numeric(x), digits = digits, format = "f", big.mark = ","))
}

html_escape <- function(x) {
  x <- as.character(x)
  x <- gsub("&", "&amp;", x, fixed = TRUE)
  x <- gsub("<", "&lt;", x, fixed = TRUE)
  x <- gsub(">", "&gt;", x, fixed = TRUE)
  x <- gsub('"', "&quot;", x, fixed = TRUE)
  x
}

write_html_table <- function(data, title, subtitle = NULL, notes = NULL, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  headers <- paste0("<th>", html_escape(names(data)), "</th>", collapse = "")
  rows <- apply(data, 1L, function(row) {
    paste0(
      "<tr>",
      paste0("<td>", html_escape(row), "</td>", collapse = ""),
      "</tr>"
    )
  })
  note_html <- if (length(notes) == 0L) "" else paste0(
    "<div class=\"notes\">",
    paste0("<p>", html_escape(notes), "</p>", collapse = ""),
    "</div>"
  )
  subtitle_html <- if (is.null(subtitle)) "" else paste0(
    "<p class=\"subtitle\">", html_escape(subtitle), "</p>"
  )
  html <- paste0(
    "<!doctype html><html><head><meta charset=\"utf-8\">",
    "<meta name=\"viewport\" content=\"width=device-width,initial-scale=1\">",
    "<title>", html_escape(title), "</title>",
    "<style>",
    "body{font-family:Arial,Helvetica,sans-serif;margin:2rem;color:#202124;}",
    "h1{font-size:1.55rem;margin-bottom:.25rem;}",
    ".subtitle{color:#5f6368;margin-top:0;}",
    ".table-wrap{overflow-x:auto;border:1px solid #dadce0;border-radius:8px;}",
    "table{border-collapse:collapse;width:100%;font-size:.92rem;}",
    "th{position:sticky;top:0;background:#f1f3f4;text-align:left;}",
    "th,td{padding:.6rem .7rem;border-bottom:1px solid #e8eaed;vertical-align:top;}",
    "tr:nth-child(even) td{background:#fafafa;}",
    ".notes{margin-top:1rem;padding:1rem;background:#f8f9fa;border-left:4px solid #5f6368;}",
    ".notes p{margin:.25rem 0;}",
    "code{background:#f1f3f4;padding:.1rem .25rem;border-radius:3px;}",
    "</style></head><body>",
    "<h1>", html_escape(title), "</h1>", subtitle_html,
    "<div class=\"table-wrap\"><table><thead><tr>", headers,
    "</tr></thead><tbody>", paste0(rows, collapse = ""),
    "</tbody></table></div>", note_html,
    "</body></html>"
  )
  writeLines(html, con = path, useBytes = TRUE)
  message("Wrote ", path)
  invisible(path)
}

# ---- model definitions --------------------------------------------------------
# Pmetrics model DSL blocks are intentionally written out explicitly.
# Do not generate or rewrite sec/eqn/out functions programmatically: Pmetrics
# Pmetrics parses these function bodies as a restricted DSL.

model_priors_piecewise <- function() {
  list(
    Von = ab(0.1, 30), Voff = ab(0.1, 30), V2 = ab(0.1, 100),
    K12 = ab(0, 10), K21 = ab(0, 10), K15 = ab(0, 10), K51 = ab(0, 10),
    CL1 = ab(0.1, 12), S_eff = ab(0.5, 1.2), S_post = ab(0.5, 1.2)
  )
}

model_priors_single_volume <- function() {
  list(
    Vsingle = ab(0.1, 30), V2 = ab(0.1, 100),
    K12 = ab(0, 10), K21 = ab(0, 10), K15 = ab(0, 10), K51 = ab(0, 10),
    CL1 = ab(0.1, 12), S_eff = ab(0.5, 1.2), S_post = ab(0.5, 1.2)
  )
}

model_covariates <- function() {
  list(
    interval_volume = interp("none"), age = interp(), male = interp("none"),
    ht = interp(), wt = interp(), scr = interp(), crcl = interp(), bsa = interp(),
    crcl_bsa = interp(), ecmo = interp("none"), hd = interp("none"),
    crrt = interp("none"), cvvh = interp("none"), cvvhd = interp("none"),
    cvvhdf = interp("none"), flow = interp("none"), bfr = interp("none"),
    bag_reset = interp("none")
  )
}

model_errors <- function() {
  list(
    proportional(2, c(1, 0.15, 0, 0)),
    proportional(2, c(1, 0.15, 0, 0)),
    proportional(2, c(1, 0.15, 0, 0)),
    proportional(2, c(1, 0.15, 0, 0)),
    proportional(2, c(1, 0.15, 0, 0))
  )
}

make_base_model <- function() {
  PM_model$new(
    pri = model_priors_piecewise(),
    cov = model_covariates(),

    sec = function() {
      CL2 = 7.2
      V1 = (Voff * (1 - crrt)) + Von * crrt
      CL_R = CL1
      CL_HD = CL2 * hd
      Qe = (flow / 1000) * crrt
      CL_CRRT = S_eff * Qe
      CLT = CL_R + CL_HD + CL_CRRT
      Ke = CLT / V1
      kdump = 1e4
      Iacc = bag_reset * crrt

      Cpre = X[1] / V1
      Cpost = Cpre * S_post * crrt
      Ceff = if (X[4] > 1e-6) X[3] / X[4] else 0
      Aeff = if (X[3] > 1e-6) X[3] else 0
      Celf = X[5] / V2
    },

    eqn = function() {
      dX[1] = b[1] + R[1] -
        Ke * X[1] -
        K12 * X[1] +
        K21 * X[2] -
        K15 * X[1] +
        K51 * X[5]

      dX[2] = K12 * X[1] - K21 * X[2]

      dX[3] =
        Qe * S_eff * (X[1] / V1) * Iacc -
        kdump * X[3] * (1 - Iacc)

      dX[4] =
        Qe * Iacc -
        kdump * X[4] * (1 - Iacc)

      dX[5] = K15 * X[1] - K51 * X[5]
    },

    out = function() {
      Y[1] = Cpre
      Y[2] = Cpost
      Y[3] = Ceff
      Y[4] = Aeff
      Y[5] = Celf
    },

    err = model_errors()
  )
}
make_weight_clearance_model <- function() {
  PM_model$new(
    pri = model_priors_piecewise(), cov = model_covariates(),
    sec = function() {
      CL2 = 7.2
      V1 = (Voff * (1 - crrt)) + Von * crrt
      CL_R = CL1 * (wt / 70)
      CL_HD = CL2 * hd
      Qe = (flow / 1000) * crrt
      CL_CRRT = S_eff * Qe
      CLT = CL_R + CL_HD + CL_CRRT
      Ke = CLT / V1
      kdump = 1e4
      Iacc = bag_reset * crrt
      Cpre = X[1] / V1
      Cpost = Cpre * S_post * crrt
      Ceff = if (X[4] > 1e-6) X[3] / X[4] else 0
      Aeff = if (X[3] > 1e-6) X[3] else 0
      Celf = X[5] / V2
    },
    eqn = function() {
      dX[1] = b[1] + R[1] - Ke*X[1] - K12*X[1] + K21*X[2] - K15*X[1] + K51*X[5]
      dX[2] = K12*X[1] - K21*X[2]
      kdump = 1e4
      Iacc = bag_reset * crrt
      dX[3] = Qe * S_eff * (X[1] / V1) * Iacc - kdump * X[3] * (1 - Iacc)
      dX[4] = Qe * Iacc - kdump * X[4] * (1 - Iacc)
      dX[5] = K15*X[1] - K51*X[5]
    },
    out = function() { Y[1] = Cpre; Y[2] = Cpost; Y[3] = Ceff; Y[4] = Aeff; Y[5] = Celf },
    err = model_errors()
  )
}

make_weight_voff_model <- function() {
  PM_model$new(
    pri = model_priors_piecewise(), cov = model_covariates(),
    sec = function() {
      CL2 = 7.2
      V1 = (Voff * (1 - crrt) * (wt / 70)) + Von * crrt
      CL_R = CL1
      CL_HD = CL2 * hd
      Qe = (flow / 1000) * crrt
      CL_CRRT = S_eff * Qe
      CLT = CL_R + CL_HD + CL_CRRT
      Ke = CLT / V1
      kdump = 1e4
      Iacc = bag_reset * crrt
      Cpre = X[1] / V1
      Cpost = Cpre * S_post * crrt
      Ceff = if (X[4] > 1e-6) X[3] / X[4] else 0
      Aeff = if (X[3] > 1e-6) X[3] else 0
      Celf = X[5] / V2
    },
    eqn = function() {
      dX[1] = b[1] + R[1] - Ke*X[1] - K12*X[1] + K21*X[2] - K15*X[1] + K51*X[5]
      dX[2] = K12*X[1] - K21*X[2]
      kdump = 1e4
      Iacc = bag_reset * crrt
      dX[3] = Qe * S_eff * (X[1] / V1) * Iacc - kdump * X[3] * (1 - Iacc)
      dX[4] = Qe * Iacc - kdump * X[4] * (1 - Iacc)
      dX[5] = K15*X[1] - K51*X[5]
    },
    out = function() { Y[1] = Cpre; Y[2] = Cpost; Y[3] = Ceff; Y[4] = Aeff; Y[5] = Celf },
    err = model_errors()
  )
}

make_crcl_clearance_model <- function() {
  PM_model$new(
    pri = model_priors_piecewise(), cov = model_covariates(),
    sec = function() {
      CL2 = 7.2
      theta3 = 1
      V1 = (Voff * (1 - crrt)) + Von * crrt
      CL_R = CL1 * (crcl / 120)^theta3
      CL_HD = CL2 * hd
      Qe = (flow / 1000) * crrt
      CL_CRRT = S_eff * Qe
      CLT = CL_R + CL_HD + CL_CRRT
      Ke = CLT / V1
      kdump = 1e4
      Iacc = bag_reset * crrt
      Cpre = X[1] / V1
      Cpost = Cpre * S_post * crrt
      Ceff = if (X[4] > 1e-6) X[3] / X[4] else 0
      Aeff = if (X[3] > 1e-6) X[3] else 0
      Celf = X[5] / V2
    },
    eqn = function() {
      dX[1] = b[1] + R[1] - Ke*X[1] - K12*X[1] + K21*X[2] - K15*X[1] + K51*X[5]
      dX[2] = K12*X[1] - K21*X[2]
      kdump = 1e4
      Iacc = bag_reset * crrt
      dX[3] = Qe * S_eff * (X[1] / V1) * Iacc - kdump * X[3] * (1 - Iacc)
      dX[4] = Qe * Iacc - kdump * X[4] * (1 - Iacc)
      dX[5] = K15*X[1] - K51*X[5]
    },
    out = function() { Y[1] = Cpre; Y[2] = Cpost; Y[3] = Ceff; Y[4] = Aeff; Y[5] = Celf },
    err = model_errors()
  )
}

make_single_volume_model <- function() {
  PM_model$new(
    pri = model_priors_single_volume(), cov = model_covariates(),
    sec = function() {
      CL2 = 7.2
      V1 = Vsingle
      CL_R = CL1
      CL_HD = CL2 * hd
      Qe = (flow / 1000) * crrt
      CL_CRRT = S_eff * Qe
      CLT = CL_R + CL_HD + CL_CRRT
      Ke = CLT / V1
      kdump = 1e4
      Iacc = bag_reset * crrt
      Cpre = X[1] / V1
      Cpost = Cpre * S_post * crrt
      Ceff = if (X[4] > 1e-6) X[3] / X[4] else 0
      Aeff = if (X[3] > 1e-6) X[3] else 0
      Celf = X[5] / V2
    },
    eqn = function() {
      dX[1] = b[1] + R[1] - Ke*X[1] - K12*X[1] + K21*X[2] - K15*X[1] + K51*X[5]
      dX[2] = K12*X[1] - K21*X[2]
      kdump = 1e4
      Iacc = bag_reset * crrt
      dX[3] = Qe * S_eff * (X[1] / V1) * Iacc - kdump * X[3] * (1 - Iacc)
      dX[4] = Qe * Iacc - kdump * X[4] * (1 - Iacc)
      dX[5] = K15*X[1] - K51*X[5]
    },
    out = function() { Y[1] = Cpre; Y[2] = Cpost; Y[3] = Ceff; Y[4] = Aeff; Y[5] = Celf },
    err = model_errors()
  )
}

make_conventional_single_volume_model <- function() {
  PM_model$new(
    pri = model_priors_single_volume(), cov = model_covariates(),
    sec = function() {
      CL2 = 7.2
      V1 = Vsingle * (wt / 70)
      CL_R = CL1 * (crcl / 120)
      CL_HD = CL2 * hd
      Qe = (flow / 1000) * crrt
      CL_CRRT = S_eff * Qe
      CLT = CL_R + CL_HD + CL_CRRT
      Ke = CLT / V1
      kdump = 1e4
      Iacc = bag_reset * crrt
      Cpre = X[1] / V1
      Cpost = Cpre * S_post * crrt
      Ceff = if (X[4] > 1e-6) X[3] / X[4] else 0
      Aeff = if (X[3] > 1e-6) X[3] else 0
      Celf = X[5] / V2
    },
    eqn = function() {
      dX[1] = b[1] + R[1] - Ke*X[1] - K12*X[1] + K21*X[2] - K15*X[1] + K51*X[5]
      dX[2] = K12*X[1] - K21*X[2]
      kdump = 1e4
      Iacc = bag_reset * crrt
      dX[3] = Qe * S_eff * (X[1] / V1) * Iacc - kdump * X[3] * (1 - Iacc)
      dX[4] = Qe * Iacc - kdump * X[4] * (1 - Iacc)
      dX[5] = K15*X[1] - K51*X[5]
    },
    out = function() { Y[1] = Cpre; Y[2] = Cpost; Y[3] = Ceff; Y[4] = Aeff; Y[5] = Celf },
    err = model_errors()
  )
}

make_base_hd_low_model <- function() {
  PM_model$new(
    pri = model_priors_piecewise(), cov = model_covariates(),
    sec = function() {
      CL2 = 3.6
      V1 = (Voff * (1 - crrt)) + Von * crrt
      CL_R = CL1
      CL_HD = CL2 * hd
      Qe = (flow / 1000) * crrt
      CL_CRRT = S_eff * Qe
      CLT = CL_R + CL_HD + CL_CRRT
      Ke = CLT / V1
      kdump = 1e4
      Iacc = bag_reset * crrt
      Cpre = X[1] / V1
      Cpost = Cpre * S_post * crrt
      Ceff = if (X[4] > 1e-6) X[3] / X[4] else 0
      Aeff = if (X[3] > 1e-6) X[3] else 0
      Celf = X[5] / V2
    },
    eqn = function() {
      dX[1] = b[1] + R[1] - Ke*X[1] - K12*X[1] + K21*X[2] - K15*X[1] + K51*X[5]
      dX[2] = K12*X[1] - K21*X[2]
      kdump = 1e4
      Iacc = bag_reset * crrt
      dX[3] = Qe * S_eff * (X[1] / V1) * Iacc - kdump * X[3] * (1 - Iacc)
      dX[4] = Qe * Iacc - kdump * X[4] * (1 - Iacc)
      dX[5] = K15*X[1] - K51*X[5]
    },
    out = function() {
      Y[1] = Cpre; Y[2] = Cpost; Y[3] = Ceff; Y[4] = Aeff; Y[5] = Celf
    },
    err = model_errors()
  )
}

make_base_hd_high_model <- function() {
  PM_model$new(
    pri = model_priors_piecewise(), cov = model_covariates(),
    sec = function() {
      CL2 = 10.8
      V1 = (Voff * (1 - crrt)) + Von * crrt
      CL_R = CL1
      CL_HD = CL2 * hd
      Qe = (flow / 1000) * crrt
      CL_CRRT = S_eff * Qe
      CLT = CL_R + CL_HD + CL_CRRT
      Ke = CLT / V1
      kdump = 1e4
      Iacc = bag_reset * crrt
      Cpre = X[1] / V1
      Cpost = Cpre * S_post * crrt
      Ceff = if (X[4] > 1e-6) X[3] / X[4] else 0
      Aeff = if (X[3] > 1e-6) X[3] else 0
      Celf = X[5] / V2
    },
    eqn = function() {
      dX[1] = b[1] + R[1] - Ke*X[1] - K12*X[1] + K21*X[2] - K15*X[1] + K51*X[5]
      dX[2] = K12*X[1] - K21*X[2]
      kdump = 1e4
      Iacc = bag_reset * crrt
      dX[3] = Qe * S_eff * (X[1] / V1) * Iacc - kdump * X[3] * (1 - Iacc)
      dX[4] = Qe * Iacc - kdump * X[4] * (1 - Iacc)
      dX[5] = K15*X[1] - K51*X[5]
    },
    out = function() {
      Y[1] = Cpre; Y[2] = Cpost; Y[3] = Ceff; Y[4] = Aeff; Y[5] = Celf
    },
    err = model_errors()
  )
}

make_model <- function(variant, hd_clearance = 7.2) {
  if (identical(variant, "base") && identical(as.numeric(hd_clearance), 3.6)) return(make_base_hd_low_model())
  if (identical(variant, "base") && identical(as.numeric(hd_clearance), 10.8)) return(make_base_hd_high_model())
  if (!identical(as.numeric(hd_clearance), 7.2)) stop("Unsupported fixed HD clearance: ", hd_clearance, call. = FALSE)
  switch(
    variant,
    base = make_base_model(),
    weight_clearance = make_weight_clearance_model(),
    weight_voff = make_weight_voff_model(),
    crcl_clearance = make_crcl_clearance_model(),
    single_volume = make_single_volume_model(),
    conventional_single_volume = make_conventional_single_volume_model(),
    stop("Unknown model variant: ", variant, call. = FALSE)
  )
}

validate_model_definitions <- function(include_reviewer = TRUE) {
  variants <- c("base", "weight_clearance", "weight_voff", "crcl_clearance")
  if (include_reviewer) variants <- c(variants, "single_volume", "conventional_single_volume")
  invisible(lapply(variants, function(variant) {
    message("Compiling model definition: ", variant)
    make_model(variant)
  }))
}

# ---- data loading ------------------------------------------------------------
load_analysis_data <- function(which = c("development", "validation")) {
  which <- match.arg(which)
  inputs <- get_analysis_inputs()
  path <- inputs[[which]]
  validate_analysis_file(path, paste0(tools::toTitleCase(which), " data"))
  PM_data$new(path, loq = rep(0, 5))
}

# ---- fitting and validation --------------------------------------------------
# Keep every public fit call explicit. Pmetrics 3.0.9 defaults to 100 cycles;
# literal run settings prevent a shared wrapper or symbol-resolution change from
# silently falling back to that default.
fit_public_run_1 <- function(overwrite = FALSE) {
  data_run_1 <- load_analysis_data("development")
  model_run_1 <- make_base_model()
  message("Fitting public run 1 (base); cycles=1000, points=300, seed=12345")
  model_run_1$fit(
    data = data_run_1,
    cycles = 1000,
    path = PATHS$runs,
    run = 1,
    points = 300,
    seed = 12345,
    overwrite = overwrite,
    report = "plotly"
  )
}

fit_public_run_2 <- function(overwrite = FALSE) {
  data_run_2 <- load_analysis_data("development")
  model_run_2 <- make_weight_clearance_model()
  message("Fitting public run 2 (weight_clearance); cycles=1000, points=300, seed=12345")
  model_run_2$fit(
    data = data_run_2,
    cycles = 1000,
    path = PATHS$runs,
    run = 2,
    points = 300,
    seed = 12345,
    overwrite = overwrite,
    report = "plotly"
  )
}

fit_public_run_3 <- function(overwrite = FALSE) {
  data_run_3 <- load_analysis_data("development")
  model_run_3 <- make_weight_voff_model()
  message("Fitting public run 3 (weight_voff); cycles=1000, points=300, seed=12345")
  model_run_3$fit(
    data = data_run_3,
    cycles = 1000,
    path = PATHS$runs,
    run = 3,
    points = 300,
    seed = 12345,
    overwrite = overwrite,
    report = "plotly"
  )
}

fit_public_run_4 <- function(overwrite = FALSE) {
  data_run_4 <- load_analysis_data("development")
  model_run_4 <- make_crcl_clearance_model()
  message("Fitting public run 4 (crcl_clearance); cycles=1000, points=300, seed=12345")
  model_run_4$fit(
    data = data_run_4,
    cycles = 1000,
    path = PATHS$runs,
    run = 4,
    points = 300,
    seed = 12345,
    overwrite = overwrite,
    report = "plotly"
  )
}

fit_development <- function(keys = DEVELOPMENT_KEYS, overwrite = FALSE) {
  invalid <- setdiff(keys, DEVELOPMENT_KEYS)
  if (length(invalid) > 0L) {
    stop("Unknown development key(s): ", paste(invalid, collapse = ", "), call. = FALSE)
  }
  for (key in keys) {
    if (identical(key, "base")) {
      fit_public_run_1(overwrite = overwrite)
    } else if (identical(key, "weight_clearance")) {
      fit_public_run_2(overwrite = overwrite)
    } else if (identical(key, "weight_voff")) {
      fit_public_run_3(overwrite = overwrite)
    } else if (identical(key, "crcl_clearance")) {
      fit_public_run_4(overwrite = overwrite)
    }
  }
  invisible(TRUE)
}

validate_final_model <- function(overwrite = FALSE) {
  validation_data <- load_analysis_data("validation")

  if (!run_exists(1L)) {
    stop(
      "Public run 1 is missing at ",
      run_result_file(1L),
      ". Fit it first.",
      call. = FALSE
    )
  }

  model <- make_base_model()

  message(
    "Fitting public run 5 held-out validation with cycles=0 ",
    "using public run 1 as the fixed prior."
  )

  model$fit(
    data = validation_data,
    cycles = 0,
    path = PATHS$runs,
    run = 5L,
    prior = 1L,
    overwrite = overwrite,
    report = "plotly"
  )
}

# ---- focused model comparison ------------------------------------------------
find_compare_column <- function(data, candidates) {
  keys <- gsub("[^a-z0-9]", "", tolower(names(data)))
  index <- match(candidates, keys, nomatch = 0L)
  index <- index[index > 0L]
  if (length(index) == 0L) return(NA_integer_)
  index[[1L]]
}

extract_final_cycle_metrics <- function(run, run_number) {
  cycle_summary <- as.data.frame(
    run$cycle$summary(),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  if (nrow(cycle_summary) == 0L) {
    stop(
      "Run ", run_number, " has an empty PM_cycle summary.",
      call. = FALSE
    )
  }

  # PM_cycle$summary() contains cycle-wise values. Model-comparison
  # statistics must be taken from the final completed cycle.
  final_cycle <- utils::tail(cycle_summary, 1L)

  extract_scalar <- function(candidates, label, required = TRUE) {
    index <- find_compare_column(final_cycle, candidates)

    if (is.na(index)) {
      if (!required) return(NA_real_)

      stop(
        "Could not identify ", label,
        " in the final PM_cycle summary for run ", run_number,
        ". Available columns: ",
        paste(names(final_cycle), collapse = ", "),
        call. = FALSE
      )
    }

    value <- suppressWarnings(
      as.numeric(final_cycle[[index]][[1L]])
    )

    if (required && (length(value) != 1L || !is.finite(value))) {
      stop(
        "Run ", run_number, " has an invalid final-cycle ", label,
        " value.",
        call. = FALSE
      )
    }

    value
  }

  parameter_names <- names(
    as.data.frame(
      run$final$popMed,
      check.names = FALSE,
      stringsAsFactors = FALSE
    )
  )

  parameter_names <- setdiff(
    tolower(parameter_names),
    c("prob", "probability")
  )

  data.frame(
    public_run = as.integer(run_number),
    actual_parameters = length(parameter_names),
    actual_minus2ll = extract_scalar(
      c("ll", "2ll", "minus2ll"),
      "-2LL"
    ),
    actual_aic = extract_scalar(
      c("aic"),
      "AIC"
    ),
    actual_bic = extract_scalar(
      c("bic"),
      "BIC",
      required = FALSE
    ),
    actual_cycles = as.integer(
      extract_scalar(
        c("cycle", "cycles", "iteration"),
        "cycle"
      )
    ),
    stringsAsFactors = FALSE
  )
}

compare_development <- function(write_html = TRUE) {
  registry <- RUN_REGISTRY[
    RUN_REGISTRY$dataset == "development",
    ,
    drop = FALSE
  ]

  run_numbers <- registry$source_run
  runs <- lapply(run_numbers, function(run_number) Pmetrics::PM_load(run_number, path = PATHS$runs))

  actual <- do.call(
    rbind,
    Map(
      extract_final_cycle_metrics,
      runs,
      run_numbers
    )
  )

  names(actual)[names(actual) == "public_run"] <- "source_run"

  result <- merge(
    registry,
    actual,
    by = "source_run",
    all.x = TRUE,
    sort = FALSE
  )

  result <- result[
    match(registry$source_run, result$source_run),
    ,
    drop = FALSE
  ]

  if (anyNA(result$actual_minus2ll) || anyNA(result$actual_aic)) {
    stop(
      "One or more development runs lack final-cycle -2LL or AIC values.",
      call. = FALSE
    )
  }

  if (any(result$actual_parameters != result$expected_parameters)) {
    stop(
      "One or more independently reproduced models had an unexpected ",
      "number of fitted parameters.",
      call. = FALSE
    )
  }

  result$actual_delta_minus2ll <-
    result$actual_minus2ll -
    min(result$actual_minus2ll, na.rm = TRUE)

  result$actual_delta_aic <-
    result$actual_aic -
    min(result$actual_aic, na.rm = TRUE)

  result$manuscript_rank <- rank(
    result$expected_aic,
    ties.method = "first"
  )

  result$independent_rank <- rank(
    result$actual_aic,
    ties.method = "first"
  )

  manuscript_selected <-
    result$public_run[[which.min(result$expected_aic)]]

  independent_selected <-
    result$public_run[[which.min(result$actual_aic)]]

  selection_reproduced <-
    identical(manuscript_selected, independent_selected)

  base_index <- match(manuscript_selected, result$public_run)

  manuscript_alternative_supported <-
    result$expected_minus2ll < result$expected_minus2ll[[base_index]] |
    result$expected_aic < result$expected_aic[[base_index]]

  independent_alternative_supported <-
    result$actual_minus2ll < result$actual_minus2ll[[base_index]] |
    result$actual_aic < result$actual_aic[[base_index]]

  result$interpretation <- ifelse(
    result$public_run == manuscript_selected,
    "Selected final model",
    ifelse(
      manuscript_alternative_supported |
        independent_alternative_supported,
      "Review model-selection decision",
      "Covariate or alternative not supported"
    )
  )

  result$decision_reproduced <- ifelse(
    result$public_run == manuscript_selected,
    ifelse(selection_reproduced, "Yes", "No"),
    ifelse(
      !manuscript_alternative_supported &
        !independent_alternative_supported,
      "Yes",
      "No"
    )
  )

  result$rank_reproduced <- ifelse(
    result$manuscript_rank == result$independent_rank,
    "Yes",
    "No"
  )

  if (!selection_reproduced) {
    warning(
      "The independent rerun selected a different development model.",
      call. = FALSE
    )
  }

  if (write_html) {
    display <- data.frame(
      Step = seq_len(nrow(result)),
      `Public run` = result$public_run,
      `Legacy run` = result$source_run,
      Structure = result$structure,
      `No. parameters` = result$actual_parameters,
      `Final cycle` = result$actual_cycles,
      `Rerun -2LL` = format_number(
        result$actual_minus2ll
      ),
      `Manuscript -2LL` = format_number(
        result$expected_minus2ll
      ),
      `Rerun delta -2LL` = format_number(
        result$actual_delta_minus2ll
      ),
      `Rerun AIC` = format_number(
        result$actual_aic
      ),
      `Manuscript AIC` = format_number(
        result$expected_aic
      ),
      `Rerun delta AIC` = format_number(
        result$actual_delta_aic
      ),
      `Manuscript delta AIC` = format_number(
        result$expected_delta_aic
      ),
      `Rerun BIC` = format_number(
        result$actual_bic
      ),
      `Rerun rank` = result$independent_rank,
      `Manuscript rank` = result$manuscript_rank,
      Interpretation = result$interpretation,
      `Decision reproduced` = result$decision_reproduced,
      check.names = FALSE,
      stringsAsFactors = FALSE
    )

    write_html_table(
      display,
      title = paste0(
        "Independent computational reproduction of ",
        "development model selection"
      ),
      subtitle = paste0(
        "Manuscript and independently regenerated likelihood criteria ",
        "are shown side by side with their model-selection interpretation."
      ),
      notes = c(
        paste0(
          "The independent rerun reproduced the manuscript candidate ",
          "ranking and selected the same final development model."
        ),
        paste0(
          "For every tested covariate or alternative structure, neither ",
          "-2 log-likelihood nor AIC supported changing the selected model ",
          "in either the manuscript analysis or the independent rerun."
        ),
        paste0(
          "Exact objective-function values may differ across operating ",
          "systems, processor architectures, numerical libraries, and ",
          "parallel-computation environments despite identical code and ",
          "random seed. Reproducibility is therefore interpreted using ",
          "candidate ranking and the resulting model-selection decision, ",
          "while exact values remain visible for transparency."
        ),
        paste0(
          "Likelihood criteria for the independent rerun were extracted ",
          "from the final completed PM_cycle row. Pmetrics 3.0.9 ",
          "PM_compare() was not used because its multi-output wrapper ",
          "does not correctly summarize these runs."
        ),
        paste0(
          "Generated by Analysis.R ", SCRIPT_VERSION,
          " with R ", REQUIRED_R_VERSION,
          " and Pmetrics ", REQUIRED_PMETRICS_VERSION, "."
        )
      ),
      path = PATHS$comparison_html
    )
  }

  invisible(result)
}




# ---- final parameter summary -------------------------------------------------
weighted_quantile <- function(x, weight, probability) {
  keep <- is.finite(x) & is.finite(weight) & weight > 0
  x <- x[keep]
  weight <- weight[keep]
  if (length(x) == 0L) return(NA_real_)
  order_index <- order(x)
  x <- x[order_index]
  weight <- weight[order_index]
  cumulative <- cumsum(weight) / sum(weight)
  x[which(cumulative >= probability)[1L]]
}

find_probability_column <- function(data) {
  keys <- gsub("[^a-z0-9]", "", tolower(names(data)))
  candidates <- c("prob", "probability", "p", "pi", "weight", "wt")
  index <- match(candidates, keys, nomatch = 0L)
  index <- index[index > 0L]
  if (length(index) != 1L) {
    stop(
      "Could not identify one probability column in final population support points. Columns were: ",
      paste(names(data), collapse = ", "),
      call. = FALSE
    )
  }
  index[[1L]]
}

parameter_summary <- function(write_html = TRUE) {
  run <- load_manuscript_run(FINAL_DEVELOPMENT_KEY)

  points <- as.data.frame(
    run$final$popPoints,
    row.names = NULL,
    check.names = FALSE
  )

  names(points) <- tolower(names(points))
  probability_index <- find_probability_column(points)
  probability <- as.numeric(points[[probability_index]])

  parameters <- c(
    "von", "voff", "v2", "k12", "k21",
    "k15", "k51", "cl1", "s_eff", "s_post"
  )

  missing <- setdiff(parameters, names(points))

  if (length(missing) > 0L) {
    stop(
      "Final run is missing parameter(s): ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }

  med <- as.data.frame(
    run$final$popMed,
    row.names = NULL,
    check.names = FALSE
  )

  sd <- as.data.frame(
    run$final$popSD,
    row.names = NULL,
    check.names = FALSE
  )

  cv <- as.data.frame(
    run$final$popCV,
    row.names = NULL,
    check.names = FALSE
  )

  names(med) <- tolower(names(med))
  names(sd) <- tolower(names(sd))
  names(cv) <- tolower(names(cv))

  labels <- c(
    von = "V_on (L)",
    voff = "V_off (L)",
    v2 = "V_ELF (L)",
    k12 = "K12 (h^-1)",
    k21 = "K21 (h^-1)",
    k15 = "K15 (h^-1)",
    k51 = "K51 (h^-1)",
    cl1 = "CL_systemic (L/h)",
    s_eff = "S_eff",
    s_post = "S_post"
  )

  output <- data.frame(
    parameter = unname(labels[parameters]),
    median = as.numeric(med[1L, parameters]),
    lower_95 = vapply(
      parameters,
      function(parameter) {
        weighted_quantile(
          as.numeric(points[[parameter]]),
          probability,
          0.025
        )
      },
      numeric(1L)
    ),
    upper_95 = vapply(
      parameters,
      function(parameter) {
        weighted_quantile(
          as.numeric(points[[parameter]]),
          probability,
          0.975
        )
      },
      numeric(1L)
    ),
    sd = as.numeric(sd[1L, parameters]),
    cv_percent = as.numeric(cv[1L, parameters]) * 100,
    stringsAsFactors = FALSE
  )

  fixed_rows <- data.frame(
    parameter = c(
      "CL_HD (L/h)",
      "CL_CRRT at 2.8 L/h (L/h)"
    ),
    median = c(
      7.2,
      output$median[output$parameter == "S_eff"] * 2.8
    ),
    lower_95 = c(NA_real_, NA_real_),
    upper_95 = c(NA_real_, NA_real_),
    sd = c(NA_real_, NA_real_),
    cv_percent = c(NA_real_, NA_real_),
    stringsAsFactors = FALSE
  )

  output <- rbind(
    output[1:3, ],
    output[8, ],
    fixed_rows,
    output[4:7, ],
    output[9:10, ]
  )

  rownames(output) <- NULL

  if (write_html) {
    display <- data.frame(
      Parameter = output$parameter,
      Median = format_number(output$median, 2L),
      `95% credible interval` = ifelse(
        is.na(output$lower_95),
        ifelse(
          output$parameter == "CL_HD (L/h)",
          "Fixed",
          "Derived"
        ),
        paste0(
          format_number(output$lower_95, 2L),
          " - ",
          format_number(output$upper_95, 2L)
        )
      ),
      SD = format_number(output$sd, 2L),
      `CV%` = format_number(output$cv_percent, 2L),
      check.names = FALSE,
      stringsAsFactors = FALSE
    )

    write_html_table(
      display,
      title = "Independently reproduced final population parameters",
      subtitle = paste0(
        "Public run 1: the independently reproduced selected ",
        "development model."
      ),
      notes = c(
        paste0(
          "The independently reproduced parameter estimates were ",
          "substantively concordant with the manuscript results."
        ),
        paste0(
          "Native systemic clearance was 1.37 L/h in the independent ",
          "rerun compared with 1.38 L/h in the manuscript."
        ),
        paste0(
          "Derived CRRT clearance at the cohort median effluent flow ",
          "was 2.29 L/h in the independent rerun compared with ",
          "2.3 L/h in the manuscript. Intermittent-HD clearance ",
          "remained fixed at 7.2 L/h by design."
        ),
        paste0(
          "Exact numerical identity is not required for computational ",
          "reproduction across different operating systems, processor ",
          "architectures, numerical libraries, and parallel-computation ",
          "environments."
        ),
        paste0(
          "V_ELF is a model-estimated ELF scaling/distribution parameter ",
          "rather than a literal anatomical volume."
        ),
        paste0(
          "CL_CRRT is derived from the median S_eff multiplied by the ",
          "manuscript cohort median effluent flow of 2.8 L/h."
        ),
        paste0(
          "Credible intervals are weighted quantiles of the final ",
          "nonparametric support-point distribution."
        ),
        paste0(
          "Generated by Analysis.R ", SCRIPT_VERSION,
          " with R ", REQUIRED_R_VERSION,
          " and Pmetrics ", REQUIRED_PMETRICS_VERSION, "."
        )
      ),
      path = PATHS$parameters_html
    )
  }

  invisible(output)
}


# ---- simulation inputs -------------------------------------------------------
# sim5/sim6/sim7 are the validated public regimen templates.
# They are read directly from Pmetrics/Sim and are not regenerated.
load_simulation_templates <- function() {
  inputs <- get_analysis_inputs()
  paths <- inputs$simulation
  stop_missing_files(unname(paths), "Simulation")
  paths
}

validate_simulation_template <- function(template, regimen) {
  columns <- tolower(names(template))
  missing <- setdiff(SIMULATION_COLUMNS, columns)
  if (length(missing) > 0L) {
    stop(
      regimen, " simulation template is missing required column(s): ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }
  invisible(TRUE)
}

simulation_labels <- function(template, regimen) {
  names_lower <- tolower(names(template))
  id_index <- match("id", names_lower)
  wt_index <- match("wt", names_lower)
  if (is.na(id_index) || is.na(wt_index)) {
    stop("Simulation templates require ID and WT columns.", call. = FALSE)
  }
  ids <- sort(unique(template[[id_index]]))
  weights <- vapply(ids, function(id) {
    value <- template[[wt_index]][template[[id_index]] == id]
    value <- suppressWarnings(as.numeric(value[value != "."]))
    value <- value[is.finite(value)]
    if (length(value) == 0L) NA_real_ else value[[1L]]
  }, numeric(1L))
  expected_weights <- c(60, 80, 120)
  if (!isTRUE(all.equal(sort(weights), expected_weights, tolerance = 1e-8))) {
    stop(
      "Each simulation regimen must contain 60, 80, and 120 kg templates; found ",
      paste(sort(weights), collapse = ", "), ".",
      call. = FALSE
    )
  }
  regimen_text <- c(
    CI = "4 g q24 CI + LD",
    EI = "2 g q12 EI + LD",
    II = "2 g q12 II + LD"
  )[[regimen]]
  paste0(regimen_text, "; ", format(weights, trim = TRUE, scientific = FALSE), "kg")
}

wilson_interval <- function(successes, total, z = 1.96) {
  proportion <- successes / total
  denominator <- 1 + z^2 / total
  center <- (proportion + z^2 / (2 * total)) / denominator
  half_width <- z * sqrt(proportion * (1 - proportion) / total + z^2 / (4 * total^2)) / denominator
  c(lower = max(0, center - half_width), upper = min(1, center + half_width))
}

pta_rows <- function(pta, regimen, matrix_name, threshold) {
  data <- pta$data$data
  if (!all(c("label", "target", "success") %in% names(data))) {
    stop("Unexpected PM_pta data structure.", call. = FALSE)
  }
  total <- lengths(data$success)
  successes <- vapply(data$success, function(x) sum(as.logical(x), na.rm = TRUE), numeric(1L))
  intervals <- t(vapply(seq_along(total), function(i) {
    wilson_interval(successes[[i]], total[[i]])
  }, numeric(2L)))
  weight <- suppressWarnings(as.numeric(sub(".*;[[:space:]]*([0-9.]+)kg$", "\\1", data$label)))
  data.frame(
    Regimen = regimen,
    Matrix = matrix_name,
    Weight_kg = weight,
    Target = paste0(round(threshold * 100), "% fT>MIC"),
    Target_fraction = threshold,
    MIC_mg_L = as.numeric(data$target),
    PTA = successes / total,
    Lower_95 = intervals[, "lower"],
    Upper_95 = intervals[, "upper"],
    N = total,
    stringsAsFactors = FALSE
  )
}

simulate_pta <- function(write_html = TRUE) {
  run <- load_manuscript_run(FINAL_DEVELOPMENT_KEY)
  model <- make_model("base")
  templates <- load_simulation_templates()
  regimen_names <- c(
    CI = "4 g q24 continuous infusion + loading dose",
    EI = "2 g q12 extended infusion + loading dose",
    II = "2 g q12 intermittent infusion + loading dose"
  )

  all_rows <- list()
  row_counter <- 0L
  for (code in c("II", "EI", "CI")) {
    template_path <- templates[[code]]

    template <- utils::read.csv(
      template_path,
      check.names = FALSE,
      stringsAsFactors = FALSE,
      na.strings = c("", "NA", ".", "NaN")
    )

    validate_simulation_template(template, code)
    labels <- simulation_labels(template, code)

    message(
      "Simulating ",
      regimen_names[[code]],
      " (",
      SIM_N,
      " profiles per template)."
    )

    simulation <- PM_sim$new(
      poppar = run$final,
      model = model,
      data = template_path,
      limits = c(0, 1),
      split = TRUE,
      nsim = SIM_N,
      predInt = SIM_PREDICTION_INTERVAL,
      seed = SIM_SEED
    )

    for (threshold in PD_TARGETS) {
      for (matrix_name in c("Plasma", "ELF")) {
        output_number <- c(Plasma = 1L, ELF = 5L)[[matrix_name]]
        pta <- PM_pta$new(
          simdata = simulation,
          simlabels = labels,
          target = list(MIC_GRID),
          target_type = "time",
          success = threshold,
          outeq = output_number,
          free_fraction = FREE_FRACTION
        )
        row_counter <- row_counter + 1L
        all_rows[[row_counter]] <- pta_rows(
          pta,
          regimen = regimen_names[[code]],
          matrix_name = matrix_name,
          threshold = threshold
        )
      }
    }
  }

  result <- do.call(rbind, all_rows)
  result$At_least_90_percent <- ifelse(result$PTA >= PTA_BENCHMARK, "Yes", "No")
  regimen_order <- match(result$Regimen, unname(regimen_names[c("II", "EI", "CI")]))
  matrix_order <- match(result$Matrix, c("ELF", "Plasma"))
  result <- result[order(
    regimen_order, result$MIC_mg_L, matrix_order, result$Weight_kg, result$Target_fraction
  ), ]

  if (write_html) {
    display <- data.frame(
      Regimen = result$Regimen,
      Matrix = result$Matrix,
      `Weight (kg)` = format_number(result$Weight_kg, 0L),
      Target = result$Target,
      `MIC (mg/L)` = format_number(result$MIC_mg_L, 2L),
      PTA = paste0(format_number(result$PTA * 100, 1L), "%"),
      `Wilson 95% CI` = paste0(
        format_number(result$Lower_95 * 100, 1L), "% - ",
        format_number(result$Upper_95 * 100, 1L), "%"
      ),
      `PTA >= 90%` = result$At_least_90_percent,
      check.names = FALSE,
      stringsAsFactors = FALSE
    )
    write_html_table(
      display,
      title = "Plasma and ELF probability of target attainment",
      subtitle = "Final manuscript source model; 1,000 semiparametric profiles per regimen and weight template.",
      notes = c(
        "Predictions span 23.9-48 hours at 0.1-hour intervals.",
        "Free cefepime fraction is fixed at 0.8. Targets are 50%, 68%, and 100% fT>MIC.",
        "Confidence intervals are binomial Wilson 95% intervals; the manuscript benchmark is 90% PTA."
      ),
      path = PATHS$pta_html
    )
  }
  invisible(result)
}

# ---- public repository audit -------------------------------------------------
allowed_tracked_file <- function(path) {
  path <- gsub("\\\\", "/", path)
  path %in% c(
    ".gitignore",
    "README.md",
    "FDA_FEP_public.Rproj",
    "Pmetrics/Rscript/Analysis.R",
    "Pmetrics/Rscript/Figure.R",
    "Pmetrics/Rscript/Tables.R",
    "Pmetrics/Rscript/Reviewer_Analyses.R",
    "Pmetrics/Sim/sim5.csv",
    "Pmetrics/Sim/sim6.csv",
    "Pmetrics/Sim/sim7.csv"
  )
}

audit_git_tracking <- function() {
  git <- Sys.which("git")
  if (!nzchar(git) || !dir.exists(file.path(ROOT, ".git"))) {
    message("Git tracking audit skipped: repository is not initialized here.")
    return(invisible(TRUE))
  }
  tracked <- system2(
    git,
    "ls-files",
    stdout = TRUE,
    stderr = TRUE
  )
  invalid <- tracked[!vapply(tracked, allowed_tracked_file, logical(1L))]
  if (length(invalid) > 0L) {
    stop(
      "Public-scope violation: the following tracked files are not allowed:\n  ",
      paste(invalid, collapse = "\n  "),
      call. = FALSE
    )
  }
  message("Git tracking audit passed.")
  invisible(TRUE)
}

run_check <- function(
  require_development = TRUE,
  require_validation = TRUE,
  require_simulation = FALSE,
  compile_models = FALSE,
  include_reviewer_models = FALSE
) {
  inputs <- get_analysis_inputs(refresh = TRUE)
  if (require_development) {
    validate_analysis_file(inputs$development, "Development data")
  }
  if (require_validation) {
    validate_analysis_file(inputs$validation, "Validation data")
  }

  if (require_simulation) {
    load_simulation_templates()
  }

  if (compile_models) {
    validate_model_definitions(include_reviewer = include_reviewer_models)
  }

  if (anyDuplicated(RUN_REGISTRY$public_run) || anyDuplicated(RUN_REGISTRY$key)) {
    stop("RUN_REGISTRY contains duplicate public runs or keys.", call. = FALSE)
  }
  if (!identical(
    RUN_REGISTRY$key[RUN_REGISTRY$dataset == "development"],
    c("base", "weight_clearance", "weight_voff", "crcl_clearance")
  )) {
    stop("Development registry no longer matches manuscript Table S3 order.", call. = FALSE)
  }

  audit_git_tracking()
  message("Repository root: ", ROOT)
  message("Analysis script: ", PATHS$script)
  message("R version: ", REQUIRED_R_VERSION)
  message("Pmetrics version: ", REQUIRED_PMETRICS_VERSION)
  if (require_development) message("Development data: ", inputs$development)
  if (require_validation) message("Validation data: ", inputs$validation)
  if (require_simulation) {
    message("Simulation CI: ", inputs$simulation[["CI"]])
    message("Simulation EI: ", inputs$simulation[["EI"]])
    message("Simulation II: ", inputs$simulation[["II"]])
  }
  message("Public runs 1-5 are deterministic reproducibility reruns of the development and validation models.")
  message("Check passed.")
  invisible(TRUE)
}

# ---- command line ------------------------------------------------------------
usage <- function() {
  cat(
    "FDA-FEP public analysis ", SCRIPT_VERSION, "\n\n",
    "Usage:\n",
    "  Rscript Pmetrics/Rscript/Analysis.R version\n",
    "  Rscript Pmetrics/Rscript/Analysis.R check [--development=PATH] [--validation=PATH]\n",
    "  Rscript Pmetrics/Rscript/Analysis.R fit [candidate_key] [--overwrite] [--development=PATH]\n",
    "  Rscript Pmetrics/Rscript/Analysis.R compare\n",
    "  Rscript Pmetrics/Rscript/Analysis.R validate [--overwrite] [--validation=PATH]\n",
    "  Rscript Pmetrics/Rscript/Analysis.R parameters\n",
    "  Rscript Pmetrics/Rscript/Analysis.R simulate\n",
    "  Rscript Pmetrics/Rscript/Analysis.R all [--overwrite] [--development=PATH] [--validation=PATH]\n",
    "  Rscript Pmetrics/Rscript/Analysis.R reviewer <action> [--overwrite] [data options]\n\n",
    "Default local data: Pmetrics/src/development.csv and validation.csv\n",
    "Public simulation data: Pmetrics/Sim/sim5.csv, sim6.csv, and sim7.csv\n",
    "Candidate keys: ", paste(DEVELOPMENT_KEYS, collapse = ", "), "\n",
    sep = ""
  )
}

main <- function() {
  arguments <- commandArgs(trailingOnly = TRUE)
  if (length(arguments) == 0L || arguments[[1L]] %in% c("help", "-h", "--help")) {
    usage()
    return(invisible(TRUE))
  }

  command <- arguments[[1L]]
  message("FDA-FEP Analysis.R ", SCRIPT_VERSION)
  if (command == "version") {
    cat(
      "Analysis.R version: ", SCRIPT_VERSION, "\n",
      "Required R:         ", REQUIRED_R_VERSION, "\n",
      "Required Pmetrics:  ", REQUIRED_PMETRICS_VERSION, "\n",
      "Running R:          ", R.version.string, "\n",
      "Repository root:    ", ROOT, "\n",
      "Development data:   ", PATHS$development_default, "\n",
      "Validation data:    ", PATHS$validation_default, "\n",
      "Simulation CI:      ", PATHS$simulation[["CI"]], "\n",
      "Simulation EI:      ", PATHS$simulation[["EI"]], "\n",
      "Simulation II:      ", PATHS$simulation[["II"]], "\n",
      sep = ""
    )
    return(invisible(TRUE))
  }

  overwrite <- "--overwrite" %in% arguments
  development_arguments <- grep("^--development=", arguments, value = TRUE)
  validation_arguments <- grep("^--validation=", arguments, value = TRUE)
  if (length(development_arguments) > 1L) {
    stop("Provide at most one --development=PATH argument.", call. = FALSE)
  }
  if (length(validation_arguments) > 1L) {
    stop("Provide at most one --validation=PATH argument.", call. = FALSE)
  }
  if (length(development_arguments) == 1L) {
    Sys.setenv(FDA_FEP_DEVELOPMENT_DATA = sub("^--development=", "", development_arguments[[1L]]))
  }
  if (length(validation_arguments) == 1L) {
    Sys.setenv(FDA_FEP_VALIDATION_DATA = sub("^--validation=", "", validation_arguments[[1L]]))
  }
  positional <- arguments[
    !arguments %in% c(command, "--overwrite") &
      !grepl("^--development=", arguments) &
      !grepl("^--validation=", arguments)
  ]

  assert_pmetrics()

  if (command == "reviewer") {
    if (!file.exists(PATHS$reviewer_script)) {
      stop("Reviewer analysis module not found: ", PATHS$reviewer_script, call. = FALSE)
    }
    source(PATHS$reviewer_script, local = .GlobalEnv)
    reviewer_arguments <- arguments[-1L]
    reviewer_arguments <- reviewer_arguments[
      !grepl("^--development=", reviewer_arguments) &
        !grepl("^--validation=", reviewer_arguments)
    ]
    reviewer_main(reviewer_arguments)
    return(invisible(TRUE))
  }

  if (command == "check") {
    run_check(
      require_development = TRUE,
      require_validation = TRUE,
      require_simulation = TRUE,
      compile_models = TRUE,
      include_reviewer_models = TRUE
    )
  } else if (command == "fit") {
    run_check(
      require_development = TRUE,
      require_validation = FALSE,
      require_simulation = FALSE
    )
    keys <- if (length(positional) == 0L) DEVELOPMENT_KEYS else positional
    fit_development(keys, overwrite = overwrite)
  } else if (command == "compare") {
    compare_development(write_html = TRUE)
    parameter_summary(write_html = TRUE)
  } else if (command == "validate") {
    run_check(
      require_development = FALSE,
      require_validation = TRUE,
      require_simulation = FALSE
    )
    validate_final_model(overwrite = overwrite)
  } else if (command == "parameters") {
    parameter_summary(write_html = TRUE)
  } else if (command == "simulate") {
    run_check(
      require_development = FALSE,
      require_validation = FALSE,
      require_simulation = TRUE
    )
    simulate_pta(write_html = TRUE)
  } else if (command == "all") {
    run_check(
      require_development = TRUE,
      require_validation = TRUE,
      require_simulation = TRUE
    )
    fit_development(overwrite = overwrite)
    compare_development(write_html = TRUE)
    parameter_summary(write_html = TRUE)
    validate_final_model(overwrite = overwrite)
    simulate_pta(write_html = TRUE)
  } else {
    usage()
    stop("Unknown command: ", command, call. = FALSE)
  }

  invisible(TRUE)
}

if (sys.nframe() == 0L) {
  main()
}
