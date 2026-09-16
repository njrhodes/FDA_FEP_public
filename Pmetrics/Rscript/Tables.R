# =============================================================================
# FDA-FEP manuscript tables
# Produces manuscript Tables 1-2 and Supplemental Tables S1-S4.
# =============================================================================

TABLE_SCRIPT_VERSION <- "1.2.0"

ROOT <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
ANALYSIS_SCRIPT <- file.path("Pmetrics", "Rscript", "Analysis.R")
TABLE_DIR <- file.path("Manuscript", "Tables")

if (!file.exists(ANALYSIS_SCRIPT)) {
  stop("Run Tables.R from the FDA_FEP_public repository root.", call. = FALSE)
}

dir.create(TABLE_DIR, recursive = TRUE, showWarnings = FALSE)
source(ANALYSIS_SCRIPT, local = .GlobalEnv)

required_packages <- c("flextable", "officer")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1L), quietly = TRUE)
]
if (length(missing_packages)) {
  stop("Install required package(s): ", paste(missing_packages, collapse = ", "), call. = FALSE)
}

assert_pmetrics()

as_numeric_safe <- function(x) suppressWarnings(as.numeric(as.character(x)))

save_table <- function(data, stem, title, footer = NULL, digits = NULL) {
  csv_path <- file.path(TABLE_DIR, paste0(stem, ".csv"))
  docx_path <- file.path(TABLE_DIR, paste0(stem, ".docx"))

  utils::write.csv(data, csv_path, row.names = FALSE, na = "")

  ft <- flextable::flextable(data)
  ft <- flextable::add_header_lines(ft, values = title)
  ft <- flextable::font(ft, fontname = "Arial", part = "all")
  ft <- flextable::fontsize(ft, size = 9, part = "all")
  ft <- flextable::bold(ft, part = "header")
  ft <- flextable::bold(ft, i = 1, part = "header")
  ft <- flextable::align(ft, align = "center", part = "all")
  ft <- flextable::align(ft, j = 1, align = "left", part = "body")
  ft <- flextable::valign(ft, valign = "center", part = "all")
  ft <- flextable::border_remove(ft)
  ft <- flextable::hline_top(
    ft, part = "header",
    border = officer::fp_border(color = "black", width = 1)
  )
  ft <- flextable::hline_bottom(
    ft, part = "header",
    border = officer::fp_border(color = "black", width = 0.8)
  )
  ft <- flextable::hline_bottom(
    ft, part = "body",
    border = officer::fp_border(color = "black", width = 1)
  )
  ft <- flextable::autofit(ft)

  if (!is.null(footer) && nzchar(footer)) {
    ft <- flextable::add_footer_lines(ft, values = footer)
    ft <- flextable::fontsize(ft, size = 8, part = "footer")
    ft <- flextable::align(ft, align = "left", part = "footer")
  }

  flextable::save_as_docx(`Manuscript table` = ft, path = docx_path)
  message("Wrote ", csv_path)
  message("Wrote ", docx_path)
  invisible(list(data = data, flextable = ft, csv = csv_path, docx = docx_path))
}

prepare_standard_data <- function(run) {
  dat <- as.data.frame(run$data$standard_data, check.names = FALSE)
  names(dat) <- tolower(names(dat))
  numeric_columns <- intersect(
    c("id", "time", "age", "ht", "wt", "scr", "crcl", "bsa", "crcl_bsa",
      "male", "ecmo", "hd", "crrt", "cvvh", "cvvhd", "cvvhdf", "flow",
      "outeq", "out"),
    names(dat)
  )
  dat[numeric_columns] <- lapply(dat[numeric_columns], as_numeric_safe)
  dat <- dat[is.finite(dat$id), , drop = FALSE]
  dat[order(dat$id, dat$time), , drop = FALSE]
}

first_subject_rows <- function(dat) {
  dat <- dat[order(dat$id, dat$time), , drop = FALSE]
  dat[!duplicated(dat$id), , drop = FALSE]
}

subject_ever <- function(dat, variable) {
  split_values <- split(dat[[variable]], dat$id)
  sum(vapply(split_values, function(x) any(x == 1, na.rm = TRUE), logical(1L)))
}

format_median_iqr <- function(x, digits = 1L) {
  x <- as_numeric_safe(x)
  x <- x[is.finite(x)]
  q <- stats::quantile(x, c(0.25, 0.50, 0.75), na.rm = TRUE, names = FALSE)
  sprintf(paste0("%.", digits, "f (%.", digits, "f - %.", digits, "f)"), q[2], q[1], q[3])
}

format_count_percent <- function(count, denominator) {
  sprintf("%d (%.1f%%)", count, 100 * count / denominator)
}

# =============================================================================
# Table 1. Baseline characteristics and clinical characteristics
# =============================================================================

make_table1 <- function(write_files = TRUE) {
  development <- prepare_standard_data(load_manuscript_run("base"))
  validation <- prepare_standard_data(load_manuscript_run("final_validation"))
  combined <- rbind(development, validation)
  baseline <- first_subject_rows(combined)
  n_subjects <- length(unique(combined$id))

  if (n_subjects != 51L) {
    stop("Table 1 expects 51 unique subjects; found ", n_subjects, ".", call. = FALSE)
  }

  first_crrt_flow <- vapply(split(combined, combined$id), function(x) {
    index <- which(x$crrt == 1 & is.finite(x$flow))
    if (!length(index)) return(NA_real_)
    x$flow[index[[1L]]] / 1000
  }, numeric(1L))

  male <- subject_ever(combined, "male")
  female <- n_subjects - male

  output <- data.frame(
    `Demographics (N = 51)` = c(
      "Age (years)",
      "Height (cm)",
      "Weight (kg)",
      "SCr (mg/dL)",
      "BSA (m^2)",
      "CrCl (mL/min)",
      "Flow (L/hr) at initial CRRT",
      "ECMO, n (%)",
      "Ever received HD, n (%)",
      "Ever received CRRT, n (%)",
      "CVVH, n (%)",
      "CVVHD, n (%)",
      "CVVHDF, n (%)",
      "Sex, n (%)",
      "  Male",
      "  Female"
    ),
    `Median (IQR); n (%)` = c(
      format_median_iqr(baseline$age, 1L),
      format_median_iqr(baseline$ht, 1L),
      format_median_iqr(baseline$wt, 1L),
      format_median_iqr(baseline$scr, 1L),
      format_median_iqr(baseline$bsa, 1L),
      format_median_iqr(baseline$crcl, 1L),
      format_median_iqr(first_crrt_flow, 1L),
      format_count_percent(subject_ever(combined, "ecmo"), n_subjects),
      format_count_percent(subject_ever(combined, "hd"), n_subjects),
      format_count_percent(subject_ever(combined, "crrt"), n_subjects),
      format_count_percent(subject_ever(combined, "cvvh"), n_subjects),
      format_count_percent(subject_ever(combined, "cvvhd"), n_subjects),
      format_count_percent(subject_ever(combined, "cvvhdf"), n_subjects),
      "",
      format_count_percent(male, n_subjects),
      format_count_percent(female, n_subjects)
    ),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  if (write_files) {
    save_table(
      output,
      "Table1_baseline_characteristics",
      "Table 1. Baseline Characteristics and Clinical Characteristics",
      paste(
        "Continuous variables are presented as median (interquartile range), and categorical values as n (%).",
        "Abbreviations: BSA, body surface area; SCr, serum creatinine; CrCl, creatinine clearance;",
        "ECMO, extracorporeal membrane oxygenation; HD, hemodialysis; CRRT, continuous renal replacement therapy;",
        "CVVH, continuous venovenous hemofiltration; CVVHD, continuous venovenous hemodialysis;",
        "CVVHDF, continuous venovenous hemodiafiltration. Renal replacement categories are not mutually exclusive."
      )
    )
  }
  invisible(output)
}

# =============================================================================
# Table 2. Median parameter values and 95% credible intervals
# =============================================================================

make_table2 <- function(write_files = TRUE) {
  summary <- parameter_summary(write_html = FALSE)
  keep <- c(
    "V_on (L)", "V_off (L)", "V_ELF (L)", "CL_systemic (L/h)",
    "CL_HD (L/h)", "CL_CRRT at 2.8 L/h (L/h)",
    "K12 (h^-1)", "K21 (h^-1)", "K15 (h^-1)", "K51 (h^-1)",
    "S_eff", "S_post"
  )
  summary <- summary[match(keep, summary$parameter), , drop = FALSE]

  output <- data.frame(
    Parameter = sub("CL_CRRT at 2.8 L/h", "CL_CRRT", summary$parameter, fixed = TRUE),
    Median = sprintf("%.2f", summary$median),
    `Non-parametric CV%` = ifelse(is.na(summary$cv_percent), "", sprintf("%.1f", summary$cv_percent)),
    `95% CrI` = ifelse(
      is.na(summary$lower_95),
      ifelse(summary$parameter == "CL_HD (L/h)", "Fixed", "Derived"),
      sprintf("%.2f - %.2f", summary$lower_95, summary$upper_95)
    ),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  if (write_files) {
    save_table(
      output,
      "Table2_population_parameters",
      "Table 2. Median Parameter Values and 95% Credible Intervals from the Final Population PK Model",
      paste(
        "Median, CV%, and 95% credible intervals were calculated from the weighted final population support-point distribution.",
        "V_on and V_off are the CRRT-on and CRRT-off central volume terms; V_ELF is the ELF scaling volume;",
        "CL_systemic is native systemic clearance; CL_HD is fixed at 7.2 L/h; CL_CRRT is derived at the cohort median",
        "initial effluent flow of 2.8 L/h; S_eff is the CRRT effluent clearance scalar; S_post is the post-filter plasma scalar."
      )
    )
  }
  invisible(output)
}

# =============================================================================
# Table S1. Samples by cohort and matrix
#
# Subject-level counts are calculated internally and then rolled up by cohort
# and sampling method. No subject-level rows are written to manuscript output.
#
# The id->sampling_group mapping is read at runtime from the current local
# development/validation CSVs through Analysis.R's get_analysis_inputs()
# interface.  It is never read from frozen Pmetrics run objects (which predate
# that column) and is never hard-coded in this script.
# =============================================================================

# Build a named character vector: id (as character) -> sampling_group label.
# Reads the local development and validation CSVs via get_analysis_inputs().
# The run objects passed to make_table_s1 are not consulted for this mapping.
read_sampling_group_lookup <- function() {
  inputs <- get_analysis_inputs()
  paths <- c(inputs$development, inputs$validation)
  frames <- lapply(paths, function(path) {
    if (!file.exists(path)) return(NULL)
    dat <- utils::read.csv(
      path,
      check.names = FALSE,
      stringsAsFactors = FALSE,
      na.strings = c("", "NA", ".")
    )
    names(dat) <- tolower(names(dat))
    if (!"sampling_group" %in% names(dat)) {
      stop(
        "Local data is missing the required 'sampling_group' column: ", path,
        call. = FALSE
      )
    }
    finite_ids <- is.finite(suppressWarnings(as.numeric(as.character(dat$id))))
    dat[finite_ids, c("id", "sampling_group"), drop = FALSE]
  })
  frames <- frames[!vapply(frames, is.null, logical(1L))]
  if (!length(frames)) stop("No local development or validation data found.", call. = FALSE)
  combined <- do.call(rbind, frames)
  tapply(
    combined$sampling_group,
    as.character(combined$id),
    function(x) x[!is.na(x)][1L]
  )
}

make_table_s1 <- function(write_files = TRUE) {
  outeq_labels <- c(
    "1" = "Pre-filter",
    "2" = "Post-filter",
    "3" = "Effluent concentration",
    "4" = "Effluent mg (volume x conc)",
    "5" = "ELF"
  )

  development <- prepare_standard_data(load_manuscript_run("base"))
  validation <- prepare_standard_data(load_manuscript_run("final_validation"))

  common_columns <- intersect(
    names(development),
    names(validation)
  )

  comb <- rbind(
    development[, common_columns, drop = FALSE],
    validation[, common_columns, drop = FALSE]
  )

  required_columns <- c("id", "outeq", "out")
  missing_columns <- setdiff(required_columns, names(comb))

  if (length(missing_columns)) {
    stop(
      "Table S1 requires column(s): ",
      paste(missing_columns, collapse = ", "),
      call. = FALSE
    )
  }

  sample_rows <- comb[
    !is.na(comb$id) &
      !is.na(comb$outeq),
    c("id", "outeq", "out"),
    drop = FALSE
  ]

  grouped <- split(
    sample_rows,
    interaction(
      sample_rows$id,
      sample_rows$outeq,
      drop = TRUE,
      lex.order = TRUE
    )
  )

  samples_per_id_outeq <- do.call(
    rbind,
    lapply(grouped, function(x) {
      data.frame(
        id = x$id[[1L]],
        outeq = x$outeq[[1L]],
        `n(samples)` = sum(!is.na(x$out)),
        check.names = FALSE,
        stringsAsFactors = FALSE
      )
    })
  )

  rownames(samples_per_id_outeq) <- NULL

  sg_lookup <- read_sampling_group_lookup()
  samples_per_id_outeq$Cohort <- unname(
    sg_lookup[as.character(samples_per_id_outeq$id)]
  )

  samples_per_id_outeq$Method <- unname(
    outeq_labels[
      as.character(samples_per_id_outeq$outeq)
    ]
  )

  unknown_method <- is.na(samples_per_id_outeq$Method)

  samples_per_id_outeq$Method[unknown_method] <- paste0(
    "Unknown outeq ",
    samples_per_id_outeq$outeq[unknown_method]
  )

  cohort_order <- c(
    "Prospective sampling",
    "Opportunistic sampling"
  )

  rollup_groups <- split(
    samples_per_id_outeq,
    interaction(
      factor(
        samples_per_id_outeq$Cohort,
        levels = cohort_order
      ),
      samples_per_id_outeq$Method,
      drop = TRUE,
      lex.order = TRUE
    )
  )

  table_s1_rollup <- do.call(
    rbind,
    lapply(rollup_groups, function(x) {
      data.frame(
        Cohort = x$Cohort[[1L]],
        `Sampling method` = x$Method[[1L]],
        `Subjects (n)` = length(unique(x$id)),
        `Samples (n)` = sum(
          x[["n(samples)"]],
          na.rm = TRUE
        ),
        check.names = FALSE,
        stringsAsFactors = FALSE
      )
    })
  )

  rownames(table_s1_rollup) <- NULL

  cohort_rank <- match(
    table_s1_rollup$Cohort,
    cohort_order
  )

  table_s1_rollup <- table_s1_rollup[
    order(
      cohort_rank,
      table_s1_rollup[["Sampling method"]]
    ),
    ,
    drop = FALSE
  ]

  if (write_files) {
    save_table(
      table_s1_rollup,
      "TableS1_samples_by_cohort",
      paste(
        "Table S1. Samples by Cohort and Matrix",
        ""
      ),
      paste(
        "Sample counts and numbers of unique contributing",
        "subjects are summarized by sampling origin and method.",
        "Totals represent individual samples and are not",
        "normalized per subject.",
        "ELF, epithelial lining fluid."
      )
    )
  }

  invisible(table_s1_rollup)
}

# =============================================================================
# Table S4. Nonparametric support points
# =============================================================================

make_table_s4 <- function(write_files = TRUE) {
  run <- load_manuscript_run("base")
  points <- as.data.frame(run$final$popPoints, check.names = FALSE)
  names(points) <- tolower(names(points))
  probability_index <- find_probability_column(points)
  probability_name <- names(points)[probability_index]
  parameters <- c("von", "voff", "v2", "k12", "k21", "k15", "k51", "cl1", "s_eff", "s_post")

  points <- points[order(as_numeric_safe(points[[probability_name]]), decreasing = TRUE), , drop = FALSE]
  probability <- as_numeric_safe(points[[probability_name]])

  body <- data.frame(
    SupportPoint = seq_len(nrow(points)),
    points[parameters],
    Probability = probability,
    `Cumulative Probability` = cumsum(probability),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  stat_row <- function(label, object, multiplier = 1) {
    x <- as.data.frame(object, check.names = FALSE)
    names(x) <- tolower(names(x))
    data.frame(
      SupportPoint = label,
      x[1L, parameters, drop = FALSE] * multiplier,
      Probability = NA_real_,
      `Cumulative Probability` = NA_real_,
      check.names = FALSE,
      stringsAsFactors = FALSE
    )
  }

  output <- rbind(
    body,
    stat_row("Mean", run$final$popMean),
    stat_row("SD", run$final$popSD),
    stat_row("CV%", run$final$popCV, 100)
  )

  names(output)[2:11] <- c(
    "V_on (L)", "V_off (L)", "V_ELF (L)", "K12 (h^-1)", "K21 (h^-1)",
    "K15 (h^-1)", "K51 (h^-1)", "CL1 (L/h)", "S_eff", "S_post"
  )
  numeric_cols <- setdiff(names(output), "SupportPoint")
  output[numeric_cols] <- lapply(output[numeric_cols], function(x) round(as_numeric_safe(x), 3L))

  if (write_files) {
    save_table(
      output,
      "TableS4_support_points",
      "Table S4. Summary of Nonparametric Population PK Model for Cefepime in Plasma and ELF",
      paste(
        "Rows are final nonparametric support points sorted by descending probability.",
        "Mean, SD, and CV% summarize the final population parameter distribution."
      )
    )
  }
  invisible(output)
}

# =============================================================================
# Table S2. CRRT intensity and effluent flow by analysis cohort
# =============================================================================

format_summary_range <- function(x, digits = 1L) {
  x <- as_numeric_safe(x); x <- x[is.finite(x)]
  if (!length(x)) return("Not available")
  q <- stats::quantile(x, c(0, 0.25, 0.50, 0.75, 1), names = FALSE, na.rm = TRUE)
  sprintf(paste0("%.", digits, "f (IQR %.", digits, "f-%.", digits, "f; range %.", digits, "f-%.", digits, "f)"),
          q[3], q[2], q[4], q[1], q[5])
}

crrt_subject_summary <- function(dat) {
  by_id <- split(dat, dat$id)
  rows <- lapply(by_id, function(x) {
    x <- x[order(x$time), , drop = FALSE]
    state <- unique(x[, c("time", "crrt", "flow", "wt"), drop = FALSE])
    state <- state[is.finite(state$time), , drop = FALSE]
    if (!nrow(state) || !any(state$crrt == 1, na.rm = TRUE)) return(NULL)
    state$crrt[!is.finite(state$crrt)] <- 0
    interval <- c(diff(state$time), 0)
    active <- state$crrt == 1 & interval > 0
    flow_ok <- active & is.finite(state$flow)
    intensity_ok <- flow_ok & is.finite(state$wt) & state$wt > 0
    data.frame(
      id = x$id[[1L]],
      flow_l_h = if (any(flow_ok)) stats::weighted.mean(state$flow[flow_ok] / 1000, interval[flow_ok]) else NA_real_,
      intensity = if (any(intensity_ok)) stats::weighted.mean(state$flow[intensity_ok] / state$wt[intensity_ok], interval[intensity_ok]) else NA_real_,
      duration_h = sum(interval[active], na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  })
  rows <- rows[!vapply(rows, is.null, logical(1L))]
  if (!length(rows)) return(data.frame())
  do.call(rbind, rows)
}

make_table_s2 <- function(write_files = TRUE) {
  build <- function(run, label) {
    dat <- prepare_standard_data(run)
    subject <- crrt_subject_summary(dat)

    data.frame(
      Dataset = label,
      `Ever on CRRT` = nrow(subject),
      `Time-averaged effluent flow, L/h` = format_summary_range(subject$flow_l_h, 1L),
      `Recorded effluent intensity, mL/kg/h` = format_summary_range(subject$intensity, 1L),
      `Observed CRRT duration, h` = format_summary_range(subject$duration_h, 1L),
      check.names = FALSE, stringsAsFactors = FALSE
    )
  }
  output <- rbind(
    build(load_manuscript_run("base"), "Development"),
    build(load_manuscript_run("final_validation"), "Validation")
  )
  if (write_files) save_table(
    output, "TableS2_crrt_intensity",
    "Table S2. CRRT Intensity and Effluent Flow Rates in Development and Validation Cohorts",
    paste("Recorded effluent intensity is recorded absolute flow divided by contemporaneous body weight.",
          "Observed CRRT duration uses the time-varying CRRT indicator over modeled follow-up intervals.")
  )
  invisible(output)
}

# =============================================================================
# Table S3. Focused model comparison
# =============================================================================

make_table_s3 <- function(write_files = TRUE) {
  comparison <- compare_development(write_html = FALSE)

  minus2ll <- as_numeric_safe(comparison$actual_minus2ll)
  aic <- as_numeric_safe(comparison$actual_aic)

  if (!length(minus2ll) || any(!is.finite(minus2ll))) {
    stop(
      "Table S3 requires finite -2LL values for all models.",
      call. = FALSE
    )
  }

  if (!length(aic) || any(!is.finite(aic))) {
    stop(
      "Table S3 requires finite AIC values for all models.",
      call. = FALSE
    )
  }

  reference_minus2ll <- minus2ll[[1L]]
  reference_aic <- aic[[1L]]

  delta_minus2ll <- minus2ll - reference_minus2ll
  delta_aic <- aic - reference_aic

  format_delta <- function(x) {
    ifelse(
      abs(x) < 0.0005,
      "\u2014",
      sprintf("%+.3f", x)
    )
  }

  output <- data.frame(
    Step = seq_len(nrow(comparison)),
    Structure = comparison$structure,
    `No. parameters` = ifelse(
      is.na(comparison$actual_parameters),
      comparison$expected_parameters,
      comparison$actual_parameters
    ),
    `-2LL` = round(minus2ll, 3L),
    AIC = round(aic, 3L),
    `Delta -2LL vs selected model` =
      format_delta(delta_minus2ll),
    `Delta AIC vs selected model` =
      format_delta(delta_aic),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  if (write_files) {
    save_table(
      output,
      "TableS3_model_comparison",
      paste(
        "Table S3. Forward-Inclusion Model Comparison",
        "for Cefepime Population Pharmacokinetic",
        "Model Development"
      ),
      paste(
        "Differences in -2 log-likelihood and AIC were",
        "calculated relative to the selected piecewise",
        "central-volume model in step 1.",
        "Candidate models were compared using -2 log-likelihood,",
        "AIC, diagnostic performance, parameter plausibility,",
        "and intended simulation purpose.",
        "Abbreviations: AIC, Akaike information criterion;",
        "CrCl, creatinine clearance;",
        "CRRT, continuous renal replacement therapy."
      )
    )
  }

  invisible(output)
}

make_main_tables <- function() {
  make_table1(TRUE)
  make_table2(TRUE)
  invisible(TRUE)
}

make_supplemental_tables <- function() {
  make_table_s1(TRUE)
  make_table_s2(TRUE)
  make_table_s3(TRUE)
  make_table_s4(TRUE)
  invisible(TRUE)
}

usage <- function() {
  cat(
    "FDA-FEP Tables.R ", TABLE_SCRIPT_VERSION, "\n\n",
    "Usage:\n",
    "  Rscript Pmetrics/Rscript/Tables.R main\n",
    "  Rscript Pmetrics/Rscript/Tables.R supplemental\n",
    "  Rscript Pmetrics/Rscript/Tables.R all\n",
    sep = ""
  )
}

main <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  if (!length(args) || args[[1L]] %in% c("help", "-h", "--help")) {
    usage()
    return(invisible(TRUE))
  }

  command <- args[[1L]]
  message("FDA-FEP Tables.R ", TABLE_SCRIPT_VERSION)

  if (command == "main") {
    make_main_tables()
  } else if (command == "supplemental") {
    make_supplemental_tables()
  } else if (command == "all") {
    make_main_tables()
    make_supplemental_tables()
  } else {
    usage()
    stop("Unknown command: ", command, call. = FALSE)
  }

  invisible(TRUE)
}

if (sys.nframe() == 0L) main()
