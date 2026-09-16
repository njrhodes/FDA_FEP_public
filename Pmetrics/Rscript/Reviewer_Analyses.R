# =============================================================================
# REVIEWER-RESPONSE ANALYSES
# =============================================================================
# This module is sourced only by:
#   Rscript Pmetrics/Rscript/Analysis.R reviewer ...
#
# It does not alter the manuscript-facing linear path (public runs 1-5).
# Reviewer-only runs use numbers >= 101 and answer specific requests for:
#   - cohort/RRT composition and CRRT exposure summaries,
#   - a single-central-volume structural comparator,
#   - exclusion of HD-only subjects,
#   - fixed intermittent-HD clearance sensitivity,
#   - exclusion of the two highest ELF observations,
#   - right-censoring of OUTEQ 1-3 concentrations above the assay ULOQ,
#   - population/posterior goodness-of-fit and weighted residual diagnostics,
#   - empirical post-filter/pre-filter and ELF/plasma paired summaries, and
#   - explicit simulation flow assumptions.
#
# Generated reviewer reports remain local and are ignored by Git; no subject
# identifiers are written to those reports. Derived data remain in memory or
# temporary files. These analyses are sensitivity/response analyses,
# not additional steps in the prespecified manuscript model-selection path.
# =============================================================================


# Reviewer 2 explicitly requested sensitivity analysis of the fixed IHD
# clearance. The reviewer did not specify alternative values.
#
# The primary value is 7.2 L/h. Lower and upper values are deterministic
# analyst-selected -50% and +50% stress-test bounds:
#   7.2 * 0.50 = 3.6 L/h
#   7.2 * 1.50 = 10.8 L/h
#
# These are not empirical confidence limits, fitted estimates, or values
# directly reported by the reviewer or supporting literature.
REVIEWER_HD_REFERENCE <- 7.2
REVIEWER_HD_LOW <- REVIEWER_HD_REFERENCE * 0.50
REVIEWER_HD_HIGH <- REVIEWER_HD_REFERENCE * 1.50

REVIEWER_MODULE_VERSION <- "1.3.0"
REVIEWER_ELF_EXCLUDE_N <- 2L
REVIEWER_ELF_PAIR_WINDOW_H <- 0.5
REVIEWER_ASSAY_ULOQ <- 100

REVIEWER_PATHS <- list(
  audit_html = file.path(PATHS$runs, "reviewer-data-audit.html"),
  sensitivity_html = file.path(PATHS$runs, "reviewer-model-sensitivity.html"),
  diagnostics_html = file.path(PATHS$runs, "reviewer-diagnostics.html"),
  comment_map_html = file.path(PATHS$runs, "reviewer-comment-map.html"),
  index_html = file.path(PATHS$runs, "reviewer-index.html")
)

REVIEWER_RUN_REGISTRY <- data.frame(
  reviewer_run = c(101L, 102L, 103L, 104L, 105L, 106L, 107L, 108L, 109L, 110L),
  key = c(
    "single_volume_development", "single_volume_validation", "no_hd_only",
    "hd_low", "hd_high", "elf_top2_removed",
    "conventional_volume_development", "conventional_volume_validation",
    "uloq_censored_development", "uloq_censored_validation"
  ),
  variant = c(
    "single_volume", "single_volume", "base", "base", "base", "base",
    "conventional_single_volume", "conventional_single_volume", "base", "base"
  ),
  dataset = c(
    "development", "validation", "development", "development", "development", "development",
    "development", "validation", "development", "validation"
  ),
  hd_clearance = c(
    REVIEWER_HD_REFERENCE, REVIEWER_HD_REFERENCE, REVIEWER_HD_REFERENCE,
    REVIEWER_HD_LOW, REVIEWER_HD_HIGH, REVIEWER_HD_REFERENCE,
    REVIEWER_HD_REFERENCE, REVIEWER_HD_REFERENCE, REVIEWER_HD_REFERENCE, REVIEWER_HD_REFERENCE
  ),
  role = c(
    "Single central volume with the same peripheral, ELF, CRRT, HD, and output structure",
    "Held-out MAP validation of the single-volume comparator",
    "Final piecewise model after excluding subjects who received HD but never CRRT",
    "Final piecewise model with fixed HD clearance reduced by 50%",
    "Final piecewise model with fixed HD clearance increased by 50%",
    "Final piecewise model after excluding the two highest development ELF observations",
    "Single central volume with fixed WT/70 volume and CrCl/120 native-clearance scaling",
    "Held-out MAP validation of the conventional single-volume covariate comparator",
    "Final piecewise model with OUTEQ 1-3 values above 100 mg/L right-censored at the assay ULOQ",
    "Held-out MAP validation of the ULOQ-censored development prior using the same censoring rule"
  ),
  stringsAsFactors = FALSE
)

reviewer_registry_row <- function(key) {
  row <- REVIEWER_RUN_REGISTRY[REVIEWER_RUN_REGISTRY$key == key, , drop = FALSE]
  if (nrow(row) != 1L) stop("Unknown reviewer run key: ", key, call. = FALSE)
  row
}

reviewer_run_exists <- function(run_number) {
  run_exists(run_number)
}

load_reviewer_run <- function(run_number) {
  if (!reviewer_run_exists(run_number)) {
    stop(
      "Reviewer run ", run_number, " is missing. Run `reviewer fit` first.",
      call. = FALSE
    )
  }
  Pmetrics::PM_load(run_number, path = PATHS$runs)
}

# ---- local data utilities ----------------------------------------------------
read_local_analysis_frame <- function(which = c("development", "validation")) {
  which <- match.arg(which)
  path <- get_analysis_inputs()[[which]]
  validate_analysis_file(path, paste0(tools::toTitleCase(which), " data"))
  data <- utils::read.csv(
    path,
    check.names = FALSE,
    stringsAsFactors = FALSE,
    na.strings = c("", "NA", ".")
  )
  names(data) <- tolower(names(data))
  data
}

as_number <- function(x) suppressWarnings(as.numeric(as.character(x)))
as_flag <- function(x) as.integer(as_number(x) == 1)

numeric_if_present <- function(data, columns) {
  for (name in intersect(columns, names(data))) data[[name]] <- as_number(data[[name]])
  data
}

normalize_analysis_frame <- function(data) {
  numeric_if_present(
    data,
    c(
      "time", "evid", "dose", "dur", "addl", "ii", "input", "out", "outeq",
      "interval_volume", "age", "male", "ht", "wt", "scr", "crcl", "bsa",
      "crcl_bsa", "ecmo", "hd", "crrt", "cvvh", "cvvhd", "cvvhdf", "flow",
      "bfr", "bag_reset"
    )
  )
}

subject_modality <- function(data) {
  data <- normalize_analysis_frame(data)
  ids <- unique(data$id)
  output <- lapply(ids, function(id) {
    rows <- data[data$id == id, , drop = FALSE]
    ever_crrt <- any(rows$crrt == 1, na.rm = TRUE)
    ever_hd <- any(rows$hd == 1, na.rm = TRUE)
    category <- if (ever_crrt && ever_hd) {
      "Both CRRT and HD"
    } else if (ever_crrt) {
      "CRRT only"
    } else if (ever_hd) {
      "HD only"
    } else {
      "Neither recorded"
    }
    data.frame(
      id = id,
      ever_crrt = ever_crrt,
      ever_hd = ever_hd,
      modality = category,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, output)
}

observation_rows <- function(data) {
  data <- normalize_analysis_frame(data)
  data[is.finite(data$out) & is.finite(data$outeq) & data$outeq %in% 1:5, , drop = FALSE]
}

matrix_label <- function(outeq) {
  labels <- c(
    `1` = "Pre-filter plasma",
    `2` = "Post-filter plasma",
    `3` = "Effluent concentration",
    `4` = "Cumulative effluent amount",
    `5` = "ELF"
  )
  unname(labels[as.character(outeq)])
}

summary_numbers <- function(x, digits = 2L) {
  x <- as_number(x)
  x <- x[is.finite(x)]
  if (length(x) == 0L) return("Not available")
  q <- stats::quantile(x, c(0, 0.25, 0.5, 0.75, 1), names = FALSE, na.rm = TRUE)
  sprintf(
    paste0("%.", digits, "f (IQR %.", digits, "f-%.", digits, "f; range %.", digits, "f-%.", digits, "f)"),
    q[[3L]], q[[2L]], q[[4L]], q[[1L]], q[[5L]]
  )
}

first_subject_rows <- function(data) {
  data <- normalize_analysis_frame(data)
  split_data <- split(data, data$id)
  output <- lapply(split_data, function(rows) {
    rows <- rows[order(rows$time), , drop = FALSE]
    rows[1L, , drop = FALSE]
  })
  do.call(rbind, output)
}

observed_state_duration <- function(data, state_column) {
  data <- normalize_analysis_frame(data)
  split_data <- split(data, data$id)
  output <- lapply(split_data, function(rows) {
    collapsed <- stats::aggregate(
      rows[[state_column]],
      by = list(time = rows$time),
      FUN = function(x) max(as_number(x), na.rm = TRUE)
    )
    names(collapsed)[[2L]] <- "state"
    collapsed <- collapsed[is.finite(collapsed$time), , drop = FALSE]
    collapsed <- collapsed[order(collapsed$time), , drop = FALSE]
    if (nrow(collapsed) == 0L) {
      return(data.frame(id = rows$id[[1L]], active_h = NA_real_, span_h = NA_real_, sessions = 0L))
    }
    collapsed$state[!is.finite(collapsed$state)] <- 0
    interval <- c(diff(collapsed$time), 0)
    active_h <- sum(interval[collapsed$state == 1], na.rm = TRUE)
    active_times <- collapsed$time[collapsed$state == 1]
    span_h <- if (length(active_times) < 2L) 0 else max(active_times) - min(active_times)
    previous <- c(0, head(collapsed$state, -1L))
    sessions <- sum(collapsed$state == 1 & previous != 1, na.rm = TRUE)
    data.frame(
      id = rows$id[[1L]],
      active_h = active_h,
      span_h = span_h,
      sessions = sessions,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, output)
}

find_source_cohort_column <- function(data) {
  candidates <- c("source_cohort", "sampling_cohort", "cohort", "source", "origin")
  found <- intersect(candidates, names(data))
  if (length(found) == 0L) NA_character_ else found[[1L]]
}

# ---- HTML helpers ------------------------------------------------------------
data_frame_html <- function(data) {
  if (is.null(data) || nrow(data) == 0L) return("<p><em>No rows available.</em></p>")
  headers <- paste0("<th>", html_escape(names(data)), "</th>", collapse = "")
  rows <- apply(data, 1L, function(row) {
    paste0("<tr>", paste0("<td>", html_escape(row), "</td>", collapse = ""), "</tr>")
  })
  paste0(
    "<div class=\"table-wrap\"><table><thead><tr>", headers,
    "</tr></thead><tbody>", paste0(rows, collapse = ""), "</tbody></table></div>"
  )
}

write_reviewer_html <- function(path, title, subtitle, sections, notes = character()) {
  section_html <- vapply(sections, function(section) {
    paste0(
      "<section><h2>", html_escape(section$title), "</h2>",
      if (!is.null(section$intro)) paste0("<p>", html_escape(section$intro), "</p>") else "",
      section$html,
      "</section>"
    )
  }, character(1L))
  note_html <- if (length(notes) == 0L) "" else paste0(
    "<aside><h2>Interpretation notes</h2>",
    paste0("<p>", html_escape(notes), "</p>", collapse = ""),
    "</aside>"
  )
  html <- paste0(
    "<!doctype html><html><head><meta charset=\"utf-8\">",
    "<meta name=\"viewport\" content=\"width=device-width,initial-scale=1\">",
    "<title>", html_escape(title), "</title><style>",
    "body{font-family:Arial,Helvetica,sans-serif;margin:2rem auto;max-width:1280px;color:#202124;line-height:1.45;padding:0 1rem}",
    "h1{font-size:1.7rem;margin-bottom:.25rem}h2{font-size:1.25rem;margin-top:2rem;border-bottom:2px solid #dadce0;padding-bottom:.3rem}",
    ".subtitle{color:#5f6368;margin-top:0}.table-wrap{overflow-x:auto;border:1px solid #dadce0;border-radius:8px;margin:.75rem 0 1.25rem}",
    "table{border-collapse:collapse;width:100%;font-size:.9rem}th{background:#f1f3f4;text-align:left}th,td{padding:.55rem .65rem;border-bottom:1px solid #e8eaed;vertical-align:top}",
    "tr:nth-child(even) td{background:#fafafa}.plots{display:grid;grid-template-columns:repeat(auto-fit,minmax(360px,1fr));gap:1rem}",
    ".plot{border:1px solid #dadce0;border-radius:8px;padding:.5rem;background:white}.plot h3{font-size:1rem;margin:.25rem .4rem}",
    "aside{margin-top:2rem;padding:1rem;background:#f8f9fa;border-left:4px solid #5f6368}code{background:#f1f3f4;padding:.1rem .25rem;border-radius:3px}",
    "a{color:#185abc}svg text{font-family:Arial,Helvetica,sans-serif}",
    "</style></head><body><h1>", html_escape(title), "</h1><p class=\"subtitle\">",
    html_escape(subtitle), "</p>", paste0(section_html, collapse = ""), note_html,
    "</body></html>"
  )
  writeLines(html, path, useBytes = TRUE)
  message("Wrote ", path)
  invisible(path)
}

scatter_svg <- function(x, y, title, xlab, ylab, reference = c("none", "identity", "zero")) {
  reference <- match.arg(reference)
  keep <- is.finite(x) & is.finite(y)
  x <- as.numeric(x[keep]); y <- as.numeric(y[keep])
  if (length(x) == 0L) return(paste0("<div class=\"plot\"><h3>", html_escape(title), "</h3><p>No finite observations.</p></div>"))
  width <- 520; height <- 360; left <- 62; right <- 18; top <- 32; bottom <- 52
  xr <- range(x, na.rm = TRUE); yr <- range(y, na.rm = TRUE)
  if (reference == "identity") {
    both <- range(c(xr, yr), na.rm = TRUE); xr <- both; yr <- both
  }
  expand_range <- function(r) {
    if (!all(is.finite(r))) return(c(0, 1))
    if (diff(r) == 0) return(r + c(-0.5, 0.5))
    r + c(-0.05, 0.05) * diff(r)
  }
  xr <- expand_range(xr); yr <- expand_range(yr)
  sx <- function(value) left + (value - xr[[1L]]) / diff(xr) * (width - left - right)
  sy <- function(value) height - bottom - (value - yr[[1L]]) / diff(yr) * (height - top - bottom)
  points <- paste0(
    "<circle cx=\"", format(sx(x), digits = 8), "\" cy=\"", format(sy(y), digits = 8),
    "\" r=\"3\" fill=\"#1a73e8\" fill-opacity=\"0.55\"/>",
    collapse = ""
  )
  ref <- ""
  if (reference == "identity") {
    low <- max(xr[[1L]], yr[[1L]]); high <- min(xr[[2L]], yr[[2L]])
    ref <- sprintf("<line x1=\"%.2f\" y1=\"%.2f\" x2=\"%.2f\" y2=\"%.2f\" stroke=\"#5f6368\" stroke-dasharray=\"6 5\"/>", sx(low), sy(low), sx(high), sy(high))
  } else if (reference == "zero" && yr[[1L]] <= 0 && yr[[2L]] >= 0) {
    ref <- sprintf("<line x1=\"%d\" y1=\"%.2f\" x2=\"%d\" y2=\"%.2f\" stroke=\"#5f6368\" stroke-dasharray=\"6 5\"/>", left, sy(0), width-right, sy(0))
  }
  xticks <- seq(xr[[1L]], xr[[2L]], length.out = 5L)
  yticks <- seq(yr[[1L]], yr[[2L]], length.out = 5L)
  x_tick_svg <- paste0(sprintf(
    "<line x1=\"%.2f\" y1=\"%d\" x2=\"%.2f\" y2=\"%d\" stroke=\"#9aa0a6\"/><text x=\"%.2f\" y=\"%d\" font-size=\"11\" text-anchor=\"middle\">%s</text>",
    sx(xticks), height-bottom, sx(xticks), height-bottom+5, sx(xticks), height-bottom+20,
    formatC(xticks, digits=2, format="fg")
  ), collapse = "")
  y_tick_svg <- paste0(sprintf(
    "<line x1=\"%d\" y1=\"%.2f\" x2=\"%d\" y2=\"%.2f\" stroke=\"#9aa0a6\"/><text x=\"%d\" y=\"%.2f\" font-size=\"11\" text-anchor=\"end\" dominant-baseline=\"middle\">%s</text>",
    left-5, sy(yticks), left, sy(yticks), left-9, sy(yticks), formatC(yticks, digits=2, format="fg")
  ), collapse = "")
  paste0(
    "<div class=\"plot\"><h3>", html_escape(title), "</h3><svg viewBox=\"0 0 ", width, " ", height, "\" role=\"img\">",
    "<line x1=\"", left, "\" y1=\"", height-bottom, "\" x2=\"", width-right, "\" y2=\"", height-bottom, "\" stroke=\"#3c4043\"/>",
    "<line x1=\"", left, "\" y1=\"", top, "\" x2=\"", left, "\" y2=\"", height-bottom, "\" stroke=\"#3c4043\"/>",
    x_tick_svg, y_tick_svg, ref, points,
    "<text x=\"", (left+width-right)/2, "\" y=\"", height-8, "\" text-anchor=\"middle\" font-size=\"12\">", html_escape(xlab), "</text>",
    "<text transform=\"translate(16 ", (top+height-bottom)/2, ") rotate(-90)\" text-anchor=\"middle\" font-size=\"12\">", html_escape(ylab), "</text>",
    "</svg></div>"
  )
}

reviewer_comment_map <- function(write_html = TRUE) {
  map <- data.frame(
    `Reviewer issue` = c(
      "RRT composition and 'CRRT-specific' framing",
      "HD-only subjects and interpretation of native clearance",
      "CRRT intensity, flow range, duration, and modality switching",
      "Simulation weight/flow assumption",
      "Limited CRRT-specific observations",
      "Two high ELF observations and influence on fit",
      "Piecewise versus conventional single-volume structure",
      "Fixed intermittent-HD clearance uncertainty",
      "Complete goodness-of-fit/residual plots",
      "Baseline characteristics by source cohort",
      "Post-filter scalar empirical support",
      "ELF/plasma relationship or approximate penetration ratio",
      "Observed concentrations above assay linearity range",
      "CLSI/EUCAST breakpoint annotations",
      "Laboratory/BAL procedural deviations and dilution integrity",
      "Continuous-infusion stability, title, abstract, prose, and grammar"
    ),
    `Code response` = c(
      "Subject and observation counts for CRRT-only, HD-only, both, and neither by development/validation cohort",
      "Reviewer run 103 excludes HD-only subjects; parameter medians are compared with the primary model",
      "Absolute flow, mL/kg/h intensity, observed active duration/span, and subjects with multiple CRRT modalities",
      "Reads actual simulation templates and reports absolute flow and mL/kg/h for each weight/regimen",
      "Counts observations and contributing subjects by matrix and RRT category",
      "Reviewer run 106 removes exactly the two highest development ELF observations and reports before/after ELF metrics",
      "Reviewer runs 101/102 replace Von/Voff with one V1; runs 107/108 add fixed WT/70 volume and CrCl/120 clearance scaling while retaining every other structural component",
      "Reviewer runs 104/105 bracket fixed CL_HD by -50%/+50% on the same development dataset",
      "Population/posterior observed-versus-predicted and Pmetrics weighted prediction error versus prediction/time for all five outputs",
      "Generated when optional de-identified source_cohort is retained in local input files",
      "Exact-time post-filter/pre-filter paired ratios and regression",
      "Nearest-time ELF/pre-filter pairs within the prespecified 0.5 h window; explicitly descriptive",
      "Reviewer runs 109/110 right-censor OUTEQ 1-3 values above 100 mg/L at the assay ULOQ and refit/validate the final model",
      "Not hard-coded: requires a prespecified agency, standard/version, organism context, and breakpoint before figure annotation",
      "Not answerable from model code; requires source laboratory, collection, and quality records",
      "Manuscript/literature response rather than a model-analysis task"
    ),
    `Output` = c(
      "reviewer-data-audit.html",
      "reviewer-model-sensitivity.html",
      "reviewer-data-audit.html",
      "reviewer-data-audit.html",
      "reviewer-data-audit.html",
      "reviewer-model-sensitivity.html",
      "reviewer-model-sensitivity.html",
      "reviewer-model-sensitivity.html",
      "reviewer-diagnostics.html",
      "reviewer-data-audit.html",
      "reviewer-data-audit.html",
      "reviewer-data-audit.html",
      "reviewer-data-audit.html",
      "Requires manuscript decision",
      "Requires source records",
      "Requires manuscript revision"
    ),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  if (write_html) {
    write_reviewer_html(
      REVIEWER_PATHS$comment_map_html,
      "Reviewer comment-to-analysis map",
      "What the public code answers, where the result appears, and what requires non-code evidence.",
      list(list(
        title = "Disposition map",
        intro = "This prevents analytical requests from being mixed with laboratory verification, literature interpretation, or editorial revision.",
        html = data_frame_html(map)
      )),
      notes = c(
        "Reviewer-only runs are sensitivity analyses and do not enter the primary model-selection sequence.",
        "A breakpoint should not be silently selected by code because the appropriate value depends on the named standard, version, organism, dosing context, and manuscript claim."
      )
    )
  }
  invisible(map)
}

# ---- reviewer data audit -----------------------------------------------------
modality_subject_table <- function(dev, val) {
  build <- function(data, label) {
    tab <- as.data.frame(table(subject_modality(data)$modality), stringsAsFactors = FALSE)
    names(tab) <- c("RRT category", "Subjects")
    tab$Dataset <- label
    tab[, c("Dataset", "RRT category", "Subjects")]
  }
  rbind(build(dev, "Development"), build(val, "Validation"))
}

modality_observation_table <- function(dev, val) {
  build <- function(data, label) {
    obs <- observation_rows(data)
    subjects <- subject_modality(data)[, c("id", "modality")]
    obs <- merge(obs, subjects, by = "id", all.x = TRUE)
    tab <- as.data.frame(table(obs$modality, matrix_label(obs$outeq)), stringsAsFactors = FALSE)
    names(tab) <- c("RRT category", "Matrix", "Observations")
    tab <- tab[tab$Observations > 0, , drop = FALSE]
    tab$Dataset <- label
    tab[, c("Dataset", "RRT category", "Matrix", "Observations")]
  }
  rbind(build(dev, "Development"), build(val, "Validation"))
}

baseline_table <- function(dev, val) {

  dev <- normalize_analysis_frame(dev)
  val <- normalize_analysis_frame(val)

  subject_table <- function(data) {
    ids <- unique(data$id)

    out <- lapply(ids, function(subject_id) {
      x <- data[data$id == subject_id, , drop = FALSE]
      x <- x[order(x$time), , drop = FALSE]
      first <- x[1L, , drop = FALSE]

      crrt_rows <- x[
        x$crrt == 1 & is.finite(x$flow),
        ,
        drop = FALSE
      ]

      initial_flow <- if (nrow(crrt_rows)) {
        crrt_rows$flow[[1L]] / 1000
      } else {
        NA_real_
      }

      data.frame(
        id = subject_id,
        age = first$age[[1L]],
        ht = first$ht[[1L]],
        wt = first$wt[[1L]],
        scr = first$scr[[1L]],
        bsa = first$bsa[[1L]],
        crcl = first$crcl[[1L]],
        initial_crrt_flow = initial_flow,
        ecmo = any(x$ecmo == 1, na.rm = TRUE),
        ever_hd = any(x$hd == 1, na.rm = TRUE),
        ever_crrt = any(x$crrt == 1, na.rm = TRUE),
        cvvh = any(x$cvvh == 1, na.rm = TRUE),
        cvvhd = any(x$cvvhd == 1, na.rm = TRUE),
        cvvhdf = any(x$cvvhdf == 1, na.rm = TRUE),
        male = first$male[[1L]],
        stringsAsFactors = FALSE
      )
    })

    do.call(rbind, out)
  }

  dev_s <- subject_table(dev)
  val_s <- subject_table(val)
  all_s <- rbind(dev_s, val_s)

  med_iqr <- function(x, digits = 1L) {
    x <- as_number(x)
    x <- x[is.finite(x)]

    if (!length(x)) return("Not available")

    q <- stats::quantile(
      x,
      probs = c(0.25, 0.50, 0.75),
      na.rm = TRUE,
      names = FALSE
    )

    sprintf(
      paste0(
        "%.", digits, "f (%.", digits, "f - %.", digits, "f)"
      ),
      q[[2L]], q[[1L]], q[[3L]]
    )
  }

  n_pct <- function(flag, denom) {
    flag <- as.logical(flag)
    n <- sum(flag, na.rm = TRUE)
    sprintf("%d (%.1f%%)", n, 100 * n / denom)
  }

  make_column <- function(x) {
    n <- nrow(x)

    c(
      med_iqr(x$age, 1L),
      med_iqr(x$ht, 1L),
      med_iqr(x$wt, 1L),
      med_iqr(x$scr, 1L),
      med_iqr(x$bsa, 1L),
      med_iqr(x$crcl, 1L),
      med_iqr(x$initial_crrt_flow, 1L),
      n_pct(x$ecmo, n),
      n_pct(x$ever_hd, n),
      n_pct(x$ever_crrt, n),
      n_pct(x$cvvh, n),
      n_pct(x$cvvhd, n),
      n_pct(x$cvvhdf, n),
      "",
      n_pct(x$male == 1, n),
      n_pct(x$male == 0, n)
    )
  }

  data.frame(
    Characteristic = c(
      "Age (years)",
      "Height (cm)",
      "Weight (kg)",
      "SCr (mg/dL)",
      "BSA (m2)",
      "CrCl (mL/min)",
      "Flow (L/h) at initial CRRT",
      "ECMO, n (%)",
      "Ever received HD, n (%)",
      "Ever received CRRT, n (%)",
      "CVVH, n (%)",
      "CVVHD, n (%)",
      "CVVHDF, n (%)",
      "Sex, n (%)",
      "Male",
      "Female"
    ),
    setNames(
      list(make_column(all_s)),
      paste0("Overall (N = ", nrow(all_s), ")")
    ),
    setNames(
      list(make_column(dev_s)),
      paste0("Development (N = ", nrow(dev_s), ")")
    ),
    setNames(
      list(make_column(val_s)),
      paste0("Validation (N = ", nrow(val_s), ")")
    ),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
}
source_cohort_table <- function(dev, val) {
  build <- function(data, dataset) {
    source_column <- find_source_cohort_column(data)
    if (is.na(source_column)) return(NULL)
    first <- first_subject_rows(data)
    tab <- as.data.frame(table(first[[source_column]], useNA = "ifany"), stringsAsFactors = FALSE)
    names(tab) <- c("Source cohort", "Subjects")
    tab$Dataset <- dataset
    tab[, c("Dataset", "Source cohort", "Subjects")]
  }
  output <- rbind(build(dev, "Development"), build(val, "Validation"))
  if (is.null(output) || nrow(output) == 0L) {
    data.frame(
      Status = "No de-identified source-cohort column found. Add optional `source_cohort` to local development/validation CSVs to generate this requested table.",
      stringsAsFactors = FALSE
    )
  } else output
}


observation_period_by_sampling_table <- function(dev, val) {
  build <- function(data, dataset_label) {
    data <- normalize_analysis_frame(data)
    by_id <- split(data, data$id)

    subject_level <- do.call(
      rbind,
      lapply(by_id, function(rows) {
        times <- rows$time[is.finite(rows$time)]

        span <- if (length(times) >= 2L) {
          diff(range(times))
        } else {
          0
        }

        doses <- sum(
          rows$evid == 1 |
            (is.finite(rows$dose) & rows$dose > 0),
          na.rm = TRUE
        )

        obs <- observation_rows(rows)

        data.frame(
          id = rows$id[[1L]],
          `Sampling group` = if ("sampling_group" %in% names(rows)) rows$sampling_group[[1L]] else NA_character_,
          `Observed follow-up, h` = span,
          `Dose records` = doses,
          `PK observations` = nrow(obs),
          check.names = FALSE,
          stringsAsFactors = FALSE
        )
      })
    )

    groups <- split(
      subject_level,
      subject_level[["Sampling group"]]
    )

    do.call(
      rbind,
      lapply(names(groups), function(group_name) {
        x <- groups[[group_name]]

        data.frame(
          Dataset = dataset_label,
          `Sampling group` = group_name,
          Subjects = nrow(x),
          `Observed follow-up, h` =
            summary_numbers(x[["Observed follow-up, h"]], 1L),
          `Dose records per subject` =
            summary_numbers(x[["Dose records"]], 0L),
          `PK observations per subject` =
            summary_numbers(x[["PK observations"]], 0L),
          check.names = FALSE,
          stringsAsFactors = FALSE
        )
      })
    )
  }

  rbind(
    build(dev, "Development"),
    build(val, "Validation")
  )
}


crrt_exposure_by_sampling_table <- function(dev, val) {
  build <- function(data, dataset_label) {
    data <- normalize_analysis_frame(data)

    duration <- observed_state_duration(data, "crrt")
    duration <- duration[
      duration$sessions > 0,
      ,
      drop = FALSE
    ]
    id_to_group <- if ("sampling_group" %in% names(data)) {
      tapply(data$sampling_group, data$id, function(x) x[!is.na(x)][1L])
    } else {
      stats::setNames(
        rep(NA_character_, length(unique(data$id))),
        as.character(unique(data$id))
      )
    }
    duration[["Sampling group"]] <- unname(id_to_group[as.character(duration$id)])

    unique_state <- unique(
      data[
        ,
        intersect(
          c("id", "time", "crrt", "flow", "wt"),
          names(data)
        ),
        drop = FALSE
      ]
    )

    active <- unique_state[
      unique_state$crrt == 1 &
        is.finite(unique_state$flow) &
        is.finite(unique_state$wt),
      ,
      drop = FALSE
    ]

    active[["Sampling group"]] <- unname(id_to_group[as.character(active$id)])
    active$intensity <- active$flow / active$wt

    groups <- unique(c(
      as.character(duration[["Sampling group"]]),
      as.character(active[["Sampling group"]])
    ))

    groups <- groups[!is.na(groups)]

    do.call(
      rbind,
      lapply(groups, function(group_name) {
        d <- duration[
          duration[["Sampling group"]] == group_name,
          ,
          drop = FALSE
        ]

        a <- active[
          active[["Sampling group"]] == group_name,
          ,
          drop = FALSE
        ]

        data.frame(
          Dataset = dataset_label,
          `Sampling group` = group_name,
          `Subjects ever CRRT` = length(unique(d$id)),
          `Absolute effluent flow, mL/h` =
            summary_numbers(a$flow, 1L),
          `Recorded effluent intensity, mL/kg/h` =
            summary_numbers(a$intensity, 1L),
          `Observed CRRT-active duration, h` =
            summary_numbers(d$active_h, 1L),
          `Observed CRRT span, h` =
            summary_numbers(d$span_h, 1L),
          check.names = FALSE,
          stringsAsFactors = FALSE
        )
      })
    )
  }

  rbind(
    build(dev, "Development"),
    build(val, "Validation")
  )
}

crrt_exposure_table <- function(dev, val) {
  build <- function(data, label) {
    data <- normalize_analysis_frame(data)
    unique_state <- unique(data[, intersect(c("id", "time", "crrt", "flow", "wt"), names(data)), drop = FALSE])
    active <- unique_state[unique_state$crrt == 1 & is.finite(unique_state$flow), , drop = FALSE]
    intensity <- active$flow / active$wt
    duration <- observed_state_duration(data, "crrt")
    duration <- duration[duration$sessions > 0, , drop = FALSE]
    data.frame(
      Dataset = label,
      `Subjects ever CRRT` = sum(subject_modality(data)$ever_crrt),
      `Absolute effluent flow, mL/h` = summary_numbers(active$flow, 1L),
      `Recorded effluent intensity, mL/kg/h` = summary_numbers(intensity, 1L),
      `Observed CRRT-active duration, h` = summary_numbers(duration$active_h, 1L),
      `Observed CRRT span, h` = summary_numbers(duration$span_h, 1L),
      check.names = FALSE,
      stringsAsFactors = FALSE
    )
  }
  rbind(build(dev, "Development"), build(val, "Validation"))
}

sample_distribution_table <- function(dev, val) {
  build <- function(data, label) {
    obs <- observation_rows(data)
    tab <- as.data.frame(table(matrix_label(obs$outeq)), stringsAsFactors = FALSE)
    names(tab) <- c("Matrix", "Observations")
    subjects <- vapply(split(obs$id, matrix_label(obs$outeq)), function(x) length(unique(x)), integer(1L))
    tab$Subjects <- subjects[match(tab$Matrix, names(subjects))]
    tab$Dataset <- label
    tab[, c("Dataset", "Matrix", "Observations", "Subjects")]
  }
  rbind(build(dev, "Development"), build(val, "Validation"))
}

observation_period_table <- function(dev, val) {
  build <- function(data, label) {
    data <- normalize_analysis_frame(data)
    by_id <- split(data, data$id)
    span <- vapply(by_id, function(rows) diff(range(rows$time, na.rm = TRUE)), numeric(1L))
    doses <- vapply(by_id, function(rows) sum(rows$evid == 1 | (is.finite(rows$dose) & rows$dose > 0), na.rm = TRUE), numeric(1L))
    data.frame(
      Dataset = label,
      Subjects = length(by_id),
      `Observed follow-up, h` = summary_numbers(span, 1L),
      `Dose records per subject` = summary_numbers(doses, 0L),
      check.names = FALSE,
      stringsAsFactors = FALSE
    )
  }
  rbind(build(dev, "Development"), build(val, "Validation"))
}

modality_switch_table <- function(dev, val) {
  build <- function(data, label) {
    data <- normalize_analysis_frame(data)
    ids <- unique(data$id)
    counts <- vapply(ids, function(id) {
      rows <- data[data$id == id, , drop = FALSE]
      sum(c(any(rows$cvvh == 1, na.rm=TRUE), any(rows$cvvhd == 1, na.rm=TRUE), any(rows$cvvhdf == 1, na.rm=TRUE)))
    }, integer(1L))
    data.frame(
      Dataset = label,
      `Subjects with >1 recorded CRRT modality` = sum(counts > 1L),
      `Maximum modalities in one subject` = if (length(counts)) max(counts) else 0L,
      check.names = FALSE,
      stringsAsFactors = FALSE
    )
  }
  rbind(build(dev, "Development"), build(val, "Validation"))
}

weight_flow_relationship_table <- function(dev, val) {
  build <- function(data, label) {
    data <- normalize_analysis_frame(data)
    rows <- unique(data[, c("id", "time", "wt", "flow", "crrt"), drop = FALSE])
    rows <- rows[rows$crrt == 1 & is.finite(rows$wt) & is.finite(rows$flow), , drop = FALSE]
    fit <- if (nrow(rows) >= 3L && stats::sd(rows$wt) > 0) stats::lm(flow ~ wt, data = rows) else NULL
    data.frame(
      Dataset = label,
      `CRRT-active time points` = nrow(rows),
      `Weight-flow correlation` = if (nrow(rows) >= 3L) format_number(stats::cor(rows$wt, rows$flow), 3L) else "",
      `Flow-on-weight R2` = if (is.null(fit)) "" else format_number(summary(fit)$r.squared, 3L),
      `Flow change per kg, mL/h` = if (is.null(fit)) "" else format_number(stats::coef(fit)[[2L]], 2L),
      check.names = FALSE,
      stringsAsFactors = FALSE
    )
  }
  rbind(build(dev, "Development"), build(val, "Validation"))
}

hd_event_alignment_table <- function(dev, val) {
  build <- function(data, label) {
    data <- normalize_analysis_frame(data)
    obs <- observation_rows(data)
    obs$HD_status <- ifelse(obs$hd == 1, "HD active", "HD inactive")
    obs_tab <- as.data.frame(table(obs$HD_status, matrix_label(obs$outeq)), stringsAsFactors = FALSE)
    names(obs_tab) <- c("HD status", "Matrix/event", "Records")
    doses <- data[(data$evid == 1 | (is.finite(data$dose) & data$dose > 0)), , drop = FALSE]
    if (nrow(doses)) {
      dose_tab <- as.data.frame(table(ifelse(doses$hd == 1, "HD active", "HD inactive")), stringsAsFactors = FALSE)
      names(dose_tab) <- c("HD status", "Records")
      dose_tab$`Matrix/event` <- "Dose record"
      dose_tab <- dose_tab[, c("HD status", "Matrix/event", "Records")]
      obs_tab <- rbind(obs_tab, dose_tab)
    }
    obs_tab <- obs_tab[obs_tab$Records > 0, , drop = FALSE]
    obs_tab$Dataset <- label
    obs_tab[, c("Dataset", "HD status", "Matrix/event", "Records")]
  }
  rbind(build(dev, "Development"), build(val, "Validation"))
}

assay_range_table <- function(dev, val) {
  build <- function(data, label) {
    obs <- observation_rows(data)

    # Output equation 4 is cumulative effluent amount, not a concentration,
    # and therefore cannot be compared with the 100 mg/L assay range.
    obs <- obs[obs$outeq %in% 1:3, , drop = FALSE]

    tab <- as.data.frame(
      table(matrix_label(obs$outeq), obs$out > 100),
      stringsAsFactors = FALSE
    )
    names(tab) <- c("Matrix", "Above 100 mg/L", "Observations")
    tab <- tab[
      tab$`Above 100 mg/L` == "TRUE" & tab$Observations > 0,
      ,
      drop = FALSE
    ]
    tab$Dataset <- label
    tab[, c("Dataset", "Matrix", "Observations")]
  }

  output <- rbind(build(dev, "Development"), build(val, "Validation"))

  if (nrow(output) == 0L) {
    data.frame(Status = "No observed concentrations above 100 mg/L.")
  } else {
    output
  }
}

paired_matrix_values <- function(data, target_eq, reference_eq, window_h = 0) {
  obs <- observation_rows(data)
  target <- obs[obs$outeq == target_eq, c("id", "time", "out"), drop = FALSE]
  reference <- obs[obs$outeq == reference_eq, c("id", "time", "out"), drop = FALSE]
  names(target)[[3L]] <- "target"
  names(reference)[[3L]] <- "reference"
  pairs <- list(); counter <- 0L
  for (i in seq_len(nrow(target))) {
    candidates <- reference[reference$id == target$id[[i]], , drop = FALSE]
    if (nrow(candidates) == 0L) next
    delta <- abs(candidates$time - target$time[[i]])
    closest <- which.min(delta)
    if (delta[[closest]] <= window_h) {
      counter <- counter + 1L
      pairs[[counter]] <- data.frame(
        target = target$target[[i]],
        reference = candidates$reference[[closest]],
        delta_h = delta[[closest]]
      )
    }
  }
  if (counter == 0L) return(data.frame(target=numeric(), reference=numeric(), delta_h=numeric()))
  do.call(rbind, pairs)
}

paired_summary_table <- function(dev, val) {
  build <- function(data, label, target_eq, reference_eq, comparison, window) {
    pairs <- paired_matrix_values(data, target_eq, reference_eq, window)
    ratio <- pairs$target / pairs$reference
    fit <- if (nrow(pairs) >= 3L) stats::lm(target ~ reference, data = pairs) else NULL
    data.frame(
      Dataset = label,
      Comparison = comparison,
      Pairs = nrow(pairs),
      `Maximum time difference, h` = if (nrow(pairs)) format_number(max(pairs$delta_h), 2L) else "",
      `Target/reference ratio` = summary_numbers(ratio, 2L),
      `R2` = if (is.null(fit)) "" else format_number(summary(fit)$r.squared, 3L),
      Intercept = if (is.null(fit)) "" else format_number(stats::coef(fit)[[1L]], 3L),
      Slope = if (is.null(fit)) "" else format_number(stats::coef(fit)[[2L]], 3L),
      check.names = FALSE,
      stringsAsFactors = FALSE
    )
  }
  rbind(
    build(dev, "Development", 2, 1, "Post-filter / pre-filter plasma", 0),
    build(val, "Validation", 2, 1, "Post-filter / pre-filter plasma", 0),
    build(dev, "Development", 5, 1, "ELF / pre-filter plasma", REVIEWER_ELF_PAIR_WINDOW_H),
    build(val, "Validation", 5, 1, "ELF / pre-filter plasma", REVIEWER_ELF_PAIR_WINDOW_H)
  )
}

simulation_assumption_table <- function() {
  templates <- load_simulation_templates()
  output <- list(); counter <- 0L
  for (code in names(templates)) {
    template_path <- templates[[code]]

    data <- utils::read.csv(
      template_path,
      check.names = FALSE,
      stringsAsFactors = FALSE,
      na.strings = c("", "NA", ".", "NaN")
    )
    data <- normalize_analysis_frame(data)

    rows <- unique(
      data[, c("id", "time", "wt", "flow", "crrt"), drop = FALSE]
    )
    ids <- unique(rows$id)
    for (id in ids) {
      one <- rows[rows$id == id & rows$crrt == 1 & is.finite(rows$flow), , drop = FALSE]
      if (nrow(one) == 0L) next
      counter <- counter + 1L
      output[[counter]] <- data.frame(
        Regimen = code,
        `Weight (kg)` = unique(one$wt[is.finite(one$wt)])[[1L]],
        `Absolute effluent flow (mL/h)` = summary_numbers(one$flow, 1L),
        `Effluent intensity (mL/kg/h)` = summary_numbers(one$flow / one$wt, 1L),
        check.names = FALSE,
        stringsAsFactors = FALSE
      )
    }
  }
  if (counter == 0L) data.frame(Status="No CRRT-active simulation rows found.") else do.call(rbind, output)
}

reviewer_data_audit <- function(write_html = TRUE) {
  dev <- read_local_analysis_frame("development")
  val <- read_local_analysis_frame("validation")
  sections <- list(
    list(title="RRT composition by analysis cohort", intro="Subjects are classified from time-varying CRRT and HD indicators across their complete local record.", html=data_frame_html(modality_subject_table(dev, val))),
    list(title="Observations by subject-level RRT category", intro="This separates CRRT-only, HD-only, and crossover subjects and shows which matrices they contributed.", html=data_frame_html(modality_observation_table(dev, val))),
    list(title="Sample distribution by matrix", intro="Counts of observations and contributing subjects in development and validation.", html=data_frame_html(sample_distribution_table(dev, val))),
    list(title="Table 1. Baseline Characteristics and Clinical Characteristics", intro="Continuous variables are presented as median (interquartile range), and categorical variables as n (%).", html=data_frame_html(baseline_table(dev, val))),
    list(title="Original sampling cohort", intro="Generated only when an optional de-identified source-cohort field is retained locally.", html=data_frame_html(source_cohort_table(dev, val))),
    list(title="Observed treatment/follow-up period by sampling group", intro="Follow-up, dose records, and PK observations are summarized separately for prospective richly sampled and sparse/opportunistically sampled subjects.", html=data_frame_html(observation_period_by_sampling_table(dev, val))),
    list(title="CRRT exposure and intensity by sampling group", intro="CRRT exposure is summarized separately for prospective richly sampled and sparse/opportunistically sampled subjects. Recorded effluent intensity is absolute flow divided by contemporaneous weight.", html=data_frame_html(crrt_exposure_by_sampling_table(dev, val))),
    list(title="Weight representation within recorded effluent flow", intro="This quantifies the extent to which weight is already embedded in the recorded CRRT flow term.", html=data_frame_html(weight_flow_relationship_table(dev, val))),
    list(title="Dose and observation alignment with HD status", intro="Counts show which dose and biological-matrix records were assigned to intradialytic versus interdialytic periods by the time-varying HD indicator.", html=data_frame_html(hd_event_alignment_table(dev, val))),
    list(title="CRRT modality switching", intro="Counts greater than the number ever receiving CRRT can occur when a subject contributes more than one modality over follow-up.", html=data_frame_html(modality_switch_table(dev, val))),
    list(title="Observed values above the stated 100 mg/L assay range", intro="This quantifies records requiring a laboratory dilution-integrity explanation; code cannot verify the laboratory procedure itself.", html=data_frame_html(assay_range_table(dev, val))),
    list(title="Empirical paired concentration summaries", intro=paste0("Post-filter pairs require exact matching time. ELF/plasma pairs use the nearest pre-filter sample within ", REVIEWER_ELF_PAIR_WINDOW_H, " h. These ratios are descriptive and do not replace dynamic modeling."), html=data_frame_html(paired_summary_table(dev, val))),
    list(title="Simulation body-weight and CRRT-flow assumptions", intro="This reads the actual local simulation templates rather than inferring the assumption from manuscript text.", html=data_frame_html(simulation_assumption_table()))
  )
  if (write_html) {
    write_reviewer_html(
      REVIEWER_PATHS$audit_html,
      "Reviewer-response data and design audit",
      "Cohort composition, CRRT/HD exposure, sampling structure, paired matrices, and simulation assumptions.",
      sections,
      notes = c(
        "Observed active duration is bounded by the first and last recorded dataset times and is not a substitute for a separately curated clinical RRT-duration variable.",
        "No subject identifiers are written to this report.",
        "Procedural issues such as dilution integrity, BAL collection deviations, and laboratory investigation require source documentation and cannot be resolved by the model code alone."
      )
    )
  }
  invisible(sections)
}

# ---- reviewer model fits -----------------------------------------------------
remove_hd_only_subjects <- function(data) {
  modality <- subject_modality(data)
  exclude <- modality$id[modality$modality == "HD only"]
  data[!data$id %in% exclude, , drop = FALSE]
}

remove_highest_elf_observations <- function(data, n = REVIEWER_ELF_EXCLUDE_N) {
  normalized <- normalize_analysis_frame(data)
  candidates <- which(normalized$outeq == 5 & is.finite(normalized$out))
  if (length(candidates) < n) stop("Fewer than ", n, " finite ELF observations are available.", call. = FALSE)
  ranked <- candidates[order(normalized$out[candidates], decreasing = TRUE)]
  data[-head(ranked, n), , drop = FALSE]
}

censor_above_assay_uloq <- function(data, uloq = REVIEWER_ASSAY_ULOQ) {
  normalized <- normalize_analysis_frame(data)
  if (!"cens" %in% names(data)) data$cens <- NA_integer_
  data$cens <- suppressWarnings(as.integer(as.character(data$cens)))

  observed <- !is.na(normalized$evid) & normalized$evid == 0 & is.finite(normalized$out)
  above <- observed & normalized$outeq %in% 1:3 & normalized$out > uloq

  # Pmetrics 3.x standard data use CENS = -1 ("aloq") for above-limit observations.
  # Right-censored rows carry the censoring boundary in OUT; the unvalidated numeric
  # value above the assay range is not treated as an exact concentration.
  data$cens[observed] <- 0L
  data$cens[above] <- -1L
  data$out[above] <- uloq
  data
}

fit_reviewer_run_101 <- function(overwrite = FALSE) {
  data_run_101 <- read_local_analysis_frame("development")
  pm_data_run_101 <- PM_data$new(data_run_101, loq = rep(0, 5))
  model_run_101 <- make_single_volume_model()
  message("Fitting reviewer run 101: single-volume development; cycles=1000, points = 300, seed=12345")
  model_run_101$fit(
    data = pm_data_run_101, cycles = 1000, path = PATHS$runs, run = 101,
    points = 300, seed = 12345, overwrite = overwrite, report = "plotly"
  )
}

fit_reviewer_run_103 <- function(overwrite = FALSE) {
  data_run_103 <- remove_hd_only_subjects(read_local_analysis_frame("development"))
  pm_data_run_103 <- PM_data$new(data_run_103, loq = rep(0, 5))
  model_run_103 <- make_base_model()
  message("Fitting reviewer run 103: exclude HD-only subjects; cycles=1000, points = 300, seed=12345")
  model_run_103$fit(
    data = pm_data_run_103, cycles = 1000, path = PATHS$runs, run = 103,
    points = 300, seed = 12345, overwrite = overwrite, report = "plotly"
  )
}

fit_reviewer_run_104 <- function(overwrite = FALSE) {
  data_run_104 <- read_local_analysis_frame("development")
  pm_data_run_104 <- PM_data$new(data_run_104, loq = rep(0, 5))
  model_run_104 <- make_base_hd_low_model()
  message("Fitting reviewer run 104: fixed CL_HD=3.6 L/h; cycles=1000, points = 300, seed=12345")
  model_run_104$fit(
    data = pm_data_run_104, cycles = 1000, path = PATHS$runs, run = 104,
    points = 300, seed = 12345, overwrite = overwrite, report = "plotly"
  )
}

fit_reviewer_run_105 <- function(overwrite = FALSE) {
  data_run_105 <- read_local_analysis_frame("development")
  pm_data_run_105 <- PM_data$new(data_run_105, loq = rep(0, 5))
  model_run_105 <- make_base_hd_high_model()
  message("Fitting reviewer run 105: fixed CL_HD=10.8 L/h; cycles=1000, points = 300, seed=12345")
  model_run_105$fit(
    data = pm_data_run_105, cycles = 1000, path = PATHS$runs, run = 105,
    points = 300, seed = 12345, overwrite = overwrite, report = "plotly"
  )
}

fit_reviewer_run_106 <- function(overwrite = FALSE) {
  data_run_106 <- remove_highest_elf_observations(read_local_analysis_frame("development"))
  pm_data_run_106 <- PM_data$new(data_run_106, loq = rep(0, 5))
  model_run_106 <- make_base_model()
  message("Fitting reviewer run 106: remove two highest ELF observations; cycles=1000, points = 300, seed=12345")
  model_run_106$fit(
    data = pm_data_run_106, cycles = 1000, path = PATHS$runs, run = 106,
    points = 300, seed = 12345, overwrite = overwrite, report = "plotly"
  )
}

fit_reviewer_run_107 <- function(overwrite = FALSE) {
  data_run_107 <- read_local_analysis_frame("development")
  pm_data_run_107 <- PM_data$new(data_run_107, loq = rep(0, 5))
  model_run_107 <- make_conventional_single_volume_model()
  message("Fitting reviewer run 107: conventional single-volume development; cycles=1000, points = 300, seed=12345")
  model_run_107$fit(
    data = pm_data_run_107, cycles = 1000, path = PATHS$runs, run = 107,
    points = 300, seed = 12345, overwrite = overwrite, report = "plotly"
  )
}

reviewer_validate_single_volume <- function(overwrite = FALSE) {
  if (!reviewer_run_exists(101L)) {
    stop(
      "Reviewer run 101 is missing. Fit it before validation.",
      call. = FALSE
    )
  }

  validation_data <- load_analysis_data("validation")
  model <- make_single_volume_model()

  message(
    "Fitting reviewer run 102 held-out validation with cycles=0 ",
    "using reviewer run 101 as the fixed prior."
  )

  model$fit(
    data = validation_data,
    cycles = 0,
    path = PATHS$runs,
    run = 102L,
    prior = 101L,
    overwrite = overwrite,
    report = "plotly"
  )
}


reviewer_validate_conventional_volume <- function(overwrite = FALSE) {
  if (!reviewer_run_exists(107L)) {
    stop(
      "Reviewer run 107 is missing. Fit it before validation.",
      call. = FALSE
    )
  }

  validation_data <- load_analysis_data("validation")
  model <- make_conventional_single_volume_model()

  message(
    "Fitting reviewer run 108 held-out validation with cycles=0 ",
    "using reviewer run 107 as the fixed prior."
  )

  model$fit(
    data = validation_data,
    cycles = 0,
    path = PATHS$runs,
    run = 108L,
    prior = 107L,
    overwrite = overwrite,
    report = "plotly"
  )
}

fit_reviewer_run_109 <- function(overwrite = FALSE) {
  data_run_109 <- censor_above_assay_uloq(read_local_analysis_frame("development"))
  pm_data_run_109 <- PM_data$new(data_run_109)
  model_run_109 <- make_base_model()
  message("Fitting reviewer run 109: OUTEQ 1-3 values >100 mg/L right-censored at 100 mg/L; cycles=1000, points = 300, seed=12345")
  model_run_109$fit(
    data = pm_data_run_109, cycles = 1000, path = PATHS$runs, run = 109L,
    points = 300, seed = 12345, overwrite = overwrite, report = "plotly"
  )
}

reviewer_validate_uloq_censored <- function(overwrite = FALSE) {
  if (!reviewer_run_exists(109L)) stop("Reviewer run 109 is missing. Fit it before validation.", call. = FALSE)
  validation_frame <- censor_above_assay_uloq(read_local_analysis_frame("validation"))
  validation_data <- PM_data$new(validation_frame)
  model <- make_base_model()
  message("Fitting reviewer run 110 held-out validation with the same ULOQ censoring rule and run 109 as the fixed prior.")
  model$fit(
    data = validation_data, cycles = 0, path = PATHS$runs, run = 110L, prior = 109L,
    overwrite = overwrite, report = "plotly"
  )
}

reviewer_fit <- function(keys = c("single_volume", "conventional_volume", "no_hd_only", "hd_low", "hd_high", "elf_top2_removed", "uloq_censored"), overwrite = FALSE) {
  valid <- c("single_volume", "conventional_volume", "no_hd_only", "hd_low", "hd_high", "elf_top2_removed", "uloq_censored")
  invalid <- setdiff(keys, valid)
  if (length(invalid)) stop("Unknown reviewer key(s): ", paste(invalid, collapse = ", "), call. = FALSE)
  for (key in keys) {
    if (identical(key, "single_volume")) {
      fit_reviewer_run_101(overwrite = overwrite); reviewer_validate_single_volume(overwrite = overwrite)
    } else if (identical(key, "conventional_volume")) {
      fit_reviewer_run_107(overwrite = overwrite); reviewer_validate_conventional_volume(overwrite = overwrite)
    } else if (identical(key, "no_hd_only")) {
      fit_reviewer_run_103(overwrite = overwrite)
    } else if (identical(key, "hd_low")) {
      fit_reviewer_run_104(overwrite = overwrite)
    } else if (identical(key, "hd_high")) {
      fit_reviewer_run_105(overwrite = overwrite)
    } else if (identical(key, "elf_top2_removed")) {
      fit_reviewer_run_106(overwrite = overwrite)
    } else if (identical(key, "uloq_censored")) {
      fit_reviewer_run_109(overwrite = overwrite); reviewer_validate_uloq_censored(overwrite = overwrite)
    }
  }
  invisible(TRUE)
}

# ---- diagnostics and sensitivity summaries ----------------------------------
prepare_op <- function(run, outeq, pred_type = c("post", "pop")) {
  pred_type <- match.arg(pred_type)
  data <- as.data.frame(run$op$data, stringsAsFactors = FALSE)
  names(data) <- tolower(names(data))
  required <- c("time", "obs", "pred", "pred.type", "icen", "outeq", "wd")
  missing <- setdiff(required, names(data))
  if (length(missing)) stop("Pmetrics OP data are missing: ", paste(missing, collapse=", "), call. = FALSE)
  data <- data[
    data$pred.type == pred_type & data$icen == "median" & data$outeq == outeq &
      is.finite(as_number(data$obs)) & is.finite(as_number(data$pred)),
    , drop = FALSE
  ]
  data$time <- as_number(data$time)
  data$obs <- as_number(data$obs)
  data$pred <- as_number(data$pred)
  data$wd <- as_number(data$wd)
  data
}

op_metrics <- function(run, run_label, cohort, outeq, pred_type) {
  data <- prepare_op(run, outeq, pred_type)
  fit <- if (nrow(data) >= 3L && stats::sd(data$pred) > 0) stats::lm(obs ~ pred, data = data) else NULL
  data.frame(
    Run = run_label,
    Cohort = cohort,
    Matrix = matrix_label(outeq),
    Prediction = if (pred_type == "pop") "Population" else "Posterior individual",
    N = nrow(data),
    R2 = if (is.null(fit)) NA_real_ else summary(fit)$r.squared,
    Intercept = if (is.null(fit)) NA_real_ else stats::coef(fit)[[1L]],
    Slope = if (is.null(fit)) NA_real_ else stats::coef(fit)[[2L]],
    Bias = if (nrow(data)) mean(data$pred - data$obs, na.rm = TRUE) else NA_real_,
    RMSE = if (nrow(data)) sqrt(mean((data$pred - data$obs)^2, na.rm = TRUE)) else NA_real_,
    `Mean weighted prediction error` = if (nrow(data)) mean(data$wd, na.rm = TRUE) else NA_real_,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
}

format_metrics_table <- function(data) {
  numeric_columns <- c("R2", "Intercept", "Slope", "Bias", "RMSE", "Mean weighted prediction error")
  for (column in numeric_columns) data[[column]] <- format_number(data[[column]], 3L)
  data
}

reviewer_diagnostics <- function(write_html = TRUE) {
  runs <- list(
    list(run=load_public_run(1L), label="Public run 1", cohort="Development"),
    list(run=load_public_run(5L), label="Public run 5", cohort="Validation")
  )
  metrics <- list(); metric_counter <- 0L; sections <- list(); section_counter <- 0L
  for (entry in runs) {
    for (outeq in 1:5) {
      pop <- prepare_op(entry$run, outeq, "pop")
      post <- prepare_op(entry$run, outeq, "post")
      metric_counter <- metric_counter + 1L
      metrics[[metric_counter]] <- op_metrics(entry$run, entry$label, entry$cohort, outeq, "pop")
      metric_counter <- metric_counter + 1L
      metrics[[metric_counter]] <- op_metrics(entry$run, entry$label, entry$cohort, outeq, "post")
      plots <- paste0(
        "<div class=\"plots\">",
        scatter_svg(pop$pred, pop$obs, "Observed vs population prediction", "Population prediction", "Observed", "identity"),
        scatter_svg(post$pred, post$obs, "Observed vs posterior prediction", "Posterior prediction", "Observed", "identity"),
        scatter_svg(post$pred, post$wd, "Individual weighted residual vs prediction", "Posterior prediction", "Weighted prediction error", "zero"),
        scatter_svg(post$time, post$wd, "Individual weighted residual vs time", "Time", "Weighted prediction error", "zero"),
        "</div>"
      )
      section_counter <- section_counter + 1L
      sections[[section_counter]] <- list(
        title = paste(entry$cohort, matrix_label(outeq)),
        intro = paste0(entry$label, "; output equation ", outeq, "."),
        html = plots
      )
    }
  }
  metrics <- do.call(rbind, metrics)
  sections <- c(list(list(
    title="Regression and error metrics",
    intro="Metrics are reported separately for population and posterior individual predictions.",
    html=data_frame_html(format_metrics_table(metrics))
  )), sections)
  if (write_html) {
    write_reviewer_html(
      REVIEWER_PATHS$diagnostics_html,
      "Reviewer-response goodness-of-fit diagnostics",
      "All five outputs in development and held-out validation.",
      sections,
      notes = c(
        "Pmetrics PM_op defines wd as (prediction - observation) divided by the observation SD. These are Pmetrics population or individual weighted prediction errors, not NONMEM CWRES terminology.",
        "The Pmetrics-generated Plotly reports remain available inside each local run directory for interactive inspection.",
        "Sparse effluent outputs should be interpreted with their observation counts."
      )
    )
  }
  invisible(metrics)
}

extract_population_medians <- function(run, label) {
  med <- as.data.frame(run$final$popMed, check.names = FALSE)
  names(med) <- tolower(names(med))
  values <- as.numeric(med[1L, , drop = TRUE])
  data.frame(Run=label, Parameter=names(med), Median=values, stringsAsFactors=FALSE)
}

same_data_model_comparison <- function() {
  run_numbers <- c(1L, 101L, 107L, 104L, 105L)

  runs <- list(
    load_public_run(1L),
    load_reviewer_run(101L),
    load_reviewer_run(107L),
    load_reviewer_run(104L),
    load_reviewer_run(105L)
  )

  labels <- c(
    `1` = "Primary piecewise V1",
    `101` = "Single V1; no fixed WT/CrCl scaling",
    `107` = "Single V1; fixed WT/CrCl scaling",
    `104` = paste0(
      "Piecewise V1; CL_HD = ",
      REVIEWER_HD_LOW,
      " L/h"
    ),
    `105` = paste0(
      "Piecewise V1; CL_HD = ",
      REVIEWER_HD_HIGH,
      " L/h"
    )
  )

  metrics <- do.call(
    rbind,
    Map(
      extract_final_cycle_metrics,
      runs,
      run_numbers
    )
  )

  output <- data.frame(
    Run = metrics$public_run,
    Model = unname(labels[as.character(metrics$public_run)]),
    Parameters = metrics$actual_parameters,
    `Final cycle` = metrics$actual_cycles,
    `-2LL` = metrics$actual_minus2ll,
    AIC = metrics$actual_aic,
    BIC = metrics$actual_bic,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  numeric_columns <- c("-2LL", "AIC", "BIC")

  for (column in numeric_columns) {
    output[[column]] <- format_number(output[[column]], 3L)
  }

  output
}


validation_structure_comparison <- function() {
  entries <- list(
    list(
      run = load_public_run(5L),
      label = "Run 5 - Primary piecewise V1 validation"
    ),
    list(
      run = load_reviewer_run(102L),
      label = "Run 102 - Single V1 validation"
    ),
    list(
      run = load_reviewer_run(108L),
      label = "Run 108 - Single V1 + fixed WT/CrCl scaling validation"
    )
  )

  output <- list()
  counter <- 0L

  for (outeq in 1:5) {
    for (entry in entries) {
      counter <- counter + 1L

      output[[counter]] <- op_metrics(
        entry$run,
        entry$label,
        "Validation",
        outeq,
        "post"
      )
    }
  }

  format_metrics_table(do.call(rbind, output))
}

sensitivity_parameter_table <- function() {
  runs <- list(
    list(
      run = load_public_run(1L),
      label = "Run 1 - Primary piecewise V1"
    ),
    list(
      run = load_reviewer_run(103L),
      label = "Run 103 - Exclude HD-only subjects"
    ),
    list(
      run = load_reviewer_run(104L),
      label = paste0(
        "Run 104 - Piecewise V1; CL_HD = ",
        REVIEWER_HD_LOW,
        " L/h"
      )
    ),
    list(
      run = load_reviewer_run(105L),
      label = paste0(
        "Run 105 - Piecewise V1; CL_HD = ",
        REVIEWER_HD_HIGH,
        " L/h"
      )
    ),
    list(
      run = load_reviewer_run(106L),
      label = "Run 106 - Exclude two highest ELF observations"
    )
    ,
    list(
      run = load_reviewer_run(109L),
      label = "Run 109 - OUTEQ 1-3 >100 mg/L right-censored at ULOQ"
    )
  )

  data <- do.call(
    rbind,
    lapply(
      runs,
      function(entry) {
        extract_population_medians(
          entry$run,
          entry$label
        )
      }
    )
  )

  wide <- reshape(
    data,
    idvar = "Parameter",
    timevar = "Run",
    direction = "wide"
  )

  names(wide) <- sub("^Median\\.", "", names(wide))

  for (name in names(wide)[-1L]) {
    wide[[name]] <- format_number(wide[[name]], 3L)
  }

  wide
}

elf_outlier_table <- function() {
  run <- load_public_run(1L)
  data <- prepare_op(run, 5L, "post")
  data <- data[order(data$obs, decreasing=TRUE), , drop=FALSE]
  data <- head(data, REVIEWER_ELF_EXCLUDE_N)
  data.frame(
    Rank = seq_len(nrow(data)),
    `Observed ELF` = format_number(data$obs, 3L),
    `Posterior prediction` = format_number(data$pred, 3L),
    `Weighted prediction error` = format_number(data$wd, 3L),
    check.names=FALSE,
    stringsAsFactors=FALSE
  )
}

elf_sensitivity_metrics <- function() {
  data <- rbind(
    op_metrics(
      load_public_run(1L),
      "Run 1 - Primary piecewise V1",
      "Development",
      5L,
      "post"
    ),
    op_metrics(
      load_reviewer_run(106L),
      "Run 106 - Exclude two highest ELF observations",
      "Development",
      5L,
      "post"
    )
  )

  format_metrics_table(data)
}

uloq_censoring_metrics <- function() {
  output <- list(); counter <- 0L
  entries <- list(
    list(run = load_public_run(5L), label = "Run 5 - Primary validation"),
    list(run = load_reviewer_run(110L), label = "Run 110 - ULOQ-censored validation")
  )
  for (outeq in 1:3) {
    for (entry in entries) {
      counter <- counter + 1L
      output[[counter]] <- op_metrics(entry$run, entry$label, "Validation", outeq, "post")
    }
  }
  format_metrics_table(do.call(rbind, output))
}

reviewer_model_sensitivity <- function(write_html = TRUE) {
  required <- REVIEWER_RUN_REGISTRY$reviewer_run
  missing <- required[!vapply(required, reviewer_run_exists, logical(1L))]
  if (length(missing)) stop("Missing reviewer runs: ", paste(missing, collapse=", "), ". Run `reviewer fit` first.", call. = FALSE)
  sections <- list(
    list(title="Same-data structural, conventional-covariate, and HD-clearance comparison", intro="Public run 1, the 9-parameter single-volume comparators, and the low/high fixed-HD-clearance fits use the same development dataset. Run 101 isolates the volume-state decision; run 107 adds fixed WT/70 volume and CrCl/120 native-clearance scaling.", html=data_frame_html(same_data_model_comparison())),
    list(title="Held-out predictive comparison: piecewise versus single central volume", intro="Both models are applied to the same validation cohort by MAP without relocating their population support points.", html=data_frame_html(validation_structure_comparison())),
    list(title="Population median sensitivity", intro="AIC should not be compared across the HD-only exclusion or ELF-observation exclusion because those analyses change the fitted data. Parameter shifts are shown instead.", html=data_frame_html(sensitivity_parameter_table())),
    list(title="Audit of the two highest development ELF observations", intro="Values are shown without subject identifiers. Procedural verification must be performed against laboratory and collection records.", html=data_frame_html(elf_outlier_table())),
    list(title="ELF fit before and after exclusion", intro="This quantifies whether the two highest observations dominate the development ELF regression.", html=data_frame_html(elf_sensitivity_metrics())),
    list(title="ULOQ right-censoring sensitivity", intro="OUTEQ 1-3 observations above 100 mg/L are encoded as Pmetrics above-limit observations (CENS=-1/ALOQ) with OUT set to the 100 mg/L censoring boundary. Run 109 refits development and run 110 applies that prior to held-out validation using the same rule.", html=data_frame_html(uloq_censoring_metrics()))
  )
  if (write_html) {
    write_reviewer_html(
      REVIEWER_PATHS$sensitivity_html,
      "Reviewer-response model sensitivity analyses",
      "Structural, HD-only, fixed-HD-clearance, ELF high-observation, and assay-ULOQ censoring analyses.",
      sections,
      notes=c(
        paste0("The HD-clearance sensitivity values ", REVIEWER_HD_LOW, " and ", REVIEWER_HD_HIGH, " L/h are transparent +/-50% brackets around the manuscript value of ", REVIEWER_HD_REFERENCE, " L/h; they are not claimed as empirical confidence limits."),
        "Run 101 is the direct isolated test of the piecewise central-volume decision. Run 107 is the conventional fixed-covariate comparator. Neither removes K12/K21.",
        "Excluding observations is a reviewer-requested influence analysis and does not redefine the primary dataset."
      )
    )
  }
  invisible(sections)
}

write_reviewer_index <- function() {
  links <- data.frame(
    Report=c("Comment-to-analysis map", "Data and design audit", "Model sensitivity", "Goodness-of-fit diagnostics"),
    File=c(
      basename(REVIEWER_PATHS$comment_map_html),
      basename(REVIEWER_PATHS$audit_html),
      basename(REVIEWER_PATHS$sensitivity_html),
      basename(REVIEWER_PATHS$diagnostics_html)
    ),
    stringsAsFactors=FALSE
  )
  links <- links[file.exists(file.path(PATHS$runs, links$File)), , drop=FALSE]
  rows <- paste0("<li><a href=\"", html_escape(links$File), "\">", html_escape(links$Report), "</a></li>", collapse="")
  html <- paste0(
    "<!doctype html><html><head><meta charset=\"utf-8\"><title>Reviewer response analyses</title>",
    "<style>body{font-family:Arial,Helvetica,sans-serif;max-width:850px;margin:3rem auto;padding:0 1rem;line-height:1.5}li{margin:.7rem 0}</style></head><body>",
    "<h1>Reviewer-response analyses</h1><p>The primary manuscript path remains public runs 1-5. Reviewer-only sensitivity runs are 101-110.</p><ul>", rows, "</ul></body></html>"
  )
  writeLines(html, REVIEWER_PATHS$index_html, useBytes=TRUE)
  message("Wrote ", REVIEWER_PATHS$index_html)
  invisible(REVIEWER_PATHS$index_html)
}

reviewer_report <- function() {
  reviewer_comment_map(TRUE)
  reviewer_data_audit(TRUE)
  reviewer_model_sensitivity(TRUE)
  reviewer_diagnostics(TRUE)
  write_reviewer_index()
  invisible(TRUE)
}

reviewer_usage <- function() {
  cat(
    "Reviewer-response analyses\n\n",
    "Usage:\n",
    "  Rscript Pmetrics/Rscript/Analysis.R reviewer audit\n",
    "  Rscript Pmetrics/Rscript/Analysis.R reviewer fit [key] [--overwrite]\n",
    "  Rscript Pmetrics/Rscript/Analysis.R reviewer report\n",
    "  Rscript Pmetrics/Rscript/Analysis.R reviewer all [--overwrite]\n\n",
    "Keys: single_volume, conventional_volume, no_hd_only, hd_low, hd_high, elf_top2_removed, uloq_censored\n",
    sep=""
  )
}

reviewer_main <- function(arguments) {
  if (length(arguments) == 0L || arguments[[1L]] %in% c("help", "-h", "--help")) {
    reviewer_usage(); return(invisible(TRUE))
  }
  action <- arguments[[1L]]
  overwrite <- "--overwrite" %in% arguments
  positional <- arguments[-1L]
  positional <- positional[positional != "--overwrite"]

  if (action == "audit") {
    run_check(TRUE, TRUE, TRUE)
    reviewer_comment_map(TRUE)
    reviewer_data_audit(TRUE)
    write_reviewer_index()
  } else if (action == "fit") {
    run_check(TRUE, TRUE, FALSE)
    keys <- if (length(positional)) positional else c("single_volume", "conventional_volume", "no_hd_only", "hd_low", "hd_high", "elf_top2_removed", "uloq_censored")
    reviewer_fit(keys, overwrite)
  } else if (action == "report") {
    run_check(TRUE, TRUE, TRUE)
    reviewer_report()
  } else if (action == "all") {
    run_check(TRUE, TRUE, TRUE)
    reviewer_comment_map(TRUE)
    reviewer_data_audit(TRUE)
    reviewer_fit(c("single_volume", "conventional_volume", "no_hd_only", "hd_low", "hd_high", "elf_top2_removed", "uloq_censored"), overwrite)
    reviewer_model_sensitivity(TRUE)
    reviewer_diagnostics(TRUE)
    write_reviewer_index()
  } else {
    reviewer_usage()
    stop("Unknown reviewer action: ", action, call. = FALSE)
  }
  invisible(TRUE)
}

# ---- sampling-density grouping -----------------------------------------------
