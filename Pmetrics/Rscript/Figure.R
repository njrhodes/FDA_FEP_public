# =============================================================================
# FDA-FEP manuscript figures
# Produces manuscript-facing figures from the final source runs and PTA simulations:
#   Figure 1, Figures 2-4, and Figures S1-S5
# =============================================================================

FIGURE_SCRIPT_VERSION <- "3.2.0"

ROOT <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
ANALYSIS_SCRIPT <- file.path("Pmetrics", "Rscript", "Analysis.R")
FIGURE_DIR <- file.path("Manuscript", "Figures")

if (!file.exists(ANALYSIS_SCRIPT)) {
  stop("Run Figure.R from the FDA_FEP_public repository root.", call. = FALSE)
}

dir.create(FIGURE_DIR, recursive = TRUE, showWarnings = FALSE)
source(ANALYSIS_SCRIPT, local = .GlobalEnv)

required_packages <- c(
  "ggplot2",
  "patchwork",
  "scales",
  "viridisLite",
  "DiagrammeR",
  "DiagrammeRsvg"
)
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1L), quietly = TRUE)
]
if (length(missing_packages)) {
  stop("Install required package(s): ", paste(missing_packages, collapse = ", "), call. = FALSE)
}

assert_pmetrics()

MAIN_MICS <- c(1, 2, 4, 8)
MIC_GRID_FIGURE <- c(0.25, 0.5, 1, 2, 4, 8, 16, 32)
TARGET_LEVELS <- c("50% fT>MIC", "68% fT>MIC", "100% fT>MIC")
WEIGHT_LEVELS <- c(60, 80, 120)

# Exact colors used in the manuscript bar figures.
TARGET_BAR_COLORS <- c(
  "50% fT>MIC" = "dodgerblue",
  "68% fT>MIC" = "seagreen3",
  "100% fT>MIC" = "indianred"
)

# Exact target-color logic used in the manuscript MIC-continuum figures.
TARGET_LINE_COLORS <- stats::setNames(
  viridisLite::viridis(3, option = "D"),
  TARGET_LEVELS
)

REGIMENS <- list(
  II = list(
    analysis = "2 g q12 intermittent infusion + loading dose",
    main_stem = "Figure2_PTA_II",
    supp_stem = "FigureS3_PTA_by_MIC_II",
    panel = "2 g q12 II + LD"
  ),
  EI = list(
    analysis = "2 g q12 extended infusion + loading dose",
    main_stem = "Figure3_PTA_EI",
    supp_stem = "FigureS4_PTA_by_MIC_EI",
    panel = "2 g q12 EI + LD"
  ),
  CI = list(
    analysis = "4 g q24 continuous infusion + loading dose",
    main_stem = "Figure4_PTA_CI",
    supp_stem = "FigureS5_PTA_by_MIC_CI",
    panel = "4 g q24 CI + LD"
  )
)

as_numeric_safe <- function(x) suppressWarnings(as.numeric(as.character(x)))

save_figure <- function(plot, stem, width, height) {
  svg <- file.path(FIGURE_DIR, paste0(stem, ".svg"))
  tif <- file.path(FIGURE_DIR, paste0(stem, ".tiff"))

  ggplot2::ggsave(svg, plot, width = width, height = height, units = "in", bg = "white")
  ggplot2::ggsave(
    tif, plot, width = width, height = height, units = "in",
    dpi = 600, compression = "lzw", bg = "white"
  )

  message("Wrote ", svg)
  message("Wrote ", tif)
  invisible(c(svg = svg, tiff = tif))
}

# =============================================================================
# Figure 1
# =============================================================================

prepare_op_data <- function(run, outputs = c(1L, 2L, 5L)) {
  dat <- as.data.frame(run$op$data, check.names = FALSE)
  names(dat) <- tolower(names(dat))

  needed <- c("id", "outeq", "obs", "pred", "pred.type", "icen")
  missing <- setdiff(needed, names(dat))
  if (length(missing)) {
    stop("Unexpected OP data; missing: ", paste(missing, collapse = ", "), call. = FALSE)
  }

  for (nm in intersect(c("id", "time", "outeq", "obs", "pred"), names(dat))) {
    dat[[nm]] <- as_numeric_safe(dat[[nm]])
  }

  dat[
    is.finite(dat$id) &
      as.character(dat$pred.type) == "post" &
      as.character(dat$icen) == "median" &
      dat$outeq %in% outputs &
      is.finite(dat$obs) &
      is.finite(dat$pred),
    , drop = FALSE
  ]
}

op_statistics <- function(dat) {
  if (nrow(dat) < 3L || length(unique(dat$pred)) < 2L) {
    return(c(r2 = NA_real_, intercept = NA_real_, slope = NA_real_))
  }
  fit <- stats::lm(obs ~ pred, data = dat)
  co <- stats::coef(fit)
  c(
    r2 = summary(fit)$r.squared,
    intercept = unname(co[["(Intercept)"]]),
    slope = unname(co[["pred"]])
  )
}

op_limit <- function(dat, step = 10) {
  values <- c(dat$obs, dat$pred)
  values <- values[is.finite(values)]
  if (!length(values)) return(step)
  ceiling(max(values) / step) * step
}

make_op_panel <- function(dat, output, title, letter, limit) {
  panel <- dat[dat$outeq == output, , drop = FALSE]
  stats <- op_statistics(panel)
  stats_label <- sprintf(
    "R² = %.3f\nIntercept = %.2f\nSlope = %.3f",
    stats[["r2"]], stats[["intercept"]], stats[["slope"]]
  )

  ggplot2::ggplot(panel, ggplot2::aes(pred, obs)) +
    ggplot2::geom_abline(
      slope = 1, intercept = 0, linetype = "dashed",
      color = "grey45", linewidth = 0.65
    ) +
    # The manuscript uses one gold point style for every displayed observation.
    ggplot2::geom_point(
      shape = 21, fill = "goldenrod2", color = "grey25",
      size = 2.7, stroke = 0.55, alpha = 0.95
    ) +
    ggplot2::annotate(
      "text", x = Inf, y = -Inf, label = stats_label,
      hjust = 1.06, vjust = -0.10, size = 3.15, lineheight = 1.15
    ) +
    ggplot2::coord_equal(
      xlim = c(0, limit), ylim = c(0, limit), expand = FALSE
    ) +
    ggplot2::scale_x_continuous(breaks = scales::breaks_pretty(n = 5)) +
    ggplot2::scale_y_continuous(breaks = scales::breaks_pretty(n = 5)) +
    ggplot2::labs(
      title = paste0(letter, ". ", title),
      x = "Predicted (mg/L)",
      y = "Observed (mg/L)"
    ) +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 12, hjust = 0),
      panel.grid.major = ggplot2::element_line(color = "grey90", linewidth = 0.35),
      panel.grid.minor = ggplot2::element_blank(),
      panel.border = ggplot2::element_rect(color = "grey35", linewidth = 0.65),
      axis.title = ggplot2::element_text(size = 10.5),
      axis.text = ggplot2::element_text(size = 9.5),
      plot.margin = ggplot2::margin(4, 5, 4, 4)
    )
}

make_figure1 <- function(write_files = TRUE) {
  dev <- prepare_op_data(load_manuscript_run("base"))
  val <- prepare_op_data(load_manuscript_run("final_validation"))

  dev_lim <- op_limit(dev)
  val_lim <- op_limit(val)

  figure <- patchwork::wrap_plots(
    make_op_panel(dev, 1L, "Pre-filter plasma", "A", dev_lim),
    make_op_panel(dev, 2L, "Post-filter plasma", "B", dev_lim),
    make_op_panel(dev, 5L, "ELF", "C", dev_lim),
    make_op_panel(val, 1L, "Pre-filter plasma", "D", val_lim),
    make_op_panel(val, 2L, "Post-filter plasma", "E", val_lim),
    make_op_panel(val, 5L, "ELF", "F", val_lim),
    ncol = 3, nrow = 2, byrow = TRUE
  )

  if (write_files) save_figure(figure, "Figure1_observed_vs_posterior_predicted", 10.5, 7.0)
  invisible(figure)
}

# =============================================================================
# PTA data
# =============================================================================

prepare_pta_data <- function() {
  dat <- simulate_pta(write_html = FALSE)
  dat$Target <- factor(dat$Target, levels = TARGET_LEVELS)
  dat$Weight_kg <- factor(dat$Weight_kg, levels = WEIGHT_LEVELS)
  dat
}

get_regimen_data <- function(dat, code) {
  out <- dat[dat$Regimen == REGIMENS[[code]]$analysis, , drop = FALSE]
  if (!nrow(out)) stop("No PTA data for regimen ", code, call. = FALSE)
  out
}

# =============================================================================
# Figures 2-4: exact 2 x 2 MIC-panel layout from manuscript pages 39-41
# =============================================================================

make_main_mic_panel <- function(dat, code, mic) {
  spec <- REGIMENS[[code]]
  panel <- dat[as_numeric_safe(dat$MIC_mg_L) == mic, , drop = FALSE]
  if (!nrow(panel)) stop("No rows for MIC ", mic, call. = FALSE)

  panel$Matrix <- factor(panel$Matrix, levels = c("ELF", "Plasma"))
  panel$Weight <- factor(
    paste0(as.character(panel$Weight_kg), "kg"),
    levels = paste0(WEIGHT_LEVELS, "kg")
  )
  panel$TargetShort <- factor(
    sub(" fT>MIC$", "", as.character(panel$Target)),
    levels = c("50%", "68%", "100%")
  )

  ggplot2::ggplot(panel, ggplot2::aes(TargetShort, PTA, fill = Target)) +
    ggplot2::geom_col(width = 0.68) +
    ggplot2::geom_errorbar(
      ggplot2::aes(ymin = Lower_95, ymax = Upper_95),
      width = 0.18, linewidth = 0.38
    ) +
    ggplot2::geom_hline(
      yintercept = 0.90, linetype = "dashed",
      color = "grey25", linewidth = 0.55
    ) +
    ggplot2::facet_grid(Matrix ~ Weight) +
    ggplot2::scale_fill_manual(values = TARGET_BAR_COLORS, drop = FALSE) +
    ggplot2::scale_y_continuous(
      limits = c(0, 1), breaks = seq(0, 1, 0.2),
      labels = scales::percent_format(accuracy = 1),
      expand = ggplot2::expansion(mult = c(0, 0.02))
    ) +
    ggplot2::labs(
      title = paste0(spec$panel, ": PTA at MIC = ", mic, " mg/L"),
      subtitle = "Bodyweight strata × 50%, 68%, 100% fT>MIC",
      x = "Pharmacodynamic target",
      y = "Probability of target attainment",
      fill = "Target"
    ) +
    ggplot2::theme_minimal(base_size = 8.6) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 9.6, hjust = 0),
      plot.subtitle = ggplot2::element_text(size = 8.0, hjust = 0),
      strip.text = ggplot2::element_text(face = "bold", size = 8.0),
      strip.background = ggplot2::element_blank(),
      panel.border = ggplot2::element_rect(color = "black", fill = NA, linewidth = 0.55),
      panel.grid.major.x = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(size = 7.3),
      axis.text.y = ggplot2::element_text(size = 7.3),
      axis.title = ggplot2::element_text(size = 8.2),
      legend.position = "bottom",
      legend.title = ggplot2::element_text(face = "bold"),
      plot.margin = ggplot2::margin(3, 4, 3, 3)
    )
}

make_main_figure <- function(dat, code, write_files = TRUE) {
  regimen <- get_regimen_data(dat, code)
  panels <- lapply(MAIN_MICS, function(mic) make_main_mic_panel(regimen, code, mic))

  figure <- patchwork::wrap_plots(
    panels,
    ncol = 2, nrow = 2, byrow = TRUE,
    guides = "collect"
  ) +
    patchwork::plot_annotation(tag_levels = "A") &
    ggplot2::theme(
      legend.position = "bottom",
      plot.tag = ggplot2::element_text(size = 11, face = "plain")
    )

  if (write_files) save_figure(figure, REGIMENS[[code]]$main_stem, 12.0, 8.3)
  invisible(figure)
}

make_figures2_to_4 <- function(write_files = TRUE) {
  dat <- prepare_pta_data()
  invisible(lapply(c("II", "EI", "CI"), function(code) {
    make_main_figure(dat, code, write_files)
  }))
}

# =============================================================================
# Figures S3-S5: plasma/ELF MIC-continuum layout
# =============================================================================

series_levels <- unlist(lapply(TARGET_LEVELS, function(target) {
  paste0(target, "; ", WEIGHT_LEVELS, "kg")
}), use.names = FALSE)

series_colors <- stats::setNames(
  rep(unname(TARGET_LINE_COLORS), each = length(WEIGHT_LEVELS)),
  series_levels
)
series_linetypes <- stats::setNames(
  rep(c("solid", "dashed", "dotted"), times = length(TARGET_LEVELS)),
  series_levels
)
series_shapes <- stats::setNames(
  rep(c(16, 15, 18), times = length(TARGET_LEVELS)),
  series_levels
)

make_supp_panel <- function(dat, code, matrix, letter) {
  spec <- REGIMENS[[code]]
  panel <- dat[dat$Matrix == matrix, , drop = FALSE]
  panel$Series <- factor(
    paste0(as.character(panel$Target), "; ", as.character(panel$Weight_kg), "kg"),
    levels = series_levels
  )

  ggplot2::ggplot(
    panel,
    ggplot2::aes(
      x = MIC_mg_L, y = PTA,
      color = Series, linetype = Series, shape = Series,
      group = Series
    )
  ) +
    ggplot2::geom_hline(
      yintercept = 0.90, linetype = "dashed",
      color = "grey25", linewidth = 0.65
    ) +
    ggplot2::geom_line(linewidth = 0.75) +
    ggplot2::geom_point(size = 2.4, stroke = 0.35) +
    ggplot2::scale_x_log10(
      breaks = MIC_GRID_FIGURE,
      labels = as.character(MIC_GRID_FIGURE),
      limits = range(MIC_GRID_FIGURE),
      expand = ggplot2::expansion(mult = c(0.02, 0.04))
    ) +
    ggplot2::scale_y_continuous(
      limits = c(-0.06, 1.02), breaks = seq(0, 1, 0.2),
      labels = scales::percent_format(accuracy = 1),
      expand = ggplot2::expansion(mult = c(0, 0))
    ) +
    ggplot2::scale_color_manual(values = series_colors, drop = FALSE) +
    ggplot2::scale_linetype_manual(values = series_linetypes, drop = FALSE) +
    ggplot2::scale_shape_manual(values = series_shapes, drop = FALSE) +
    ggplot2::labs(
      title = paste0(letter, ". ", matrix, " PTA: ", spec$panel),
      x = "FEP MIC (mg/L)",
      y = "Proportion Achieving Goal (%)",
      color = NULL, linetype = NULL, shape = NULL
    ) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 10.8, hjust = 0),
      panel.grid.major = ggplot2::element_line(color = "grey91", linewidth = 0.35),
      panel.grid.minor = ggplot2::element_blank(),
      axis.line = ggplot2::element_line(color = "grey45", linewidth = 0.45),
      axis.title = ggplot2::element_text(size = 10.5),
      axis.text = ggplot2::element_text(size = 9.5),
      legend.position = "bottom",
      legend.text = ggplot2::element_text(size = 7.1),
      legend.key.width = grid::unit(1.25, "lines"),
      plot.margin = ggplot2::margin(4, 5, 3, 4)
    )
}

make_supp_figure <- function(dat, code, write_files = TRUE) {
  regimen <- get_regimen_data(dat, code)

  figure <- patchwork::wrap_plots(
    make_supp_panel(regimen, code, "Plasma", "A"),
    make_supp_panel(regimen, code, "ELF", "B"),
    nrow = 1,
    guides = "collect"
  ) &
    ggplot2::theme(legend.position = "bottom")

  if (write_files) save_figure(figure, REGIMENS[[code]]$supp_stem, 12.0, 5.6)
  invisible(figure)
}

make_figures_s3_to_s5 <- function(write_files = TRUE) {
  dat <- prepare_pta_data()
  invisible(lapply(c("II", "EI", "CI"), function(code) {
    make_supp_figure(dat, code, write_files)
  }))
}

# =============================================================================
# Figure S1. Structural cefepime population pharmacokinetic model
#
# Structural schematic corresponding to the final piecewise model.
# =============================================================================

save_grviz <- function(grviz_obj, stem, out_dir = FIGURE_DIR) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  svg_file <- file.path(out_dir, paste0(stem, ".svg"))
  svg_text <- DiagrammeRsvg::export_svg(grviz_obj)

  writeLines(svg_text, svg_file, useBytes = TRUE)
  message("Wrote ", svg_file)

  invisible(svg_file)
}

make_figure_s1_structure <- function(write_file = TRUE) {
  figure <- DiagrammeR::grViz("
    digraph cef_struct_minimal {
      graph [
        layout=neato,
        overlap=false,
        splines=true,
        outputorder=edgesfirst,
        bgcolor='white',
        margin=0.05,
        pad=0.20
      ]

      node [
        shape=box,
        fontsize=12,
        fontname='Arial',
        style='filled',
        penwidth=1.1,
        color='black'
      ]

      edge [
        fontsize=10,
        fontname='Arial',
        arrowsize=0.7,
        penwidth=1,
        color='black'
      ]

      dose [
        label='Dose\\n(IV infusion)',
        style='rounded,filled',
        fillcolor='#D3ECE4',
        pos='0,2.6!'
      ]

      central [
        label='Pre-filter\\nplasma\\nObs: C_pre',
        fillcolor='#FADADD',
        pos='0,0!'
      ]

      peripheral [
        label='Peripheral\\n(no Obs)',
        fillcolor='#D4E4E8',
        pos='-2.6,1.6!'
      ]

      elf [
        label='ELF\\nObs: C_elf',
        fillcolor='#D4E4E8',
        pos='-2.6,-1.6!'
      ]

      post [
        label='Post-filter\\nplasma\\nObs: C_post',
        fillcolor='#FADADD',
        pos='2.6,1.6!'
      ]

      eff [
        label='Effluent\\nObs: C_eff',
        fillcolor='#FFF2C6',
        pos='0,-2.6!'
      ]

      out [
        label='Non-CRRT\\nElimination',
        fillcolor='#D4E4E8',
        pos='2.6,-1.6!'
      ]

      node [
        shape=plaintext,
        style='',
        fontsize=10,
        fontname='Arial'
      ]

      lab_k12 [
        label='k12, k21',
        pos='-1.4,1.1!'
      ]

      lab_k15 [
        label='k15, k51',
        pos='-1.4,-1.1!'
      ]

      lab_spost [
        label='S_post',
        pos='1.4,1.1!'
      ]

      lab_clsys [
        label='CL_systemic',
        pos='1.4,-1.1!'
      ]

      dose:s     -> central:n     [label='Rate_in']
      central:nw -> peripheral:se [dir=both]
      central:sw -> elf:ne        [dir=both]
      central:ne -> post:sw       [dir=both]
      central:s  -> eff:n         [label='CL_CRRT']
      central:se -> out:nw
    }
  ")

  if (write_file) {
    save_grviz(
      figure,
      "FigureS1_structural_model"
    )
  }

  invisible(figure)
}


# =============================================================================
# Supplemental final-model posterior residual diagnostics
# Pmetrics OP filtering: pred.type = "post", icen = "median"
# OUTEQ 1 = pre-filter plasma; 2 = post-filter plasma; 5 = ELF
# =============================================================================

prepare_post_residual_data <- function(run, output) {
  if (!output %in% c(1L, 2L, 5L)) {
    stop("OUTEQ must be 1, 2, or 5.", call. = FALSE)
  }

  dat <- as.data.frame(run$op$data, check.names = FALSE)
  names(dat) <- tolower(names(dat))

  needed <- c(
    "id", "time", "outeq", "obs", "pred",
    "pred.type", "icen", "wd"
  )
  missing <- setdiff(needed, names(dat))

  if (length(missing)) {
    stop(
      "Unexpected OP data; missing: ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }

  for (nm in c("id", "time", "outeq", "obs", "pred", "wd")) {
    dat[[nm]] <- as_numeric_safe(dat[[nm]])
  }

  dat[
    is.finite(dat$id) &
      as.character(dat$pred.type) == "post" &
      as.character(dat$icen) == "median" &
      dat$outeq == output &
      is.finite(dat$obs) &
      is.finite(dat$pred) &
      is.finite(dat$wd),
    ,
    drop = FALSE
  ]
}


residual_output_title <- function(output) {
  switch(
    as.character(output),
    `1` = "Pre-filter plasma",
    `2` = "Post-filter plasma",
    `5` = "ELF",
    stop("Unsupported OUTEQ.", call. = FALSE)
  )
}


make_supplemental_model_figures <- function(write_files = TRUE) {
  make_figure_s_validation_residuals(write_files)
  make_figure_s1_structure(write_files)
  invisible(TRUE)
}


# =============================================================================
# Supplemental posterior residual diagnostics — validation cohort only
#
# A-C: WPE versus posterior individual prediction
# D-F: WPE versus time
#
# pred.type = "post"
# icen      = "median"
# OUTEQ 1   = pre-filter plasma
# OUTEQ 2   = post-filter plasma
# OUTEQ 5   = ELF
# =============================================================================

validation_residual_limit <- function(run) {
  dat <- do.call(
    rbind,
    lapply(
      c(1L, 2L, 5L),
      function(output) prepare_post_residual_data(run, output)
    )
  )

  vals <- dat$wd[is.finite(dat$wd)]

  if (!length(vals)) return(1)

  lim <- ceiling(max(abs(vals), na.rm = TRUE))

  max(1, lim)
}


make_validation_resid_pred_panel <- function(
  run,
  output,
  letter,
  y_limit
) {
  panel <- prepare_post_residual_data(run, output)

  ggplot2::ggplot(
    panel,
    ggplot2::aes(x = pred, y = wd)
  ) +
    ggplot2::geom_hline(
      yintercept = 0,
      linetype = "dashed",
      color = "grey45",
      linewidth = 0.65
    ) +
    ggplot2::geom_point(
      shape = 21,
      fill = "goldenrod2",
      color = "grey25",
      size = 2.7,
      stroke = 0.55,
      alpha = 0.95
    ) +
    ggplot2::scale_x_continuous(
      breaks = scales::breaks_pretty(n = 5)
    ) +
    ggplot2::scale_y_continuous(
      limits = c(-y_limit, y_limit),
      breaks = scales::breaks_pretty(n = 5)
    ) +
    ggplot2::labs(
      title = paste0(
        letter,
        ". ",
        residual_output_title(output)
      ),
      x = "Posterior predicted (mg/L)",
      y = "Weighted prediction error"
    ) +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        face = "bold",
        size = 12,
        hjust = 0
      ),
      panel.grid.major = ggplot2::element_line(
        color = "grey90",
        linewidth = 0.35
      ),
      panel.grid.minor = ggplot2::element_blank(),
      panel.border = ggplot2::element_rect(
        color = "grey35",
        linewidth = 0.65
      ),
      axis.title = ggplot2::element_text(size = 10.5),
      axis.text = ggplot2::element_text(size = 9.5),
      plot.margin = ggplot2::margin(4, 5, 4, 4)
    )
}


make_validation_resid_time_panel <- function(
  run,
  output,
  letter,
  y_limit
) {
  panel <- prepare_post_residual_data(run, output)

  ggplot2::ggplot(
    panel,
    ggplot2::aes(x = time, y = wd)
  ) +
    ggplot2::geom_hline(
      yintercept = 0,
      linetype = "dashed",
      color = "grey45",
      linewidth = 0.65
    ) +
    ggplot2::geom_point(
      shape = 21,
      fill = "goldenrod2",
      color = "grey25",
      size = 2.7,
      stroke = 0.55,
      alpha = 0.95
    ) +
    ggplot2::scale_x_continuous(
      breaks = scales::breaks_pretty(n = 5)
    ) +
    ggplot2::scale_y_continuous(
      limits = c(-y_limit, y_limit),
      breaks = scales::breaks_pretty(n = 5)
    ) +
    ggplot2::labs(
      title = paste0(
        letter,
        ". ",
        residual_output_title(output)
      ),
      x = "Time (h)",
      y = "Weighted prediction error"
    ) +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        face = "bold",
        size = 12,
        hjust = 0
      ),
      panel.grid.major = ggplot2::element_line(
        color = "grey90",
        linewidth = 0.35
      ),
      panel.grid.minor = ggplot2::element_blank(),
      panel.border = ggplot2::element_rect(
        color = "grey35",
        linewidth = 0.65
      ),
      axis.title = ggplot2::element_text(size = 10.5),
      axis.text = ggplot2::element_text(size = 9.5),
      plot.margin = ggplot2::margin(4, 5, 4, 4)
    )
}


make_figure_s_validation_residuals <- function(
  write_files = TRUE
) {
  val <- load_manuscript_run("final_validation")

  y_limit <- validation_residual_limit(val)

  figure <- patchwork::wrap_plots(
    make_validation_resid_pred_panel(
      val, 1L, "A", y_limit
    ),
    make_validation_resid_pred_panel(
      val, 2L, "B", y_limit
    ),
    make_validation_resid_pred_panel(
      val, 5L, "C", y_limit
    ),
    make_validation_resid_time_panel(
      val, 1L, "D", y_limit
    ),
    make_validation_resid_time_panel(
      val, 2L, "E", y_limit
    ),
    make_validation_resid_time_panel(
      val, 5L, "F", y_limit
    ),
    ncol = 3,
    nrow = 2,
    byrow = TRUE
  ) +
    patchwork::plot_annotation(
      title = paste(
        "Final-model posterior residual diagnostics:",
        "validation cohort"
      ),
      subtitle = paste(
        "A-C, weighted prediction error versus posterior prediction;",
        "D-F, weighted prediction error versus time"
      ),
      theme = ggplot2::theme(
        plot.title = ggplot2::element_text(
          face = "bold",
          size = 13
        ),
        plot.subtitle = ggplot2::element_text(
          size = 10.5
        )
      )
    )

  if (write_files) {
    save_figure(
      figure,
      "FigureS2_validation_posterior_residual_diagnostics",
      10.5,
      7.0
    )
  }

  invisible(figure)
}

usage <- function() {
  cat(
    "FDA-FEP Figure.R ", FIGURE_SCRIPT_VERSION, "\n\n",
    "Usage:\n",
    "  Rscript Pmetrics/Rscript/Figure.R figure1\n",
    "  Rscript Pmetrics/Rscript/Figure.R main-pta\n",
    "  Rscript Pmetrics/Rscript/Figure.R supplemental-pta\n",
    "  Rscript Pmetrics/Rscript/Figure.R supplemental-model\n",
    "  Rscript Pmetrics/Rscript/Figure.R all\n",
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
  message("FDA-FEP Figure.R ", FIGURE_SCRIPT_VERSION)

  if (command == "figure1") {
    make_figure1(TRUE)
  } else if (command == "main-pta") {
    make_figures2_to_4(TRUE)
  } else if (command == "supplemental-pta") {
    make_figures_s3_to_s5(TRUE)
  } else if (command == "supplemental-model") {
    make_supplemental_model_figures(TRUE)
  } else if (command == "all") {
    make_figure1(TRUE)
    dat <- prepare_pta_data()
    for (code in c("II", "EI", "CI")) {
      make_main_figure(dat, code, TRUE)
      make_supp_figure(dat, code, TRUE)
    }
    make_supplemental_model_figures(TRUE)
  } else {
    usage()
    stop("Unknown command: ", command, call. = FALSE)
  }

  invisible(TRUE)
}

if (sys.nframe() == 0L) main()
