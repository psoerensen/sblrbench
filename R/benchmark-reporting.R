# Shared, table-only presentation helpers for benchmark reports.

sblrbench_method_levels <- function() c(
  "st_bed_bayesc", "st_bed_bayesr", "st_csr_sbayesc", "st_csr_sbayesr"
)

sblrbench_method_labels <- function() c(
  st_bed_bayesc = "ST-BED BayesC", st_bed_bayesr = "ST-BED BayesR",
  st_csr_sbayesc = "ST-CSR SBayesC", st_csr_sbayesr = "ST-CSR SBayesR"
)

sblrbench_method_colours <- function() c(
  st_bed_bayesc = "#0072B2", st_bed_bayesr = "#D55E00",
  st_csr_sbayesc = "#009E73", st_csr_sbayesr = "#CC79A7"
)

sblrbench_method_shapes <- function() c(
  st_bed_bayesc = 16, st_bed_bayesr = 16,
  st_csr_sbayesc = 17, st_csr_sbayesr = 17
)

sblrbench_method_linetypes <- function() c(
  st_bed_bayesc = "solid", st_bed_bayesr = "solid",
  st_csr_sbayesc = "dashed", st_csr_sbayesr = "dashed"
)

sblrbench_architecture_levels <- function() c("sparse_homogeneous", "sparse_mixture")
sblrbench_architecture_labels <- function() c(
  sparse_homogeneous = "Sparse homogeneous",
  sparse_mixture = "Sparse variance mixture"
)

sblrbench_method_factor <- function(x) factor(x, levels = sblrbench_method_levels(),
  labels = unname(sblrbench_method_labels()[sblrbench_method_levels()]))
sblrbench_architecture_factor <- function(x) factor(x,
  levels = sblrbench_architecture_levels(),
  labels = unname(sblrbench_architecture_labels()[sblrbench_architecture_levels()]))

sblrbench_model_specification_labels <- function() c(
  matched = "Matched prior class", misspecified = "Misspecified prior class"
)

format_sblrbench_number <- function(x, digits = 3L) {
  ifelse(is.na(x), "—", formatC(x, format = "f", digits = digits, big.mark = ","))
}
format_sblrbench_proportion <- function(x, digits = 3L) format_sblrbench_number(x, digits)
format_sblrbench_runtime <- function(x) {
  ifelse(is.na(x), "—", ifelse(x < 60, sprintf("%.1f s", x), sprintf("%.1f min", x / 60)))
}
format_sblrbench_interval <- function(mean, lower, upper, digits = 3L) paste0(
  format_sblrbench_number(mean, digits), " [",
  format_sblrbench_number(lower, digits), ", ",
  format_sblrbench_number(upper, digits), "]"
)
format_sblrbench_status <- function(x) {
  labels <- c(ok = "Passed", pass = "Passed", complete = "Complete",
    failed = "Failed", fail = "Failed", unavailable = "Unavailable",
    indeterminate = "Indeterminate")
  out <- unname(labels[as.character(x)]); out[is.na(out)] <- as.character(x)[is.na(out)]; out
}

theme_sblrbench <- function(base_size = 11) {
  ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_line(colour = "#E5E7EB", linewidth = .3),
      plot.title.position = "plot", plot.title = ggplot2::element_text(face = "bold"),
      plot.caption.position = "plot", plot.caption = ggplot2::element_text(hjust = 0),
      strip.text = ggplot2::element_text(face = "bold"),
      strip.background = ggplot2::element_rect(fill = "#F3F4F6", colour = NA),
      panel.spacing = grid::unit(.8, "lines"), legend.position = "bottom",
      plot.margin = ggplot2::margin(8, 10, 8, 8)
    )
}

sblrbench_method_scales <- function() list(
  ggplot2::scale_colour_manual(values = sblrbench_method_colours(),
    breaks = sblrbench_method_levels(), labels = sblrbench_method_labels()),
  ggplot2::scale_shape_manual(values = sblrbench_method_shapes(),
    breaks = sblrbench_method_levels(), labels = sblrbench_method_labels()),
  ggplot2::scale_linetype_manual(values = sblrbench_method_linetypes(),
    breaks = sblrbench_method_levels(), labels = sblrbench_method_labels())
)

prepare_sblrbench_replicates <- function(data, group_cols, value_col = "value",
                                         replicate_col = "replicate") {
  needed <- c(group_cols, value_col, replicate_col)
  if (!all(needed %in% names(data))) stop("Replicate data are missing required columns.", call. = FALSE)
  observations <- data[needed]
  key <- interaction(observations[group_cols], drop = TRUE, lex.order = TRUE)
  pieces <- split(observations, key)
  summary <- do.call(rbind, lapply(pieces, function(z) {
    values <- z[[value_col]]
    row <- z[1L, group_cols, drop = FALSE]
    row$n_replicates <- length(unique(z[[replicate_col]]))
    row$mean <- mean(values); row$sd <- if (row$n_replicates > 1L) stats::sd(values) else NA_real_
    row$minimum <- min(values); row$maximum <- max(values); row
  }))
  rownames(summary) <- NULL
  list(observations = observations, summary = summary)
}

.sblrbench_script_text <- function(path) {
  if (length(path) != 1L || !file.exists(path)) stop("Capsule script does not exist: ", path, call. = FALSE)
  paste(readLines(path, warn = FALSE, encoding = "UTF-8"), collapse = "\n")
}

display_capsule_script <- function(path) {
  cat("```r\n", .sblrbench_script_text(path), "\n```\n", sep = "")
  invisible(path)
}

display_capsule_script_collapsed <- function(path, summary = "Show developer contract smoke test") {
  cat("<details><summary>", summary, "</summary>\n\n```r\n",
    .sblrbench_script_text(path), "\n```\n\n</details>\n", sep = "")
  invisible(path)
}
