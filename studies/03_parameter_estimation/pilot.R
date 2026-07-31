.study03_paths <- function() list(
  glist_path = Sys.getenv("SBLR_BENCH_GLIST", ""),
  data_dir = Sys.getenv("SBLR_BENCH_DATA_DIR", file.path("results", "local", "03_parameter_estimation", "data")),
  output_dir = Sys.getenv("SBLR_BENCH_LD_DIR",
    file.path("results", "local", "03_parameter_estimation", "genotype_setup")))

.study03_summary_stats <- function(simulation, Glist, config) {
  stats <- sblr::make_summary_stats(Glist = Glist, y = simulation$truth$phenotypes,
    chr = config$chr, rows = seq_len(nrow(simulation$truth$phenotypes)), scale = TRUE, nthreads = 1L)
  if (!identical(stats$marker_names, simulation$data$marker_ids) ||
      !identical(stats$trait_names, config$trait) ||
      !identical(as.integer(stats$n), nrow(simulation$truth$phenotypes)))
    stop("Study 03 summary statistics are not aligned.", call. = FALSE)
  stats
}

.study03_write_csv <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(x, path, row.names = FALSE, na = "")
  path
}

.study03_paired_summary <- function(x) do.call(rbind, lapply(split(x,
  interaction(x$architecture, x$comparison_id, x$estimand_id, drop = TRUE)), function(z) {
    v <- z[z$complete_pair, , drop = FALSE]
    data.frame(architecture = z$architecture[1], comparison_id = z$comparison_id[1],
      comparison_type = z$comparison_type[1], estimand_id = z$estimand_id[1],
      replicate_count = length(unique(z$replicate)), successful_replicates = nrow(v),
      mean_estimate_difference = if (nrow(v)) mean(v$estimate_difference) else NA_real_,
      sd_estimate_difference = if (nrow(v) > 1L) stats::sd(v$estimate_difference) else NA_real_,
      median_estimate_difference = if (nrow(v)) stats::median(v$estimate_difference) else NA_real_,
      minimum_estimate_difference = if (nrow(v)) min(v$estimate_difference) else NA_real_,
      maximum_estimate_difference = if (nrow(v)) max(v$estimate_difference) else NA_real_,
      mean_absolute_error_difference = if (nrow(v)) mean(v$absolute_error_difference) else NA_real_,
      mean_interval_width_difference = if (nrow(v)) mean(v$interval_width_difference) else NA_real_,
      coverage_agreement_observation = if (nrow(v) == 1L) v$coverage_agreement else NA,
      stringsAsFactors = FALSE)
  }))
