.study03_paths <- function() list(
  glist_path = Sys.getenv("SBLR_BENCH_GLIST", ""),
  data_dir = Sys.getenv("SBLR_BENCH_DATA_DIR", file.path("results", "local", "03_parameter_estimation", "data")),
  output_dir = file.path("results", "local", "03_parameter_estimation", "genotype_setup"))

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
  interaction(x$architecture, x$comparison_id, x$estimand_id, drop = TRUE)), function(z)
    data.frame(architecture = z$architecture[1], comparison_id = z$comparison_id[1],
      comparison_type = z$comparison_type[1], estimand_id = z$estimand_id[1],
      replicate_count = length(unique(z$replicate)), mean_estimate_difference = mean(z$estimate_difference),
      mean_absolute_error_difference = mean(z$absolute_error_difference),
      mean_interval_width_difference = mean(z$interval_width_difference),
      coverage_agreement_observation = if (nrow(z) == 1L) z$coverage_agreement else NA,
      stringsAsFactors = FALSE)))
