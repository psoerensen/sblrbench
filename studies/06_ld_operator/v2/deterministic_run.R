.study06v2_atomic_rds <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  tmp <- tempfile(".study06v2-", dirname(path), ".rds")
  saveRDS(x, tmp, compress = FALSE)
  if (!file.rename(tmp, path)) {
    unlink(tmp)
    stop("Atomic Study 06 v2 checkpoint write failed.", call. = FALSE)
  }
  path
}

.study06v2_write_csv <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  if (is.data.frame(x) && nrow(x))
    x <- x[do.call(order, unname(x)), , drop = FALSE]
  tmp <- tempfile(".study06v2-", dirname(path), ".csv")
  write.csv(x, tmp, row.names = FALSE, na = "")
  if (!file.rename(tmp, path)) {
    unlink(tmp)
    stop("Atomic Study 06 v2 CSV write failed.", call. = FALSE)
  }
  path
}

.study06v2_safe_v1_design <- function(store) {
  required <- c("study06_operator_stats", "study06_ld_bundle",
    "study06_scaled_genotypes", "study06_operator_simulation")
  out <- stats::setNames(lapply(required, targets::tar_read_raw,
    store = store), required)
  stats <- out$study06_operator_stats
  scaled <- out$study06_scaled_genotypes
  sim <- out$study06_operator_simulation
  if (!identical(stats$marker_names, colnames(scaled$train)) ||
      !identical(stats$marker_names, sim$marker_ids) ||
      !identical(stats$rows, scaled$train_rows) ||
      !identical(sim$simulation_seed, 61100L))
    stop("Safe v1-to-v2 design cache contract failed.", call. = FALSE)
  out
}

.study06v2_inspect <- function(stats, blocks, effects, eigen_prop) {
  sblr:::stblr_block_low_rank_contract_internal(
    bed_files = stats$bed_files, n_bed = as.integer(stats$n_bed),
    cls = stats$cls, rows = stats$rows,
    af = unlist(stats$af, use.names = FALSE),
    block_start = as.integer(blocks$start - 1L),
    wy = do.call(rbind, stats$wy), effects = as.numeric(effects),
    eigen_prop = eigen_prop)
}

.study06v2_run_deterministic <- function(config) {
  output <- file.path(config$local_dir, "deterministic")
  dir.create(output, recursive = TRUE, showWarnings = FALSE)
  store <- file.path("results", "local", "study06_ld_operator", "_targets")
  design <- .study06v2_safe_v1_design(store)
  stats <- design$study06_operator_stats
  scaled <- design$study06_scaled_genotypes
  sim <- design$study06_operator_simulation
  source(file.path("studies", "06_ld_operator", "blocks.R"), local = TRUE)
  blocks <- .study06_blocks(stats$marker_names, config$block_size)
  .study06_validate_blocks(blocks, stats$marker_names)
  historical <- read.csv(file.path(config$historical_capsules[1L],
    "selected_block_definitions.csv"), stringsAsFactors = FALSE)
  if (!identical(blocks$start, as.integer(historical$start)) ||
      !identical(blocks$end, as.integer(historical$end)))
    stop("V2 block boundaries differ from historical v1.", call. = FALSE)
  y <- as.numeric(sim$phenotype[rownames(scaled$train), 1L])
  direct_score <- as.numeric(crossprod(scaled$train, y))
  if (max(abs(direct_score - unlist(stats$wy, use.names = FALSE))) >
      config$operator_tolerance$product_absolute)
    stop("Direct training score differs from frozen summary statistics.",
      call. = FALSE)
  source_blocks <- lapply(seq_len(nrow(blocks)), function(b) {
    idx <- blocks$start[b]:blocks$end[b]
    crossprod(scaled$train[, idx, drop = FALSE])
  })
  specifications <- data.frame(configuration = c("low_rank_full",
    "low_rank_0999", "low_rank_0995"), eigen_prop = c(
      config$eigen_prop_full, config$eigen_props[["low_rank_0999"]],
      config$eigen_props[["low_rank_0995"]]), stringsAsFactors = FALSE)
  summaries <- list(); block_rows <- list()
  for (i in seq_len(nrow(specifications))) {
    specification <- specifications[i, ]
    controls <- list(representation = "low_rank",
      eigen_prop = specification$eigen_prop)
    .study06v2_assert_fit_spec(specification$configuration, controls, config)
    started <- proc.time()[["elapsed"]]
    inspect <- .study06v2_inspect(stats, blocks, sim$effects,
      specification$eigen_prop)
    elapsed <- proc.time()[["elapsed"]] - started
    source_oracle <- if (specification$configuration == "low_rank_full")
      source_blocks else NULL
    score_oracle <- if (specification$configuration == "low_rank_full")
      direct_score else NULL
    gate <- .study06v2_low_rank_identity_gate(inspect,
      as.numeric(sim$effects), as.numeric(stats$yy),
      source_blocks = source_oracle, source_score = score_oracle,
      tolerance = list(absolute = config$operator_tolerance$absolute,
        relative = config$operator_tolerance$relative, sse = 1e-6,
        update = config$operator_tolerance$absolute,
        source_matrix = config$operator_tolerance$absolute,
        source_score = config$operator_tolerance$product_absolute))
    if (specification$configuration == "low_rank_full")
      .study06v2_validate_full_positive_rank(inspect,
        specification$eigen_prop)
    diagnostics <- .study06v2_operator_diagnostics(inspect,
      specification$eigen_prop, specification$configuration)
    immutable_bytes <- sum(vapply(inspect$factor, object.size, numeric(1))) +
      sum(vapply(inspect$transformed_score, object.size, numeric(1)))
    chain_bytes <- object.size(inspect$residual)
    summaries[[i]] <- cbind(data.frame(
      configuration = specification$configuration,
      representation = "low_rank", operator_contract = "block_low_rank_v1",
      eigen_prop = specification$eigen_prop,
      block_count = length(inspect$factor),
      positive_rank = sum(diagnostics$positive_rank),
      retained_rank = sum(diagnostics$retained_rank),
      rank_fraction = sum(diagnostics$retained_rank) /
        sum(diagnostics$block_size),
      retained_mass_fraction = sum(diagnostics$retained_eigenvalue_mass) /
        sum(diagnostics$positive_eigenvalue_mass),
      work_storage_proxy = unique(diagnostics$work_storage_proxy)[1L],
      immutable_operator_bytes = as.numeric(immutable_bytes),
      per_chain_residual_bytes = as.numeric(chain_bytes),
      construction_seconds = elapsed, stringsAsFactors = FALSE), gate$summary)
    block_rows[[i]] <- cbind(data.frame(
      configuration = specification$configuration,
      eigen_prop = specification$eigen_prop, stringsAsFactors = FALSE),
      diagnostics)
    if (specification$configuration == "low_rank_full")
      .study06v2_atomic_rds(list(inspect = inspect, blocks = blocks,
        marker_ids = stats$marker_names, source_sha = config$required_sblr_sha,
        package_version = config$required_sblr_version,
        representation = "low_rank", eigen_prop = specification$eigen_prop),
        file.path(output, "near_full_operator_checkpoint.rds"))
    rm(inspect); gc()
  }
  .study06v2_write_csv(do.call(rbind, summaries),
    file.path(output, "deterministic_identity_summary.csv"))
  .study06v2_write_csv(do.call(rbind, block_rows),
    file.path(output, "low_rank_block_diagnostics.csv"))
  .study06v2_write_csv(blocks,
    file.path(output, "block_definitions.csv"))
  .study06v2_write_csv(.study06v2_design_crosswalk(config),
    file.path(output, "v1_to_v2_design_crosswalk.csv"))
  invisible(TRUE)
}
