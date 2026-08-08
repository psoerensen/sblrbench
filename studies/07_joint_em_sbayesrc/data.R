study07_assert <- function(condition, message) {
  if (!isTRUE(condition)) stop(message, call. = FALSE)
  invisible(TRUE)
}

study07_sha256 <- function(path) {
  unname(digest::digest(file = path, algo = "sha256", serialize = FALSE))
}

study07_write_block_csr <- function(truth, output_dir) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  prefix <- file.path(output_dir, "study06_exact_block")
  suffix <- c(".row_ptr.u64.bin", ".col_idx.u32.0based.bin",
    ".values.f32.bin", ".meta.txt")
  paths <- paste0(prefix, suffix)
  if (!all(file.exists(paths))) {
    X <- as.matrix(truth$training_x)
    block <- as.integer(truth$block)
    m <- ncol(X)
    rows <- cols <- integer()
    values <- numeric()
    for (b in unique(block)) {
      index <- which(block == b)
      correlation <- crossprod(X[, index, drop = FALSE]) / nrow(X)
      edge <- which(upper.tri(correlation), arr.ind = TRUE)
      rows <- c(rows, index[edge[, 1L]])
      cols <- c(cols, index[edge[, 2L]])
      values <- c(values, correlation[edge])
    }
    order <- order(rows, cols)
    rows <- rows[order]; cols <- cols[order]; values <- values[order]
    row_ptr <- c(0, cumsum(tabulate(rows, nbins = m)))
    sblr:::.stblr_write_uint64_file(paths[[1L]], row_ptr)
    sblr:::.stblr_write_uint32_file(paths[[2L]], cols - 1L)
    connection <- file(paths[[3L]], open = "wb")
    tryCatch(writeBin(as.numeric(values), connection, size = 4L,
      endian = .Platform$endian), finally = close(connection))
    writeLines(c(
      "format=sparse_ld_csr", "storage=streamed_upper_triangle",
      "source=frozen_Study06_training_x_exact_registered_blocks",
      paste0("n_bed=", nrow(X)), paste0("n_used=", nrow(X)),
      paste0("n_samples_used=", nrow(X)), paste0("n_variants=", m),
      paste0("nnz=", length(values)), "triangle=upper",
      "diagonal=implicit_1", paste0("row_ptr_file=", paths[[1L]]),
      paste0("col_idx_file=", paths[[2L]]),
      paste0("values_file=", paths[[3L]]), "row_ptr_type=uint64",
      "col_idx_type=uint32", "values_type=float32", "index_base=0",
      "value=r", "cross_block_edges=0"), paths[[4L]])
  }
  study07_assert(all(file.exists(paths)),
    "Failed to construct the deterministic Study 07 block CSR input.")
  list(prefix = prefix, paths = paths,
    sha256 = vapply(paths, study07_sha256, character(1)))
}

study07_load_inputs <- function(spec) {
  paths <- unlist(spec$source$paths, use.names = TRUE)
  study07_assert(all(file.exists(paths)),
    paste("A frozen Study 06 input is missing:",
      paste(paths[!file.exists(paths)], collapse = ", ")))
  observed <- vapply(paths, study07_sha256, character(1))
  expected <- spec$source$sha256[names(paths)]
  study07_assert(identical(unname(observed), unname(expected)),
    "A frozen Study 06 input hash changed.")

  truth <- readRDS(paths[["truth"]])
  reference <- readRDS(paths[["learned_block"]])
  glist <- readRDS(paths[["human_glist"]])
  ld_glist <- readRDS(paths[["training_ld"]])
  rows <- as.integer(ld_glist$sparseLD$rows)
  marker_ids <- as.character(truth$marker_ids)
  study07_assert(length(rows) == spec$source$training_count,
    "Study 06 training-row count changed.")
  full_ld_names <- as.character(ld_glist$sparseLD$marker_names)
  selected_index <- match(marker_ids, full_ld_names)
  study07_assert(!anyNA(selected_index) && !anyDuplicated(selected_index),
    "Study 06 selected markers are absent or duplicated in frozen LD provenance.")
  study07_assert(identical(rownames(truth$annotations), marker_ids),
    "Study 06 annotation marker order changed.")
  study07_assert(identical(colnames(truth$annotations),
    spec$model$annotation_columns), "Study 06 annotation columns changed.")

  y <- matrix(as.numeric(truth$training_y), ncol = 1L,
    dimnames = list(glist$ids[rows], "trait1"))
  working_glist <- benchmark_set_training_af(glist, 1L, marker_ids,
    as.numeric(truth$gwas$freq))
  stats <- sblr::make_summary_stats(Glist = working_glist, y = y, chr = 1L,
    rows = rows, scale = TRUE, nthreads = 1L)
  study07_assert(identical(stats$marker_names, marker_ids),
    "Reconstructed summary-statistic marker order changed.")
  study07_assert(identical(as.integer(stats$n), spec$source$training_count),
    "Reconstructed summary-statistic sample count changed.")
  native_reference <- reference$result$native_fit
  score_error <- max(abs(as.numeric(stats$wy[[1L]]) -
    as.numeric(native_reference$wy[, 1L])))
  study07_assert(is.finite(score_error) && score_error < 1e-5,
    "Reconstructed summary-statistic scores differ from Study 06.")

  study06 <- new.env(parent = environment())
  sys.source(file.path(.study07_root, "studies", "06_annotation_models",
    "spec.R"), envir = study06)
  sys.source(file.path(.study07_root, "studies", "06_annotation_models",
    "annotation-design.R"), envir = study06)
  alpha_truth <- study06$construct_annotation_truth(truth$annotations,
    study06$spec)$informative_annotations
  prior_truth <- study06$annotation_marker_probabilities(truth$annotations,
    alpha_truth, spec$model$gamma)
  recorded_prior <- as.matrix(truth$marker_truth[paste0(
    "true_prior_component_", 0:3)])
  study07_assert(max(abs(prior_truth - recorded_prior)) < 1e-10,
    "Study 06 alpha-to-prior truth identity changed.")

  input <- native_reference$input
  study07_assert(max(abs(as.numeric(input$B) - spec$model$initial_B)) < 1e-12 &&
    max(abs(as.numeric(input$E) - spec$model$initial_E)) < 1e-12,
    "Registered Study 06 B/E initialization changed.")
  intercept_spec <- input$annotation_intercept_prior
  intercept_native <- rbind(
    type = rep(0, ncol(alpha_truth)),
    mean = as.numeric(intercept_spec$mean),
    precision = as.numeric(intercept_spec$precision))
  colnames(intercept_native) <- colnames(alpha_truth)

  raw <- truth$marker_truth$raw_effect
  active <- truth$marker_truth$true_nonnull & is.finite(raw) & raw != 0
  effect_scale <- stats::median(truth$effects[active] / raw[active])
  residual <- truth$training_y - as.numeric(truth$training_x %*% truth$effects)
  truth_variance <- c(B = effect_scale^2,
    E = stats::var(residual),
    Vg = stats::var(as.numeric(truth$training_x %*% truth$effects)))
  truth_variance <- c(truth_variance,
    h2 = truth_variance[["Vg"]] /
      (truth_variance[["Vg"]] + truth_variance[["E"]]))

  csr <- study07_write_block_csr(truth,
    file.path(spec$output$local_dir, "inputs"))
  list(
    truth = truth,
    marker_truth = truth$marker_truth,
    alpha_truth = alpha_truth,
    prior_truth = prior_truth,
    variance_truth = truth_variance,
    stats = stats,
    Glist = working_glist,
    training_rows = rows,
    block_start = which(!duplicated(truth$block)),
    csr_prefix = csr$prefix,
    csr_hashes = csr$sha256,
    alpha_init = input$alpha_init,
    intercept_spec = intercept_spec,
    intercept_native = intercept_native,
    source_hashes = observed,
    score_crosscheck_max_abs = score_error)
}

study07_alpha_starts <- function(inputs) {
  baseline <- as.matrix(inputs$alpha_init)
  positive <- negative <- mixed <- baseline
  positive[-1L, ] <- 0.5
  negative[-1L, ] <- -0.5
  mixed["enriched_binary", ] <- c(0.5, -0.5, 0.5)
  mixed["continuous_signal", ] <- c(-0.5, 0.5, -0.5)
  mixed["null_annotation", ] <- c(0.25, -0.25, 0.25)
  list(baseline = baseline, positive = positive, negative = negative,
    mixed = mixed)
}

study07_delta_starts <- function() list(
  baseline = c(0L, 0L, 0L),
  positive = c(1L, 1L, 1L),
  negative = c(1L, 1L, 1L),
  mixed = c(1L, 0L, 1L))

study07_validate_design <- function(spec, inputs) {
  provenance <- benchmark_resolve_package_provenance("sblr", "../sblr",
    installation_command = "R CMD INSTALL ../sblr --no-multiarch")
  benchmark_assert_experiment_package(provenance, spec$packages$sblr$sha)
  study07_assert(identical(provenance$version, spec$packages$sblr$version),
    "Installed sblr version differs from the frozen Study 07 design.")
  study07_assert(identical(provenance$installed_tree_sha256,
    spec$packages$sblr$installed_tree_sha256),
    "Installed sblr tree differs from the frozen Study 07 design.")
  study07_assert(identical(spec$scope$simulated_data_replicates, 1L),
    "Study 07 must contain exactly one simulated-data replicate.")
  study07_assert(identical(spec$methods, c("SBayesR", "SBayesRC",
    "SBayesRC-EM", "SBayesRC-S", "SBayesRC-S-EM")),
    "Study 07 method registry changed.")
  study07_assert(file.exists(paste0(inputs$csr_prefix, ".meta.txt")),
    "The frozen Study 06 CSR prefix is unavailable.")
  identity_spec <- spec
  identity_spec$status <- "design_frozen_execution_pending"
  identity_spec$design_amendments <- NULL
  identity_spec$result <- NULL
  identity_spec$model$selection$joint_delta_init <- NULL
  list(provenance = provenance,
    design_hash = benchmark_hash_object(identity_spec),
    input_hash = benchmark_hash_object(list(
      source_hashes = inputs$source_hashes,
      marker_ids = inputs$truth$marker_ids,
      annotations = inputs$truth$annotations,
      effects = inputs$truth$effects,
      block = inputs$truth$block,
      score = inputs$stats$wy,
      rows = inputs$training_rows,
      selection_csr = inputs$csr_hashes)))
}
