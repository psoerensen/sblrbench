# Study 06 external-reference parity helpers for official zhilizheng/SBayesRC.
# This profile is diagnostic only and cannot update qualification decisions.

study06_gctb_profile <- function() list(
  id = "gctb_parity",
  output = file.path("results", "local", "06_annotation_models",
    "gctb_parity"),
  official_repository = "https://github.com/zhilizheng/SBayesRC",
  official_tag = "v0.2.6",
  official_sha = "b95d3fcbad8ff358290922a58fff893439296138",
  spec_hash = "241c15afab8fefc571e38e625130de6e4ab58b958c30cff369e26998ee30fa56",
  truth_hash = "169d52bef390022a9106d7e61b200493869b40cbb76f1c8d911eebbc80fea1eb",
  gamma_matched = c(0, .01, .1, 1),
  start_pi_matched = c(.88, .06, .036, .024),
  gamma_native = c(0, .001, .01, .1, 1),
  start_pi_native = c(.990, .005, .003, .001, .001),
  niter = 3000L, burn = 1000L, out_freq = 1L,
  smoke_niter = 40L, smoke_burn = 20L,
  tune_iter = 150L, tune_burn = 100L,
  tune_step = c(.995, .99, .95, .9),
  chain_seed_offsets = c(G0 = 10000L, G1 = 20000L, G2 = 30000L))

study06_gctb_registry <- function(chain_seeds, profile = study06_gctb_profile()) {
  chain_seeds <- as.integer(chain_seeds)
  if (length(chain_seeds) != 4L || anyNA(chain_seeds) ||
      anyDuplicated(chain_seeds))
    stop("Four unique Study 06 chain seeds are required.", call. = FALSE)
  condition <- rep(c("G0", "G1", "G2"), each = 4L)
  chain <- rep(seq_len(4L), 3L)
  offset <- unname(profile$chain_seed_offsets[condition])
  data.frame(
    fit_id = paste(condition, paste0("chain", chain), sep = "--"),
    condition = condition, chain = chain,
    seed_source = chain_seeds[chain],
    requested_seed = chain_seeds[chain] + offset,
    annotations = condition != "G0",
    gamma = vapply(condition, function(x) paste(if (x == "G2")
      profile$gamma_native else profile$gamma_matched, collapse = ";"),
      character(1)),
    start_pi = vapply(condition, function(x) paste(if (x == "G2")
      profile$start_pi_native else profile$start_pi_matched, collapse = ";"),
      character(1)),
    bTune = condition == "G2", niter = profile$niter,
    burn = profile$burn, out_freq = profile$out_freq,
    stringsAsFactors = FALSE)
}

study06_gctb_sha256 <- function(path) {
  digest::digest(file = path, algo = "sha256")
}

study06_gctb_chr_field <- function(x, chromosome) {
  if (is.list(x)) x[[as.integer(chromosome)]] else x
}

study06_gctb_marker_metadata <- function(data, spec) {
  chromosome <- as.integer(spec$data$chromosome)
  glist <- data$working_glist
  ids <- as.character(data$markers$marker_ids)
  source_ids <- as.character(study06_gctb_chr_field(glist$rsids, chromosome))
  index <- match(ids, source_ids)
  if (anyNA(index)) stop("Official export marker metadata are unaligned.",
    call. = FALSE)
  get_field <- function(name) {
    value <- study06_gctb_chr_field(glist[[name]], chromosome)
    value[index]
  }
  panel <- data$marker_panel
  panel_index <- match(ids, panel$marker_id)
  if (anyNA(panel_index)) stop("Registered block panel is unaligned.",
    call. = FALSE)
  data.frame(
    SNP = ids, Chrom = chromosome,
    PhysPos = as.numeric(panel$position_bp[panel_index]),
    A1 = as.character(get_field("a1")),
    A2 = as.character(get_field("a2")),
    freq = as.numeric(data$scaled$allele_frequency[ids]),
    Block = as.integer(panel$block_id[panel_index]),
    stringsAsFactors = FALSE)
}

study06_gctb_gwas <- function(x, y, stats, metadata) {
  x <- as.matrix(x)
  y <- as.numeric(y)
  n <- nrow(x)
  if (n != 1400L || ncol(x) != 1500L || length(y) != n)
    stop("Study 06 official GWAS dimensions are invalid.", call. = FALSE)
  xpx <- colSums(x * x)
  xpy <- as.numeric(crossprod(x, y))
  yy <- sum(y * y)
  b <- xpy / xpx
  rss <- pmax(0, yy - 2 * b * xpy + b * b * xpx)
  se <- sqrt((rss / (n - 2L)) / xpx)
  z <- b / se
  p <- 2 * stats::pt(abs(z), df = n - 2L, lower.tail = FALSE)
  if (!identical(colnames(x), metadata$SNP) ||
      !identical(as.character(stats$marker_names), metadata$SNP) ||
      max(abs(xpy - as.numeric(stats$wy[[1L]]))) > 1e-8 ||
      max(abs(xpx - as.numeric(stats$ww[[1L]]))) > 1e-8)
    stop("Official GWAS reconstruction disagrees with Study 06 stats.",
      call. = FALSE)
  data.frame(SNP = metadata$SNP, A1 = metadata$A1, A2 = metadata$A2,
    freq = metadata$freq, b = b, se = se, p = p, N = n,
    Block = metadata$Block, stringsAsFactors = FALSE)
}

study06_gctb_annotation_export <- function(annotations, marker_ids) {
  if (!identical(rownames(annotations), marker_ids) ||
      !all(annotations[, 1L] == 1))
    stop("Study 06 annotations are unaligned or lack an exact intercept.",
      call. = FALSE)
  data.frame(SNP = marker_ids, Intercept = 1,
    enriched_binary = annotations[, "enriched_binary"],
    continuous_signal = annotations[, "continuous_signal"],
    null_annotation = annotations[, "null_annotation"], check.names = FALSE,
    row.names = NULL)
}

study06_gctb_write_eigen <- function(r, path) {
  eig <- eigen(r, symmetric = TRUE)
  tolerance <- max(abs(eig$values)) * nrow(r) * .Machine$double.eps
  keep <- which(eig$values > tolerance)
  if (!length(keep)) stop("Official LD block has no positive eigenmodes.",
    call. = FALSE)
  lambda <- eig$values[keep]
  u <- eig$vectors[, keep, drop = FALSE]
  connection <- file(path, "wb")
  on.exit(close(connection), add = TRUE)
  writeBin(as.integer(nrow(r)), connection, size = 4L)
  writeBin(as.integer(length(lambda)), connection, size = 4L)
  writeBin(as.numeric(sum(lambda)), connection, size = 4L)
  writeBin(as.numeric(1), connection, size = 4L)
  writeBin(as.numeric(lambda), connection, size = 4L)
  writeBin(as.numeric(u), connection, size = 4L)
  list(rank = length(lambda), eigenvalues = lambda, vectors = u,
    tolerance = tolerance,
    reconstruction_error = max(abs(r - u %*% (lambda * t(u)))))
}

study06_gctb_export <- function(data, simulation, stats, annotations, spec,
                                output, official_library) {
  export <- file.path(output, "export")
  ld_dir <- file.path(export, "ld")
  dir.create(ld_dir, recursive = TRUE, showWarnings = FALSE)
  metadata <- study06_gctb_marker_metadata(data, spec)
  marker_ids <- metadata$SNP
  x <- data$scaled$train[, marker_ids, drop = FALSE]
  y <- simulation$truth$phenotypes[data$split$train_ids, 1L]
  gwas <- study06_gctb_gwas(x, y, stats, metadata)

  ma_path <- file.path(export, "study06_informative.ma")
  data.table::fwrite(gwas[c("SNP", "A1", "A2", "freq", "b", "se",
    "p", "N")], ma_path, sep = "\t", quote = FALSE)

  # Official v0.2.6 treats the second input column as the intercept slot,
  # overwrites it with ones, and reads requested annotations from column 3.
  annotation_export <- study06_gctb_annotation_export(annotations, marker_ids)
  annotation_path <- file.path(export, "study06_informative.annot")
  data.table::fwrite(annotation_export, annotation_path, sep = "\t",
    quote = FALSE)

  block_ids <- sort(unique(metadata$Block))
  block_rows <- vector("list", length(block_ids))
  block_audit <- vector("list", length(block_ids))
  start <- 0L
  for (i in seq_along(block_ids)) {
    index <- which(metadata$Block == block_ids[i])
    if (length(index) != 100L || !identical(index,
        seq.int(min(index), max(index))))
      stop("Official LD blocks are not contiguous 100-marker blocks.",
        call. = FALSE)
    r <- crossprod(x[, index, drop = FALSE]) / nrow(x)
    r <- (r + t(r)) / 2
    path <- file.path(ld_dir, paste0("block", block_ids[i], ".eigen.bin"))
    written <- study06_gctb_write_eigen(r, path)
    block_rows[[i]] <- data.frame(Block = block_ids[i], Chrom = 1L,
      StartSnpIdx = start, StartSnpID = marker_ids[min(index)],
      EndSnpIdx = start + length(index) - 1L,
      EndSnpID = marker_ids[max(index)], NumSnps = length(index))
    block_audit[[i]] <- data.frame(Block = block_ids[i], marker_count =
      length(index), rank = written$rank, diagonal_min = min(diag(r)),
      diagonal_max = max(diag(r)), symmetry_error = max(abs(r - t(r))),
      reconstruction_error = written$reconstruction_error,
      positive_eigen_min = min(written$eigenvalues),
      positive_eigen_max = max(written$eigenvalues),
      retained_mass = sum(written$eigenvalues) / sum(eigen(r,
        symmetric = TRUE, only.values = TRUE)$values[
          eigen(r, symmetric = TRUE, only.values = TRUE)$values > 0]),
      stringsAsFactors = FALSE)
    start <- start + length(index)
  }
  ldm <- do.call(rbind, block_rows)
  data.table::fwrite(ldm, file.path(ld_dir, "ldm.info"), sep = "\t",
    quote = FALSE)
  snp_info <- data.frame(Chrom = metadata$Chrom, ID = metadata$SNP,
    Index = seq_len(nrow(metadata)) - 1L, GenPos = 0,
    PhysPos = metadata$PhysPos, A1 = metadata$A1, A2 = metadata$A2,
    A1Freq = metadata$freq, N = 1400L, Block = metadata$Block)
  data.table::fwrite(snp_info, file.path(ld_dir, "snp.info"), sep = "\t",
    quote = FALSE)

  old <- .libPaths()
  on.exit(.libPaths(old), add = TRUE)
  .libPaths(c(official_library, old))
  reader_rows <- lapply(block_ids, function(i) {
    value <- SBayesRC::readEig(ld_dir, i)
    data.frame(Block = i, reader_m = value$m, reader_k = value$k,
      reader_threshold = value$thresh,
      reader_reconstruction_error = max(abs(
        crossprod(x[, metadata$Block == i, drop = FALSE]) / nrow(x) -
        value$U %*% (value$lambda * t(value$U)))))
  })
  block_audit <- merge(do.call(rbind, block_audit),
    do.call(rbind, reader_rows), by = "Block", sort = FALSE)
  if (any(block_audit$reader_m != 100L) ||
      any(block_audit$reader_k != block_audit$rank) ||
      any(block_audit$reader_threshold != 1) ||
      max(block_audit$reader_reconstruction_error) > 1e-5)
    stop("Official LD reader reconstruction failed.", call. = FALSE)
  data.table::fwrite(block_audit, file.path(export, "ld_audit.csv"))

  truth <- list(marker_ids = marker_ids,
    marker_truth = simulation$extras$marker_truth,
    effects = as.numeric(simulation$truth$effects[, 1L]),
    validation_x = data$scaled$test[, marker_ids, drop = FALSE],
    validation_genetic = as.numeric(simulation$truth$genetic_values[
      data$split$test_ids, 1L]),
    validation_phenotype = as.numeric(simulation$truth$phenotypes[
      data$split$test_ids, 1L]), annotations = annotations,
    gwas = gwas, training_x = x, training_y = as.numeric(y),
    block = metadata$Block)
  saveRDS(truth, file.path(export, "truth.rds"), compress = FALSE)

  files <- sort(list.files(export, recursive = TRUE, full.names = TRUE))
  files <- files[!grepl("manifest\\.json$", files)]
  relative_paths <- substring(normalizePath(files, winslash = "/"),
    nchar(normalizePath(export, winslash = "/")) + 2L)
  file_manifest <- data.frame(path = unname(relative_paths),
    bytes = file.info(files)$size,
    sha256 = unname(vapply(files, study06_gctb_sha256, character(1))),
    row.names = NULL)
  manifest <- list(schema = "sblrbench-study06-gctb-export-v1",
    profile = "gctb_parity", specification_hash =
      study06_gctb_profile()$spec_hash,
    truth_hash = study06_gctb_profile()$truth_hash,
    marker_order_hash = benchmark_hash_object(marker_ids),
    gwas_hash = benchmark_hash_object(gwas),
    annotation_hash = benchmark_hash_object(annotations),
    ld_block_hash = benchmark_hash_object(list(metadata$Block,
      block_audit)), official_package_sha =
      study06_gctb_profile()$official_sha,
    export_script_sha256 = study06_gctb_sha256(file.path(
      "studies", "06_annotation_models", "gctb-parity.R")),
    n = 1400L, marker_count = 1500L, block_count = 15L,
    annotations_include_explicit_official_intercept_slot = TRUE,
    files = file_manifest)
  jsonlite::write_json(manifest, file.path(export, "manifest.json"),
    pretty = TRUE, auto_unbox = TRUE, digits = NA)
  list(manifest = manifest, gwas = gwas, annotations = annotation_export,
    block_audit = block_audit, truth = truth)
}

study06_gctb_bayesian_fdr <- function(pip, causal, level) {
  ordering <- order(-pip, seq_along(pip))
  expected_fdp <- cumsum(1 - pip[ordering]) / seq_along(ordering)
  eligible <- which(expected_fdp <= level)
  selected <- if (length(eligible)) ordering[seq_len(max(eligible))] else
    integer()
  c(selected = length(selected), true_causal = sum(causal[selected]))
}
