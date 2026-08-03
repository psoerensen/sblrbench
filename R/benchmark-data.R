# Shared data mechanics established by the Study 02 migration. Supported input
# remains the pinned qgdata PLINK panel or an explicitly supplied compatible
# Glist; this is not a general data-ingestion layer.

benchmark_data_paths <- function(output_dir) {
  .benchmark_scalar_string(output_dir, "output_dir")
  list(
    glist_path = Sys.getenv("SBLR_BENCH_GLIST", ""),
    data_dir = file.path(output_dir, "checkpoints", "data"),
    ld_dir = file.path(output_dir, "checkpoints", "ld")
  )
}

benchmark_example_files <- function(data_dir, data_spec) {
  paths <- download_sblrbench_example_data(data_dir)
  expected <- data_spec$example_data
  if (!identical(attr(paths, "repository"), expected$repository) ||
      !identical(attr(paths, "commit"), expected$commit) ||
      !identical(attr(paths, "subdirectory"), expected$subdirectory) ||
      !identical(names(paths), expected$files))
    stop("Example-data files do not match the pinned benchmark specification.",
      call. = FALSE)
  paths
}

benchmark_load_glist <- function(paths, example_files = NULL) {
  if (!requireNamespace("qgg", quietly = TRUE))
    stop("Data preparation requires the suggested package `qgg`.",
      call. = FALSE)
  if (nzchar(paths$glist_path)) {
    if (!file.exists(paths$glist_path))
      stop("SBLR_BENCH_GLIST does not exist: ", paths$glist_path,
        call. = FALSE)
    return(readRDS(paths$glist_path))
  }
  cache <- file.path(paths$data_dir, "human_glist.rds")
  if (file.exists(cache)) return(readRDS(cache))
  glist <- qgg::gprep(study = "sblrbench Study 02 prediction",
    bedfiles = unname(example_files["human.bed"]),
    bimfiles = unname(example_files["human.bim"]),
    famfiles = unname(example_files["human.fam"]))
  benchmark_atomic_save_rds(glist, cache, compress = FALSE,
    temporary_prefix = ".glist-")
  glist
}

benchmark_selected_ids <- function(glist, sample_limit = NULL) {
  ids <- glist$ids
  if (!is.character(ids) || anyNA(ids) || anyDuplicated(ids))
    stop("Glist sample IDs must be unique and non-missing.", call. = FALSE)
  if (!is.null(sample_limit)) {
    sample_limit <- as.integer(sample_limit)
    if (length(sample_limit) != 1L || is.na(sample_limit) || sample_limit < 2L)
      stop("sample_limit must be NULL or an integer of at least two.",
        call. = FALSE)
    ids <- ids[seq_len(min(sample_limit, length(ids)))]
  }
  ids
}

benchmark_filter_markers <- function(glist, chromosome, qc, sparse_ld) {
  if (!requireNamespace("qgg", quietly = TRUE))
    stop("Marker filtering requires the suggested package `qgg`.",
      call. = FALSE)
  retained <- do.call(qgg::gfilter, c(list(Glist = glist), qc))
  chromosome <- as.integer(chromosome)
  chromosome_ids <- glist$rsids[[chromosome]]
  if (!is.character(chromosome_ids) || anyNA(chromosome_ids) ||
      anyDuplicated(chromosome_ids))
    stop("Chromosome marker IDs must be unique and non-missing.",
      call. = FALSE)
  marker_ids <- chromosome_ids[!is.na(match(chromosome_ids, retained))]
  if (!length(marker_ids)) stop("No markers remain after QC.", call. = FALSE)
  idx <- match(marker_ids, chromosome_ids)
  af <- glist$af[[chromosome]][idx]
  maf <- glist$maf[[chromosome]][idx]
  if (any(!is.finite(af)) || any(!is.finite(maf)) ||
      any(maf < qc$excludeMAF))
    stop("Retained marker frequencies do not satisfy the configured QC threshold.",
      call. = FALSE)
  if (sparse_ld$max_distance_variants > 0L && length(marker_ids) < 2L)
    stop("At least two markers are required for the configured LD window.",
      call. = FALSE)
  list(marker_ids = marker_ids, af = af, maf = maf,
    marker_count_before = length(chromosome_ids),
    marker_count_after = length(marker_ids),
    marker_order_id = paste(chromosome, length(marker_ids), marker_ids[[1L]],
      marker_ids[[length(marker_ids)]], sep = ":"))
}

benchmark_set_glist_marker_order <- function(glist, chromosome, marker_ids) {
  chromosome <- as.integer(chromosome)
  marker_ids <- .assert_ids(as.character(marker_ids), "marker_ids")
  if (is.null(glist$rsidsLD)) glist$rsidsLD <- vector("list", length(glist$rsids))
  if (length(glist$rsidsLD) < length(glist$rsids))
    length(glist$rsidsLD) <- length(glist$rsids)
  glist$rsidsLD[[chromosome]] <- marker_ids
  glist
}

benchmark_extract_raw_genotypes <- function(glist, chromosome, sample_ids,
                                             marker_ids) {
  if (!requireNamespace("qgg", quietly = TRUE))
    stop("Genotype extraction requires the suggested package `qgg`.",
      call. = FALSE)
  sample_ids <- .assert_ids(as.character(sample_ids), "sample_ids")
  marker_ids <- .assert_ids(as.character(marker_ids), "marker_ids")
  x <- qgg::getG(Glist = glist, chr = as.integer(chromosome), ids = sample_ids,
    rsids = marker_ids, impute = FALSE, scale = FALSE)
  if (!identical(rownames(x), sample_ids) ||
      !identical(colnames(x), marker_ids))
    stop("Raw genotype extraction lost canonical sample or marker order.",
      call. = FALSE)
  x
}

benchmark_extract_scaled_genotypes <- function(glist, chromosome, sample_ids,
                                                marker_ids) {
  if (!requireNamespace("qgg", quietly = TRUE))
    stop("Genotype extraction requires the suggested package `qgg`.",
      call. = FALSE)
  sample_ids <- .assert_ids(as.character(sample_ids), "sample_ids")
  marker_ids <- .assert_ids(as.character(marker_ids), "marker_ids")
  x <- qgg::getG(Glist = glist, chr = as.integer(chromosome), ids = sample_ids,
    rsids = marker_ids, impute = TRUE, scale = TRUE)
  if (!identical(rownames(x), sample_ids) ||
      !identical(colnames(x), marker_ids) || any(!is.finite(x)))
    stop("Scaled genotype extraction lost canonical order or produced invalid values.",
      call. = FALSE)
  x
}

benchmark_make_full_sample_ld <- function(glist, marker_info, data_spec,
                                          output_dir) {
  chromosome <- as.integer(data_spec$chromosome)
  working <- benchmark_set_glist_marker_order(glist, chromosome,
    marker_info$marker_ids)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  prefix <- file.path(output_dir, "ld_sparse_bed")
  cache <- paste0(prefix, "_glist.rds")
  if (file.exists(cache)) {
    cached <- readRDS(cache)
    ids <- cached$rsids[[chromosome]][cached$sparseLD$cls[[1L]]]
    if (identical(ids, marker_info$marker_ids) &&
        length(Sys.glob(paste0(cached$sparseLD$prefix, "*")))) return(cached)
  }
  working <- do.call(sblr::make_sparse_ld, c(list(Glist = working,
    out_prefix = prefix, chr = chromosome, pos_bp = NULL),
    data_spec$sparse_ld))
  ids <- working$rsids[[chromosome]][working$sparseLD$cls[[1L]]]
  if (!identical(ids, marker_info$marker_ids))
    stop("Sparse-LD marker order does not match canonical marker order.",
      call. = FALSE)
  benchmark_atomic_save_rds(working, cache, compress = FALSE,
    temporary_prefix = ".full-sample-ld-")
  working
}

prepare_parameter_estimation_data <- function(spec, output_dir) {
  validate_benchmark_spec(spec)
  paths <- benchmark_data_paths(output_dir)
  files <- if (nzchar(paths$glist_path)) NULL else
    benchmark_example_files(paths$data_dir, spec$data)
  glist <- benchmark_load_glist(paths, files)
  markers <- benchmark_filter_markers(glist, spec$data$chromosome,
    spec$markers$qc, spec$data$sparse_ld)
  sample_ids <- benchmark_selected_ids(glist, spec$data$sample_limit)
  working <- benchmark_set_glist_marker_order(glist, spec$data$chromosome,
    markers$marker_ids)
  scaled <- benchmark_extract_scaled_genotypes(working,
    spec$data$chromosome, sample_ids, markers$marker_ids)
  ld_glist <- benchmark_make_full_sample_ld(working, markers, spec$data,
    paths$ld_dir)
  list(paths=paths,files=files,glist=glist,markers=markers,
    sample_ids=sample_ids,scaled=scaled,working_glist=working,
    ld_glist=ld_glist)
}

#' Create a deterministic prediction train/test split
#'
#' Sampling determines membership, while original sample order is retained
#' within each returned set.
#' @param sample_ids Unique sample identifiers.
#' @param train_fraction Fraction assigned to training.
#' @param seed Deterministic sampling seed.
#' @export
make_prediction_split <- function(sample_ids, train_fraction = 0.7, seed = 1L) {
  sample_ids <- .assert_ids(as.character(sample_ids), "sample_ids")
  if (!is.numeric(train_fraction) || length(train_fraction) != 1L ||
      !is.finite(train_fraction) || train_fraction <= 0 || train_fraction >= 1)
    stop("train_fraction must be a finite scalar strictly between zero and one.",
      call. = FALSE)
  n_train <- floor(length(sample_ids) * train_fraction)
  if (n_train < 1L || n_train >= length(sample_ids))
    stop("The split must contain training and test samples.", call. = FALSE)
  if (length(seed) != 1L || is.na(seed))
    stop("seed must be a non-missing scalar.", call. = FALSE)
  set.seed(as.integer(seed))
  chosen <- sort(sample.int(length(sample_ids), n_train, replace = FALSE))
  test_rows <- setdiff(seq_along(sample_ids), chosen)
  list(train_ids = sample_ids[chosen], test_ids = sample_ids[test_rows],
    train_rows = chosen, test_rows = test_rows, split_seed = as.integer(seed),
    train_fraction = train_fraction,
    split_id = paste0("train", n_train, "-test", length(test_rows),
      "-seed", as.integer(seed)))
}

#' Learn and apply training-only genotype scaling
#'
#' @param raw_genotypes Numeric sample-by-marker dosage matrix.
#' @param train_rows Integer training row indices.
#' @return Training frequencies and consistently scaled matrices.
#' @export
training_scaled_genotypes <- function(raw_genotypes, train_rows) {
  if (!is.matrix(raw_genotypes) || !is.numeric(raw_genotypes) ||
      is.null(rownames(raw_genotypes)) || is.null(colnames(raw_genotypes)))
    stop("raw_genotypes must be a named numeric matrix.", call. = FALSE)
  .assert_ids(rownames(raw_genotypes), "genotype sample IDs")
  .assert_ids(colnames(raw_genotypes), "genotype marker IDs")
  if (any(!is.finite(raw_genotypes[!is.na(raw_genotypes)])))
    stop("Observed genotypes must be finite.", call. = FALSE)
  train_rows <- as.integer(train_rows)
  if (!length(train_rows) || anyNA(train_rows) || anyDuplicated(train_rows) ||
      any(train_rows < 1L | train_rows > nrow(raw_genotypes)))
    stop("train_rows must be unique valid row indices.", call. = FALSE)
  train <- raw_genotypes[train_rows, , drop = FALSE]
  means <- colMeans(train, na.rm = TRUE)
  if (any(!is.finite(means)))
    stop("Every marker must have an observed training genotype.", call. = FALSE)
  af <- means / 2
  scale <- sqrt(2 * af * (1 - af))
  if (any(!is.finite(af)) || any(af <= 0 | af >= 1) || any(scale <= 0))
    stop("Training allele frequencies must lie strictly between zero and one.",
      call. = FALSE)
  filled <- raw_genotypes
  missing <- which(is.na(filled), arr.ind = TRUE)
  if (nrow(missing)) filled[missing] <- means[missing[, 2L]]
  scaled <- sweep(sweep(filled, 2L, means, "-"), 2L, scale, "/")
  test_rows <- setdiff(seq_len(nrow(raw_genotypes)), train_rows)
  list(allele_frequency = stats::setNames(af, colnames(raw_genotypes)),
    center = stats::setNames(means, colnames(raw_genotypes)),
    scale = stats::setNames(scale, colnames(raw_genotypes)), all = scaled,
    train = scaled[train_rows, , drop = FALSE],
    test = scaled[test_rows, , drop = FALSE], train_rows = train_rows,
    test_rows = test_rows)
}

benchmark_set_training_af <- function(glist, chromosome, marker_ids, af) {
  marker_ids <- .assert_ids(as.character(marker_ids), "marker_ids")
  idx <- match(marker_ids, glist$rsids[[as.integer(chromosome)]])
  if (anyNA(idx) || length(af) != length(idx) || any(!is.finite(af)))
    stop("Training allele frequencies are not aligned.", call. = FALSE)
  out <- glist
  out$af[[as.integer(chromosome)]][idx] <- unname(af)
  out$maf[[as.integer(chromosome)]][idx] <- pmin(unname(af), 1 - unname(af))
  benchmark_set_glist_marker_order(out, chromosome, marker_ids)
}

benchmark_make_training_ld <- function(glist, split, marker_ids, data_spec,
                                        output_dir) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  prefix <- file.path(output_dir, paste0("training_ld_", split$split_id))
  cache <- paste0(prefix, "_glist.rds")
  chromosome <- as.integer(data_spec$chromosome)
  if (file.exists(cache)) {
    x <- readRDS(cache)
    ids <- x$rsids[[chromosome]][x$sparseLD$cls[[1L]]]
    if (identical(ids, marker_ids) &&
        identical(x$sparseLD$rows, split$train_rows) &&
        length(Sys.glob(paste0(x$sparseLD$prefix, "*")))) return(x)
  }
  x <- do.call(sblr::make_sparse_ld, c(list(Glist = glist,
    rows = split$train_rows, out_prefix = prefix, chr = chromosome),
    data_spec$sparse_ld))
  ids <- x$rsids[[chromosome]][x$sparseLD$cls[[1L]]]
  if (!identical(ids, marker_ids) ||
      !identical(x$sparseLD$rows, split$train_rows) ||
      !identical(as.integer(x$sparseLD$reference_n), length(split$train_rows)))
    stop("Training sparse-LD provenance is invalid.", call. = FALSE)
  benchmark_atomic_save_rds(x, cache, compress = FALSE,
    temporary_prefix = ".training-ld-")
  x
}

benchmark_summary_stats <- function(simulation, glist, split, data_spec) {
  chromosome <- as.integer(data_spec$chromosome)
  y <- simulation$truth$phenotypes[split$train_ids, , drop = FALSE]
  stats <- sblr::make_summary_stats(Glist = glist, y = y, chr = chromosome,
    rows = split$train_rows, scale = TRUE, nthreads = 1L)
  if (!identical(stats$marker_names, simulation$data$marker_ids) ||
      !identical(stats$trait_names, simulation$data$trait_names) ||
      !identical(as.integer(stats$n), length(split$train_ids)) ||
      !isTRUE(all.equal(unname(stats$af[[1L]]),
        unname(glist$sparseLD$af[[1L]]), tolerance = 0)))
    stop("Training summary statistics failed sample, marker, trait, or frequency checks.",
      call. = FALSE)
  stats
}

prepare_prediction_data <- function(spec, output_dir) {
  validate_benchmark_spec(spec)
  paths <- benchmark_data_paths(output_dir)
  files <- if (nzchar(paths$glist_path)) NULL else
    benchmark_example_files(paths$data_dir, spec$data)
  glist <- benchmark_load_glist(paths, files)
  markers <- benchmark_filter_markers(glist, spec$data$chromosome,
    spec$markers$qc, spec$data$sparse_ld)
  sample_ids <- benchmark_selected_ids(glist, spec$data$sample_limit)
  split <- make_prediction_split(sample_ids, spec$split$train_fraction,
    spec$split$seed)
  raw <- benchmark_extract_raw_genotypes(glist, spec$data$chromosome,
    sample_ids, markers$marker_ids)
  scaled <- training_scaled_genotypes(raw, split$train_rows)
  working <- benchmark_set_training_af(glist, spec$data$chromosome,
    markers$marker_ids, scaled$allele_frequency)
  ld_glist <- benchmark_make_training_ld(working, split, markers$marker_ids,
    spec$data, paths$ld_dir)
  list(paths = paths, files = files, glist = glist, markers = markers,
    sample_ids = sample_ids, split = split, scaled = scaled,
    working_glist = working, ld_glist = ld_glist)
}
