.study01_require <- function(package) {
  if (!requireNamespace(package, quietly = TRUE)) {
    stop("Study 01 requires the optional package '", package,
         "'. Install the declared Suggests dependency before running it.",
         call. = FALSE)
  }
}

.study01_paths <- function() {
  glist_path <- Sys.getenv("SBLR_BENCH_GLIST", "")
  data_dir <- Sys.getenv(
    "SBLR_BENCH_DATA_DIR",
    file.path("results", "local", "01_finemapping", "data")
  )
  output_dir <- file.path("results", "local", "01_finemapping",
                          "genotype_setup")
  list(glist_path = glist_path, data_dir = data_dir,
       output_dir = output_dir,
       source = if (nzchar(glist_path)) "existing_glist" else "qgg_example")
}

.study01_example_files <- function(data_dir, example_data) {
  paths <- sblrbench::download_sblrbench_example_data(data_dir)
  if (!identical(attr(paths, "repository"), example_data$repository) ||
      !identical(attr(paths, "commit"), example_data$commit) ||
      !identical(attr(paths, "subdirectory"), example_data$subdirectory) ||
      !identical(names(paths), example_data$files)) {
    stop("Study 01 example-data configuration does not match the pinned package manifest.", call. = FALSE)
  }
  paths
}

.study01_load_glist <- function(paths, example_files = NULL) {
  .study01_require("qgg")
  if (nzchar(paths$glist_path)) {
    if (!file.exists(paths$glist_path)) {
      stop("SBLR_BENCH_GLIST does not exist: ", paths$glist_path,
           call. = FALSE)
    }
    return(readRDS(paths$glist_path))
  }
  cache <- file.path(paths$data_dir, "human_glist.rds")
  if (file.exists(cache)) return(readRDS(cache))
  Glist <- qgg::gprep(
    study = "Example",
    bedfiles = unname(example_files["human.bed"]),
    bimfiles = unname(example_files["human.bim"]),
    famfiles = unname(example_files["human.fam"])
  )
  saveRDS(Glist, cache)
  Glist
}

.study01_filter_markers <- function(Glist, chr, retained_ids, maf_threshold,
                                    max_distance_variants) {
  chr <- as.integer(chr)
  chromosome_ids <- Glist$rsids[[chr]]
  if (is.null(chromosome_ids)) stop("Chromosome is absent from Glist$rsids.", call. = FALSE)
  if (!is.character(chromosome_ids) || anyNA(chromosome_ids) ||
      anyDuplicated(chromosome_ids)) {
    stop("Chromosome marker IDs must be unique and non-missing.", call. = FALSE)
  }
  if (!is.character(retained_ids) || anyNA(retained_ids) ||
      anyDuplicated(retained_ids)) {
    stop("Retained marker IDs must be unique and non-missing.", call. = FALSE)
  }
  absent <- retained_ids[is.na(match(retained_ids, unlist(Glist$rsids,
                                                          use.names = FALSE)))]
  if (length(absent)) {
    stop("Retained marker IDs are absent from Glist: ",
         paste(absent, collapse = ", "), call. = FALSE)
  }
  marker_ids <- chromosome_ids[!is.na(match(chromosome_ids, retained_ids))]
  if (!length(marker_ids)) stop("No markers remain on the selected chromosome.", call. = FALSE)
  if (max_distance_variants > 0L && length(marker_ids) < 2L) {
    stop("At least two markers are required for the configured LD window.",
         call. = FALSE)
  }
  idx <- match(marker_ids, chromosome_ids)
  af <- Glist$af[[chr]][idx]
  maf <- Glist$maf[[chr]][idx]
  if (length(af) != length(marker_ids) || any(!is.finite(af))) {
    stop("Aligned allele frequencies must be finite.", call. = FALSE)
  }
  if (length(maf) != length(marker_ids) || any(!is.finite(maf)) ||
      any(maf < maf_threshold)) {
    stop("Retained MAF values do not satisfy the configured QC threshold.",
         call. = FALSE)
  }
  list(marker_ids = marker_ids, af = af, maf = maf,
       marker_count_before = length(chromosome_ids),
       marker_count_after = length(marker_ids),
       marker_order_id = paste(chr, length(marker_ids), marker_ids[1L],
                               marker_ids[length(marker_ids)], sep = ":"))
}

.study01_set_rsids_ld <- function(Glist, chr, marker_ids) {
  chr <- as.integer(chr)
  if (is.null(Glist$rsidsLD)) {
    Glist$rsidsLD <- vector("list", length(Glist$rsids))
  } else if (length(Glist$rsidsLD) < length(Glist$rsids)) {
    length(Glist$rsidsLD) <- length(Glist$rsids)
  }
  Glist$rsidsLD[[chr]] <- marker_ids
  Glist
}

.study01_run_qc <- function(Glist, config) {
  .study01_require("qgg")
  retained <- do.call(qgg::gfilter, c(list(Glist = Glist), config$qc))
  .study01_filter_markers(
    Glist, config$chr, retained, config$qc$excludeMAF,
    config$sparse_ld$max_distance_variants
  )
}

.study01_make_sparse_ld <- function(Glist, marker_info, config, output_dir) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  working <- .study01_set_rsids_ld(Glist, config$chr,
                                   marker_info$marker_ids)
  prefix <- file.path(output_dir, "ld_sparse_bed")
  cache <- paste0(prefix, "_glist.rds")
  if (file.exists(cache)) {
    cached <- readRDS(cache)
    cached_ids <- cached$rsids[[config$chr]][cached$sparseLD$cls[[1L]]]
    if (identical(cached_ids, marker_info$marker_ids) &&
        length(Sys.glob(paste0(cached$sparseLD$prefix, "*")))) return(cached)
  }
  working <- do.call(sblr::make_sparse_ld, c(
    list(Glist = working, out_prefix = prefix, chr = config$chr,
         pos_bp = NULL), config$sparse_ld
  ))
  if (is.null(working$sparseLD)) stop("sblr::make_sparse_ld() did not return sparseLD metadata.", call. = FALSE)
  sparse_ids <- working$rsids[[config$chr]][working$sparseLD$cls[[1L]]]
  if (!identical(sparse_ids, marker_info$marker_ids)) {
    stop("Sparse-LD marker order does not match the canonical marker order.", call. = FALSE)
  }
  saveRDS(working, cache)
  working
}

.study01_extract_genotypes <- function(Glist, chr, selected_ids, marker_ids) {
  .study01_require("qgg")
  Z <- qgg::getG(Glist = Glist, chr = chr, ids = selected_ids,
                 rsids = marker_ids, impute = TRUE, scale = TRUE)
  if (!identical(rownames(Z), selected_ids)) stop("qgg::getG() did not preserve selected sample order.", call. = FALSE)
  if (!identical(colnames(Z), marker_ids)) stop("qgg::getG() did not preserve canonical marker order.", call. = FALSE)
  if (anyDuplicated(rownames(Z)) || anyDuplicated(colnames(Z)) ||
      any(!is.finite(Z)) ||
      !identical(dim(Z), c(length(selected_ids), length(marker_ids)))) {
    stop("qgg::getG() returned an invalid genotype matrix.", call. = FALSE)
  }
  Z
}

.study01_selected_ids <- function(Glist, sample_limit = NULL) {
  ids <- Glist$ids
  if (!is.character(ids) || anyNA(ids) || anyDuplicated(ids)) {
    stop("Glist sample IDs must be unique and non-missing.", call. = FALSE)
  }
  if (!is.null(sample_limit)) ids <- ids[seq_len(min(as.integer(sample_limit), length(ids)))]
  ids
}

.study01_summary <- function(config, paths, marker_info, Z, simulation,
                             oracle) {
  list(
    study = config$study, task = config$task, data_source = paths$source,
    chromosome = config$chr, sample_count = nrow(Z),
    marker_count_before_qc = marker_info$marker_count_before,
    marker_count_after_qc = marker_info$marker_count_after,
    marker_order_id = marker_info$marker_order_id,
    qc = config$qc, sparse_ld = config$sparse_ld,
    simulation_seed = config$simulation$seed,
    number_of_causal_markers = length(simulation$truth$causal$all),
    heritability = config$simulation$h2, standardize_W = FALSE,
    genotype_extraction_function = "qgg::getG",
    impute = TRUE, scale = TRUE,
    oracle_tolerance = oracle$tolerance,
    oracle_maximum_absolute_error = oracle$max_abs_error,
    oracle_ok = oracle$ok, sblr_version = as.character(utils::packageVersion("sblr")),
    sblr_commit_inspected = "fd76b18828bc77756948aa3138de07ae4dc75513",
    sblrbench_commit = suppressWarnings(sblrbench::sblrbench_git_commit(".")),
    qgg_version = as.character(utils::packageVersion("qgg")),
    created_at = format(Sys.time(), tz = "UTC", usetz = TRUE)
  )
}
