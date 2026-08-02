.five_replicate_recommendation_path <- function() {
  configured <- Sys.getenv("SBLR_BENCH_RECOMMENDATIONS", "")
  if (nzchar(configured)) return(configured)
  file.path("results", "reference", "04_convergence", "current-selection",
    "method_recommendations.csv")
}

.current_refresh_local_root <- function() {
  configured <- Sys.getenv("SBLR_BENCH_REFRESH_ROOT", "")
  if (nzchar(configured)) return(configured)
  file.path("results", "local", "current_benchmark_refresh")
}

.five_replicate_recommendations <- function(path = .five_replicate_recommendation_path()) {
  if (!file.exists(path)) stop("Frozen Study 04 recommendations are unavailable.", call. = FALSE)
  x <- utils::read.csv(path, stringsAsFactors = FALSE)
  expected <- data.frame(
    method = c("st_bed_bayesc", "st_bed_bayesr", "st_csr_sbayesc", "st_csr_sbayesr"),
    recommended_nburn = c(250L, 250L, 250L, 250L),
    recommended_post_burnin_draws = c(250L, 2000L, 250L, 2000L),
    recommended_nit_argument = c(250L, 2000L, 250L, 2000L),
    recommended_nthin = rep(1L, 4L), recommended_nchains = rep(4L, 4L),
    recommended_ncores = rep(4L, 4L), stringsAsFactors = FALSE
  )
  cols <- names(expected)
  if (!all(cols %in% names(x)) || anyDuplicated(x$method) ||
      !setequal(x$method, expected$method))
    stop("Frozen Study 04 recommendation grid is invalid.", call. = FALSE)
  z <- x[match(expected$method, x$method), cols, drop = FALSE]
  rownames(z) <- NULL
  if (!isTRUE(all.equal(z, expected, check.attributes = FALSE)))
    stop("Frozen Study 04 recommendations do not match the validated settings.", call. = FALSE)
  x[match(expected$method, x$method), , drop = FALSE]
}

.five_replicate_mcmc <- function(method, recommendations = .five_replicate_recommendations()) {
  z <- recommendations[recommendations$method == method, , drop = FALSE]
  if (nrow(z) != 1L) stop("No unique recommendation for method: ", method, call. = FALSE)
  list(nit = as.integer(z$recommended_nit_argument),
    nburn = as.integer(z$recommended_nburn), nthin = as.integer(z$recommended_nthin),
    nchains = as.integer(z$recommended_nchains), ncores = as.integer(z$recommended_ncores),
    convergence = "core", keep_chains = TRUE,
    convergence_control = list(warn = FALSE, keep_traces = TRUE))
}

.five_replicate_chain_seeds <- function(fit_seed, nchains = 4L) {
  seeds <- as.integer(fit_seed + seq_len(nchains) * 100000L)
  if (anyDuplicated(seeds)) stop("Chain seeds must be unique.", call. = FALSE)
  seeds
}

.five_replicate_sblr_provenance <- function() {
  d <- utils::packageDescription("sblr")
  commit <- d$RemoteSha %||% d$GithubSHA1 %||% NA_character_
  list(version = as.character(utils::packageVersion("sblr")), commit = commit,
    source_status = if (is.na(commit)) "installed package; source commit unavailable" else
      "installed package; commit from package metadata")
}

.five_replicate_thread_settings <- function() {
  vars <- c(OMP_NUM_THREADS = "4", OMP_THREAD_LIMIT = "4", OPENBLAS_NUM_THREADS = "1",
    MKL_NUM_THREADS = "1", VECLIB_MAXIMUM_THREADS = "1")
  observed <- Sys.getenv(names(vars), unset = NA_character_)
  if (!identical(unname(observed), unname(vars)))
    stop("Conservative numerical-library thread settings are not active.", call. = FALSE)
  as.list(observed)
}

`%||%` <- function(x, y) if (is.null(x) || !length(x) || is.na(x[[1L]]) || !nzchar(x[[1L]])) y else x
