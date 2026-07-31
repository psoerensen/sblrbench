.study04_method_spec <- function(id) {
  map <- list(st_bed_bayesc = c("stblr_bed", "bayesc"), st_bed_bayesr = c("stblr_bed", "bayesr"),
    st_csr_sbayesc = c("stblr_csr", "sbayesc"), st_csr_sbayesr = c("stblr_csr", "sbayesr"))
  if (!id %in% names(map)) stop("Unknown Study 04 method.", call. = FALSE)
  list(id = id, interface = map[[id]][1], native_method = map[[id]][2],
    representation = if (grepl("bed", id)) "BED" else "CSR",
    prior_class = if (grepl("bayesr", id)) "BayesR" else "BayesC")
}

.study04_chain_seeds <- function(architecture, method, config, replicate = 1L) {
  a <- match(architecture, names(config$simulation$architectures)); m <- match(method, config$methods)
  as.integer(config$seeds$base_fit + a * 1000000L + as.integer(replicate) * 10000L +
    m * 1000L + seq_len(4L) * config$seeds$chain_stride)
}

.study04_fit <- function(spec, simulation, stats, Glist, config) {
  method <- .study04_method_spec(spec$method)
  p <- if (identical(config$profile, "five_replicate_validation"))
    .five_replicate_mcmc(spec$method) else config$profiles$development
  seeds <- .study04_chain_seeds(spec$architecture, spec$method, config,
    spec$replicate %||% 1L)
  controls <- p; controls$seed <- seeds[1]; controls$chain_seeds <- seeds
  controls$verbose <- FALSE; controls$h2 <- config$priors$h2
  if (method$prior_class == "BayesR") {
    active <- config$priors$bayesr_active_probability
    controls$pi <- c(1 - active, rep(active / 3, 3)); controls$mixture_var <- config$priors$bayesr_mixture_var
  } else controls$pi_init <- config$priors$bayesc_inclusion_probability
  native <- sblrbench::new_sblr_native_method(method$id, method$id, method$interface,
    method$native_method, capabilities = c("scalar_trait", "multichain"))
  input <- if (method$representation == "BED") list(y = simulation$truth$phenotypes,
    Glist = Glist, rows = seq_len(nrow(simulation$truth$phenotypes))) else list(stats = stats, Glist = Glist)
  start <- proc.time()[["elapsed"]]
  tryCatch({
    result <- sblrbench::run_sblrbench_method(native, fit_inputs = input, controls = controls)
    list(status = "ok", reason = "", method = method, result = result,
      fit = result$native_fit, chain_seeds = seeds, runtime = proc.time()[["elapsed"]] - start)
  }, error = function(e) list(status = "failed", reason = conditionMessage(e), method = method,
    result = NULL, fit = NULL, chain_seeds = seeds, runtime = proc.time()[["elapsed"]] - start))
}

.study04_marker_agreement <- function(run, architecture) {
  if (run$status != "ok" || is.null(run$fit$chains)) return(data.frame(architecture = architecture,
    method = run$method$id, status = "unavailable", reason = "chain marker summaries not retained"))
  ch <- run$fit$chains; pairs <- utils::combn(seq_along(ch), 2)
  do.call(rbind, lapply(seq_len(ncol(pairs)), function(i) {
    a <- ch[[pairs[1, i]]]; b <- ch[[pairs[2, i]]]; ea <- as.numeric(a$bm); eb <- as.numeric(b$bm)
    pa <- as.numeric(a$dm); pb <- as.numeric(b$dm)
    data.frame(architecture = architecture, method = run$method$id, chain_a = pairs[1, i], chain_b = pairs[2, i],
      effect_correlation = cor(ea, eb), effect_rmse = sqrt(mean((ea-eb)^2)),
      maximum_absolute_effect_difference = max(abs(ea-eb)), pip_correlation = cor(pa,pb),
      mean_absolute_pip_difference = mean(abs(pa-pb)), top50_pip_overlap = length(intersect(order(pa,decreasing=TRUE)[1:50],order(pb,decreasing=TRUE)[1:50]))/50,
      status = "ok", reason = "")
  }))
}
