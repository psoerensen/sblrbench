# Deterministic prediction simulation mechanics moved intact from Study 02.

benchmark_prediction_simulation_coordinates <- function(spec,
                                                         profile = "benchmark") {
  coordinates <- unique(benchmark_seeds(spec, profile)[c("scenario",
    "replicate", "architecture_seed", "simulation_seed")])
  rownames(coordinates) <- NULL
  coordinates
}

simulate_prediction_architecture <- function(coordinate, scaled_genotypes,
                                             spec) {
  validate_benchmark_spec(spec)
  scenario <- .benchmark_scalar_string(coordinate$scenario, "coordinate$scenario")
  if (!scenario %in% names(spec$scenarios))
    stop("Unknown prediction scenario: ", scenario, call. = FALSE)
  architecture <- spec$scenarios[[scenario]]
  marker_ids <- colnames(scaled_genotypes)
  sample_ids <- rownames(scaled_genotypes)
  if (is.null(marker_ids) || is.null(sample_ids))
    stop("scaled_genotypes must have canonical sample and marker names.",
      call. = FALSE)
  n_causal <- as.integer(spec$controls$simulation$n_causal)
  if (n_causal > length(marker_ids))
    stop("n_causal exceeds the marker count.", call. = FALSE)
  set.seed(as.integer(coordinate$simulation_seed))
  causal_index <- sort(sample.int(length(marker_ids), n_causal,
    replace = FALSE))
  causal_ids <- marker_ids[causal_index]
  if (identical(architecture$effect_distribution, "single_normal")) {
    component <- rep("single_normal", length(causal_index))
    raw_effect <- stats::rnorm(length(causal_index))
  } else if (identical(architecture$effect_distribution,
      "variance_mixture")) {
    component_index <- sample.int(length(architecture$mixture_var),
      length(causal_index), replace = TRUE, prob = architecture$mixture_prob)
    component <- paste0("variance_", architecture$mixture_var[component_index])
    raw_effect <- stats::rnorm(length(causal_index),
      sd = sqrt(architecture$mixture_var[component_index]))
  } else {
    stop("Unknown effect distribution for scenario `", scenario, "`.",
      call. = FALSE)
  }
  trait <- spec$data$trait
  effects <- matrix(0, nrow = length(marker_ids), ncol = 1L,
    dimnames = list(marker_ids, trait))
  effects[causal_index, 1L] <- raw_effect
  genetic_values <- scaled_genotypes %*% effects
  h2 <- spec$controls$simulation$h2
  target_vg <- h2 / (1 - h2)
  effect_scale <- sqrt(target_vg / stats::var(genetic_values[, 1L]))
  effects[, 1L] <- effects[, 1L] * effect_scale
  genetic_values <- scaled_genotypes %*% effects
  # Study 05 historically used an explicit residual stream. Other migrated
  # studies omit residual_offset and therefore retain their exact RNG order.
  if (!is.null(spec$seeds$residual_offset))
    set.seed(as.integer(coordinate$simulation_seed + spec$seeds$residual_offset))
  residual <- stats::rnorm(nrow(scaled_genotypes))
  residual <- residual - mean(residual)
  residual <- residual / stats::sd(residual)
  residuals <- matrix(residual, ncol = 1L,
    dimnames = list(sample_ids, trait))
  phenotypes <- genetic_values + residuals
  observed_h2 <- stats::var(genetic_values[, 1L]) /
    (stats::var(genetic_values[, 1L]) + stats::var(residuals[, 1L]))
  raw <- list(y = phenotypes, W = scaled_genotypes, B = effects,
    G = genetic_values, E = residuals,
    causal = list(shared = causal_ids,
      specific = stats::setNames(list(character()), trait), all = causal_ids),
    rsids = marker_ids, ids = sample_ids, h2_target = h2,
    h2_observed = observed_h2, shared_idx = causal_index,
    specific_idx = stats::setNames(list(integer()), trait),
    causal_rsids = causal_ids)
  simulation <- as_sblrbench_simulation(raw, study = spec$study,
    architecture = scenario, replicate = as.integer(coordinate$replicate),
    seed = as.integer(coordinate$simulation_seed))
  simulation$extras$effect_components <- data.frame(marker = causal_ids,
    component = component, raw_effect = raw_effect,
    final_effect = effects[causal_ids, 1L], stringsAsFactors = FALSE)
  simulation$extras$effect_distribution <- architecture$effect_distribution
  simulation$extras$effect_scale <- effect_scale
  validate_sblrbench_simulation(simulation)
  simulation
}

prepare_prediction_simulations <- function(spec, profile, data) {
  coordinates <- benchmark_prediction_simulation_coordinates(spec, profile)
  lapply(seq_len(nrow(coordinates)), function(i) {
    coordinate <- as.list(coordinates[i, , drop = FALSE])
    simulation <- simulate_prediction_architecture(coordinate,
      data$scaled$all, spec)
    simulation$data$train_ids <- data$split$train_ids
    simulation$data$test_ids <- data$split$test_ids
    validate_sblrbench_simulation(simulation)
    oracle <- check_oracle_genetic_values(simulation,
      tolerance = spec$validation$oracle_tolerance)
    test_simulation <- subset_sblrbench_simulation_samples(simulation,
      data$split$test_ids)
    stats <- benchmark_summary_stats(simulation, data$ld_glist, data$split,
      spec$data)
    list(coordinate = coordinate, simulation = simulation,
      test_simulation = test_simulation, oracle = oracle, stats = stats)
  })
}

parameter_estimand_truth <- function(simulation, spec) {
  effects <- simulation$truth$effects[, 1L]
  causal <- effects != 0
  vg <- stats::var(simulation$truth$genetic_values[, 1L])
  ve <- stats::var(simulation$truth$residuals[, 1L])
  values <- c(sum(causal)/length(effects), mean(effects[causal]^2),
    sum(effects^2), vg, ve, vg/(vg+ve))
  data.frame(study=spec$study,scenario=simulation$scenario$architecture,
    architecture=simulation$scenario$architecture,
    replicate=simulation$scenario$replicate,
    estimand_id=spec$estimands$estimand_id,truth=values,
    truth_type="realized_quantity",
    truth_definition=spec$estimands$truth_source,
    source="validated simulation on the fitted analysis scale",status="ok",
    reason="",stringsAsFactors=FALSE)
}

parameter_summary_stats <- function(simulation,glist,spec) {
  stats <- sblr::make_summary_stats(Glist=glist,
    y=simulation$truth$phenotypes,chr=spec$data$chromosome,
    rows=seq_len(nrow(simulation$truth$phenotypes)),scale=TRUE,nthreads=1L)
  if(!identical(stats$marker_names,simulation$data$marker_ids) ||
     !identical(stats$trait_names,spec$data$trait) ||
     !identical(as.integer(stats$n),nrow(simulation$truth$phenotypes)))
    stop("Parameter summary statistics are not aligned.",call.=FALSE)
  stats
}

prepare_parameter_estimation_simulations <- function(spec, profile, data) {
  coordinates <- benchmark_prediction_simulation_coordinates(spec, profile)
  lapply(seq_len(nrow(coordinates)), function(i) {
    coordinate <- as.list(coordinates[i,,drop=FALSE])
    simulation <- simulate_prediction_architecture(coordinate, data$scaled,
      spec)
    oracle <- check_oracle_genetic_values(simulation,
      tolerance=spec$validation$oracle_tolerance)
    if (!oracle$ok) stop("Parameter-estimation simulation oracle failed.",
      call.=FALSE)
    stats <- parameter_summary_stats(simulation,data$ld_glist,spec)
    list(coordinate=coordinate,simulation=simulation,oracle=oracle,stats=stats,
      truth=parameter_estimand_truth(simulation,spec))
  })
}

prepare_finemapping_simulations <- function(spec, profile, data) {
  logic <- new.env(parent = baseenv())
  sys.source(benchmark_spec_path(spec, spec$locus_design$implementation),
    envir = logic)
  coordinates <- unique(benchmark_seeds(spec, profile)[c("scenario",
    "replicate", "causal_seed", "simulation_seed")])
  lapply(seq_len(nrow(coordinates)), function(i) {
    coordinate <- as.list(coordinates[i, , drop = FALSE])
    design <- spec$causal_design
    selection <- logic$select_separated_causal_markers(data$working_glist,
      spec$data$chromosome, data$markers$marker_ids,
      spec$controls$simulation$n_causal, design$min_distance_bp,
      coordinate$causal_seed, design$min_maf, design$max_maf)
    raw <- sblr::mtsim(W=data$scaled,standardize_W=FALSE,nt=1L,
      n_shared=length(selection$marker_ids),n_specific=0L,
      causal_rsids=selection$marker_ids,h2=spec$controls$simulation$h2,
      seed=coordinate$simulation_seed)
    if(!identical(sort(raw$causal$all),sort(selection$marker_ids)))
      stop("mtsim() did not preserve the selected causal markers.",call.=FALSE)
    simulation <- as_sblrbench_simulation(raw,study=spec$study,
      architecture=coordinate$scenario,replicate=coordinate$replicate,
      seed=coordinate$simulation_seed)
    if(is.null(dim(simulation$truth$phenotypes)))
      simulation$truth$phenotypes <- matrix(simulation$truth$phenotypes,ncol=1L,
        dimnames=list(simulation$data$sample_ids,simulation$data$trait_names))
    simulation$extras$causal_selection <- selection
    validate_sblrbench_simulation(simulation)
    oracle <- check_oracle_genetic_values(simulation,
      tolerance=spec$validation$oracle_tolerance)
    stats <- sblr::make_summary_stats(Glist=data$ld_glist,
      y=simulation$truth$phenotypes,chr=spec$data$chromosome,scale=TRUE,
      nthreads=1L)
    list(coordinate=coordinate,simulation=simulation,selection=selection,
      loci=logic$study01_locus_table(selection,spec),oracle=oracle,stats=stats)
  })
}
