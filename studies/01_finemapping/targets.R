source(file.path("studies", "01_finemapping", "setup_example_data.R"),
       local = TRUE)

list(
  targets::tar_target(
    config,
    source(file.path("studies", "01_finemapping", "config.R"),
           local = TRUE)$value
  ),
  targets::tar_target(data_paths, .study01_paths()),
  targets::tar_target(
    example_files,
    if (identical(data_paths$source, "qgg_example")) {
      .study01_example_files(data_paths$data_dir)
    } else character()
  ),
  targets::tar_target(Glist, .study01_load_glist(data_paths, example_files)),
  targets::tar_target(filtered_markers, .study01_run_qc(Glist, config)),
  targets::tar_target(
    working_glist,
    .study01_set_rsids_ld(Glist, config$chr, filtered_markers$marker_ids)
  ),
  targets::tar_target(
    sparse_ld_glist,
    .study01_make_sparse_ld(working_glist, filtered_markers, config,
                            data_paths$output_dir)
  ),
  targets::tar_target(
    selected_ids,
    .study01_selected_ids(sparse_ld_glist, config$sample_limit)
  ),
  targets::tar_target(
    standardized_genotypes,
    .study01_extract_genotypes(
      sparse_ld_glist, config$chr, selected_ids,
      filtered_markers$marker_ids
    )
  ),
  targets::tar_target(
    raw_simulation,
    {
      required <- config$simulation$n_shared +
        config$simulation$nt * config$simulation$n_specific
      if (ncol(standardized_genotypes) < required) {
        stop("Too few QC-retained markers for the configured causal counts.")
      }
      tryCatch(
        sblr::mtsim(
          W = standardized_genotypes,
          standardize_W = FALSE,
          nt = config$simulation$nt,
          n_shared = config$simulation$n_shared,
          n_specific = config$simulation$n_specific,
          h2 = config$simulation$h2,
          seed = config$simulation$seed
        ),
        error = function(e) {
          stop(
            "Installed sblr could not execute the configured one-trait ",
            "mtsim() call: ", conditionMessage(e),
            ". No oracle summary was written.", call. = FALSE
          )
        }
      )
    }
  ),
  targets::tar_target(
    simulation,
    sblrbench::as_sblrbench_simulation(
      raw_simulation, study = config$study,
      architecture = "scale_validation", replicate = 0L,
      seed = config$simulation$seed, keep_genotypes = TRUE
    )
  ),
  targets::tar_target(
    oracle,
    sblrbench::check_oracle_genetic_values(
      simulation, tolerance = config$oracle_tolerance,
      stop_on_failure = TRUE
    )
  ),
  targets::tar_target(
    compact_summary,
    .study01_summary(config, data_paths, filtered_markers,
                     standardized_genotypes, simulation, oracle)
  ),
  targets::tar_target(
    summary_file,
    {
      dir.create(data_paths$output_dir, recursive = TRUE,
                 showWarnings = FALSE)
      path <- file.path(data_paths$output_dir, "genotype_setup_summary.json")
      jsonlite::write_json(compact_summary, path, pretty = TRUE,
                           auto_unbox = TRUE, null = "null", digits = NA)
      path
    },
    format = "file"
  )
)
