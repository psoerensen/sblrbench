.study07_detect_memory <- function() {
  value <- tryCatch({
    x <- system2("powershell", c("-NoProfile", "-Command",
      "(Get-CimInstance Win32_ComputerSystem).TotalPhysicalMemory"),
      stdout = TRUE, stderr = FALSE)
    as.numeric(tail(x, 1L))
  }, error = function(e) NA_real_)
  if (!is.finite(value) || value <= 0) NA_real_ else value
}

.study07_runtime_row <- function(run, marker_count, sample_count,
                                 controls, physical_memory = NA_real_) {
  fit_bytes <- if (identical(run$status, "ok"))
    as.numeric(object.size(run$fit)) else NA_real_
  estimate <- if (identical(run$status, "ok"))
    run$fit$memory_estimate$execution_estimated_total_bytes else NA_real_
  if (!is.finite(estimate) && identical(run$status, "ok"))
    estimate <- run$fit$memory_estimate$estimated_total_bytes
  iterations <- controls$nit + controls$nburn
  data.frame(implementation = run$implementation$id,
    marker_count = marker_count, sample_count = sample_count,
    chain_count = controls$nchains, measured_iterations = iterations,
    elapsed_seconds = run$runtime,
    seconds_per_iteration = run$runtime / iterations,
    seconds_per_marker_iteration = run$runtime /
      (iterations * marker_count), fit_object_bytes = fit_bytes,
    estimated_memory_bytes = estimate,
    physical_memory_bytes = physical_memory,
    warnings = paste(run$warnings, collapse = " | "),
    status = run$status, error_message = run$error,
    stringsAsFactors = FALSE)
}

.study07_runtime_projection <- function(timing, config,
                                        assumed_iterations = 1250L) {
  one <- timing[timing$chain_count == 1L & timing$status == "ok", ]
  if (!nrow(one)) stop("No successful Study 07 timing fits.")
  one$projected_four_chain_seconds <- one$seconds_per_iteration *
    assumed_iterations * 4
  one$projected_10_fit_seconds <- one$projected_four_chain_seconds * 10
  memory_limit <- if (any(is.finite(one$physical_memory_bytes)))
    max(one$physical_memory_bytes, na.rm = TRUE) *
      config$runtime_limits$maximum_memory_fraction else
        config$runtime_limits$conservative_memory_limit_bytes
  one$memory_limit_bytes <- memory_limit
  one$per_fit_time_pass <- one$projected_four_chain_seconds <=
    config$runtime_limits$maximum_projected_four_chain_seconds
  one$total_time_pass <- TRUE
  one$memory_pass <- !is.finite(one$estimated_memory_bytes) |
    one$estimated_memory_bytes * 4 <= memory_limit
  one$fit_object_pass <- one$fit_object_bytes <=
    config$runtime_limits$maximum_fit_object_bytes
  one$feasible <- one$per_fit_time_pass & one$total_time_pass &
    one$memory_pass & one$fit_object_pass
  one
}

.study07_select_marker_count <- function(projection, config) {
  required <- config$runtime_implementations
  candidates <- sort(unique(projection$marker_count), decreasing = TRUE)
  for (m in candidates) {
    z <- projection[projection$marker_count == m, ]
    total <- sum(z$projected_10_fit_seconds[
      match(config$implementations, z$implementation)])
    if (all(required %in% z$implementation) && is.finite(total) &&
        total <= config$runtime_limits$maximum_projected_benchmark_seconds &&
        all(z$feasible[match(required, z$implementation)]))
      return(data.frame(marker_count = m,
        selection_status = "frozen_feasible",
        projected_30_fit_seconds = total,
        rationale = "largest common nested marker count within frozen per-fit, total-grid, memory, and object-size limits",
        stringsAsFactors = FALSE))
  }
  stop("No common marker count of at least 2,000 is feasible.", call. = FALSE)
}
