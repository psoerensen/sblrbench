# Deferred multi-trait Study 02 extension
#
# These specifications preserve the experimentally validated public interfaces
# without placing them in the active targets graph. They are deliberately not
# sourced by targets.R while config$multitrait$enabled is FALSE.
.study02_deferred_multitrait_specs <- function() {
  list(
    list(id = "mt_bed_bayesr", interface = "mtblr_bed", method = "bayesr",
      controls = list(residual_covariance = "diagonal")),
    list(id = "mt_csr_sbayesr", interface = "mtblr_csr", method = "sbayesr",
      controls = list(sample_overlap = "not_modeled"))
  )
}

.study02_deferred_multitrait_status <- function() {
  list(enabled = FALSE, status = "deferred",
    observed_development_runtime_seconds = c(mt_bed_bayesr = 3783.14,
      mt_csr_sbayesr = 2909.62),
    reason = "Excluded pending separate computational profiling and optimization.")
}
