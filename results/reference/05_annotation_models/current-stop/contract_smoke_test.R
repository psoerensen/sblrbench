run_study05_contract_smoke_test <- function() {
  required <- c("stblr_bed", "stblr_csr", "stblr_csr_annot",
    "sbayesrc_marker_pi", "make_sbayesrc_alpha_init")
  exports <- getNamespaceExports("sblr")
  if (!all(required %in% exports))
    stop("Installed sblr lacks required Study 05 exports.", call. = FALSE)
  bed <- formals(sblr::stblr_bed)
  csr <- formals(sblr::stblr_csr_annot)
  if (!"bayesrc" %in% eval(bed$method) ||
      !"annotation_probit_stick" %in% eval(csr$annotation_model))
    stop("Installed Study 05 model strings differ from the audited contract.",
      call. = FALSE)
  internal_bed <- formals(get("mtblr_bed", asNamespace("sblr")))
  required_controls <- c("mixture_var", "annotations", "alpha_init",
    "sigmaSqAlpha_init", "sigmaSqAlpha_a", "sigmaSqAlpha_b",
    "alpha_update_every", "nit", "nburn", "nthin", "nchains",
    "ncores", "chain_seeds", "keep_chains", "convergence",
    "convergence_control")
  if (!all(required_controls %in% names(internal_bed)))
    stop("Installed BayesRC control contract is incomplete.", call. = FALSE)
  A <- cbind(Intercept = 1, signal = c(-1, 0, 1))
  rownames(A) <- paste0("m", 1:3)
  alpha <- matrix(c(-3, -.25, -.67, .5, 0, 0), 2L, 3L,
    byrow = TRUE, dimnames = list(colnames(A), paste0("step_", 1:3)))
  p <- sblr::sbayesrc_marker_pi(A, alpha, c(0, .01, .1, 1))
  if (any(!is.finite(p)) || any(p < 0 | p > 1) ||
      max(abs(rowSums(p) - 1)) > 1e-12)
    stop("Installed stick-breaking transform contract failed.", call. = FALSE)
  invisible(TRUE)
}

if (sys.nframe() == 0L) {
  run_study05_contract_smoke_test()
  cat("Study 05 sampler-free contract smoke test passed.\n")
}
