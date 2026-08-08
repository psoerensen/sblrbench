study07_save_plot <- function(plot, path, width = 9, height = 6) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  ggplot2::ggsave(path, plot, width = width, height = height, dpi = 180,
    bg = "white")
  path
}

study07_make_figures <- function(summary, inputs, spec) {
  figures <- file.path(spec$output$local_dir, "figures")
  joint <- summary$payload$SBayesRC$result
  descriptor <- joint$convergence_traces$quantities
  keep <- which(descriptor$parameter_name == "alpha" &
    descriptor$stick_index == 1L & descriptor$annotation_name %in%
      c("enriched_binary", "continuous_signal"))
  values <- joint$convergence_traces$values[, , keep, drop = FALSE]
  retained <- seq.int(1L, dim(values)[1L], by = 20L)
  joint_trace <- do.call(rbind, lapply(seq_len(dim(values)[2L]), function(chain)
    do.call(rbind, lapply(seq_along(keep), function(k) data.frame(
      source = "joint retained chain", run = paste0("chain ", chain),
      step = retained, quantity = descriptor$annotation_name[keep[k]],
      value = values[retained, chain, k])))))
  em_trace <- do.call(rbind, lapply(names(summary$payload[["SBayesRC-EM"]]),
    function(start) {
      history <- summary$payload[["SBayesRC-EM"]][[start]]$result$mcem$history$alpha
      annotation_rows <- c(enriched_binary = 2L, continuous_signal = 3L)
      do.call(rbind, lapply(names(annotation_rows), function(annotation) data.frame(
        source = "MCEM outer trajectory", run = start,
        step = seq_len(dim(history)[3L]) - 1L,
        quantity = annotation,
        value = history[annotation_rows[[annotation]], 1L, ])))
    }))
  alpha_plot <- ggplot2::ggplot(rbind(joint_trace, em_trace),
      ggplot2::aes(step, value, colour = run, group = run)) +
    ggplot2::geom_line(linewidth = 0.45, alpha = 0.85) +
    ggplot2::facet_grid(quantity ~ source, scales = "free_x") +
    ggplot2::labs(x = "Retained draw (thinned) or outer iteration",
      y = "First-stick alpha", colour = NULL) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(legend.position = "bottom")

  joint_prior <- summary$joint_prior$pooled
  prior_records <- list(data.frame(method = "joint SBayesRC",
    truth = 1 - inputs$prior_truth[, 1L], estimate = 1 - joint_prior[, 1L]))
  for (start in names(summary$payload[["SBayesRC-EM"]])) {
    estimate <- summary$payload[["SBayesRC-EM"]][[start]]$result$mcem$component_prior
    prior_records[[length(prior_records) + 1L]] <- data.frame(
      method = paste("SBayesRC-EM", start),
      truth = 1 - inputs$prior_truth[, 1L], estimate = 1 - estimate[, 1L])
  }
  prior_plot <- ggplot2::ggplot(do.call(rbind, prior_records),
      ggplot2::aes(truth, estimate)) +
    ggplot2::geom_abline(slope = 1, intercept = 0, colour = "grey55") +
    ggplot2::geom_point(alpha = 0.25, size = 0.7, colour = "#0072B2") +
    ggplot2::facet_wrap(~method) +
    ggplot2::labs(x = "True active prior", y = "Estimated active prior") +
    ggplot2::theme_minimal(base_size = 11)

  truth_beta <- inputs$truth$effects
  beta_records <- list()
  for (method in c("SBayesR", "SBayesRC")) {
    estimate <- study07_genomic_values(summary$payload[[method]], method,
      inputs)$beta
    beta_records[[length(beta_records) + 1L]] <- data.frame(method = method,
      truth = truth_beta, estimate = estimate)
  }
  for (start in names(summary$payload[["SBayesRC-EM"]])) {
    estimate <- study07_genomic_values(
      summary$payload[["SBayesRC-EM"]][[start]], "SBayesRC-EM", inputs)$beta
    beta_records[[length(beta_records) + 1L]] <- data.frame(
      method = paste("SBayesRC-EM", start), truth = truth_beta,
      estimate = estimate)
  }
  beta_plot <- ggplot2::ggplot(do.call(rbind, beta_records),
      ggplot2::aes(truth, estimate)) +
    ggplot2::geom_abline(slope = 1, intercept = 0, colour = "grey55") +
    ggplot2::geom_point(alpha = 0.30, size = 0.7, colour = "#009E73") +
    ggplot2::facet_wrap(~method) +
    ggplot2::labs(x = "True beta", y = "Posterior mean beta") +
    ggplot2::theme_minimal(base_size = 11)

  selection <- summary$selection
  selection_long <- rbind(
    data.frame(annotation = selection$annotation, method = "joint PIP",
      start = "posterior", value = selection$joint_annotation_pip),
    do.call(rbind, lapply(names(summary$payload[["SBayesRC-S-EM"]]),
      function(start) data.frame(annotation = selection$annotation,
        method = "conditional EB-PIP", start = start,
        value = summary$payload[["SBayesRC-S-EM"]][[start]]$result$mcem$annotation_pip_eb))))
  truth_role <- data.frame(annotation = selection$annotation,
    value = selection$true_role)
  selection_plot <- ggplot2::ggplot(selection_long,
      ggplot2::aes(annotation, value, colour = start)) +
    ggplot2::geom_point(position = ggplot2::position_dodge(width = 0.45),
      size = 2.2) +
    ggplot2::geom_point(data = truth_role,
      ggplot2::aes(annotation, value), inherit.aes = FALSE,
      shape = 4, size = 3, stroke = 1.1, colour = "black") +
    ggplot2::facet_wrap(~method) +
    ggplot2::coord_cartesian(ylim = c(0, 1)) +
    ggplot2::labs(x = NULL, y = "Inclusion probability",
      colour = "Chain/start", caption = "Crosses show binary annotation truth.") +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 25, hjust = 1),
      legend.position = "bottom")

  c(
    annotation_stability = study07_save_plot(alpha_plot,
      file.path(figures, "annotation_stability.png"), 10, 6.5),
    active_prior_truth = study07_save_plot(prior_plot,
      file.path(figures, "active_prior_truth.png"), 10, 7),
    genomic_beta_truth = study07_save_plot(beta_plot,
      file.path(figures, "genomic_beta_truth.png"), 10, 7),
    annotation_selection = study07_save_plot(selection_plot,
      file.path(figures, "annotation_selection.png"), 9, 5.5))
}
