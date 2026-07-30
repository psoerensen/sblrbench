# Study 04 multichain convergence contract

The four public fits are `stblr_bed(method = "bayesc"|"bayesr")` and
`stblr_csr(method = "sbayesc"|"sbayesr")`. Their shared chain-control parser
accepts `nchains`, `ncores`, an optional vector of exactly one signed 32-bit
`chain_seeds` value per chain, and `keep_chains`. Logical chains are independent;
workers are limited to `min(nchains, ncores)`. Resolved task seeds are recorded
in `fit$input$task_seeds_resolved`.

`nit` is the number of post-warm-up sampler draws. `nburn` is additional
warm-up. Native chain records may contain `nburn + nit` values, while the
retained convergence bundle always contains exactly `nit` values per chain.
`nthin` is a sampler/finalizer control and must not be applied again by Study 04.
For the maximum-history design Study 04 therefore sets `nburn = 0`, `nit =
3000`, and `nthin = 1`, then applies candidate burn-in windows once downstream.

With `convergence = "core"` and `convergence_control$keep_traces = TRUE`, the
fit retains `fit$convergence_traces$values` in proven dimensions iterations ×
chains × quantities, plus descriptors. Core quantities are `vbs`, `vgs`,
`ves`, `vle`, and `vld`. Study 04 uses the first three and derives heritability
draw by draw. The core bundle does not retain `pi_trace`; marker-level BayesR
component probabilities are posterior summaries, not draws, and are never
substituted. Chain-specific marker summaries are available only when
`keep_chains = TRUE`; they are secondary diagnostics.

`sblr` already provides rank-normalized diagnostics internally, but Study 04
recalculates R-hat, bulk/tail ESS, and MCSE with `posterior` for every requested
burn-in/checkpoint window. This preserves chain identity and supports the
recommendation algorithm. Accuracy, truth recovery, prediction, and universal
convergence claims are outside this study.
