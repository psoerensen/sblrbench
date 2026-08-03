# Supplemental Study 06: SBayesR LD-operator robustness

The scheduler and exact-versus-sparse SBayesR diagnostics are consolidated as
additional Study 06 package-validation evidence:

[`studies/06_ld_operator/sbayesr_ld_robustness/report.qmd`](../../studies/06_ld_operator/sbayesr_ld_robustness/report.qmd)

The compact report capsule is under
`results/reference/06_ld_operator/sbayesr_ld_robustness/current/`. Large genotype inputs,
operators, spectra, and fit checkpoints remain ignored under `results/local/`.

The bounded result is that exact CSR reproduced BED within one deterministic
LD-rich 1,500-marker window, while hard sparsification changed marker scores,
quadratic and residual expressions, posterior variance, component allocation,
and recovery. The configured full-rank block operator faithfully reproduced its
block target but not global exact LD because it omitted cross-block entries.
Retaining 0.995 positive spectral mass gave modest improvement but did not
remove that block-boundary limitation.

All fits use or were validated against installed `sblr` SHA
`02e8c74baa906e83c4a08d42a9cc6339b4e81072`. Principal limitations are the
single deterministic window, one reduced phenotype, one fixed block policy,
and the inability to infer a universal LD construction from this diagnostic.
