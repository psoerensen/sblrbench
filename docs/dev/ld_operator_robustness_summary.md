# Study 05 SBayesR LD-sensitivity evidence

The former supplemental Study 06 evidence is now an explicit component of the
single **Study 05 — LD-operator validation** contract in
[`spec.R`](../../studies/05_ld_operator/spec.R). Its scientific narrative is
integrated into the sole [`report.qmd`](../../studies/05_ld_operator/report.qmd),
and its compact tables use `sbayesr_` prefixes within the sole capsule at
`results/reference/05_ld_operator/current/`.

The nested source, report, website page, and capsule paths were retired. Unique
block, spectral, corrected-score, quadratic, and residual audit functions were
consolidated into
[`operator-design.R`](../../studies/05_ld_operator/operator-design.R). No fit
was rerun during integration.

The bounded result is unchanged: exact CSR reproduced full-sweep BED BayesR in
one deterministic LD-rich 1,500-marker window. Hard-sparse and block-diagonal
operators omitted weak or cross-block LD, changing corrected marker scores,
component probabilities, residual expressions, and posterior recovery. The
hard-sparse operator was positive definite, and retaining 99.5% of within-block
spectral mass gave only modest improvement.

The evidence does not establish a universal LD window, threshold, block rule,
or retained-rank policy, and it does not independently justify an `sblr`
default change. The validated `sblr` SHA remains
`02e8c74baa906e83c4a08d42a9cc6339b4e81072`.
