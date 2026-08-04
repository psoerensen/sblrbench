# Reporting harmonization audit

## Before harmonization

Study 01 used base-R boxplots, hard-coded colors, and bespoke formatting. Study
02 used ggplot defaults and `theme_minimal()`, and linked rather than displayed
its scripts. Study 03 had minimal front matter, default plot mappings, no
interpretation section, link-only reproduction material, and described Study 04
as future work. Study 04 displayed its worked example but combined smoke-test
and reproduction links and used different chain/estimand colors and formatting.

The catalogue was ordered 04, 03, 01, 02 and mixed wide and two-column entry
formats. The homepage described only Study 01. Terminology varied between
“runtime” and “computational summary”, “variance mixture” and other forms, and
BED/CSR comparison headings. Development notices used block quotes rather than
a common callout.

## Harmonized contract

Reports source the authoritative `R/benchmark-reporting.R`. The former
compatibility loader was removed in the final cleanup. The canonical method order is
BED BayesC, BED BayesR, CSR SBayesC, CSR SBayesR. Named Okabe–Ito-derived colors
are paired with circles/solid lines for BED and triangles/dashed lines for CSR.
Architectures are ordered sparse homogeneous then sparse variance mixture.

All reports use parallel title/description front matter, a development callout,
shared terminology and theme, readable `knitr::kable()` tables, structured
downloads/provenance, visible frozen worked examples and run scripts, and
collapsed frozen developer smoke tests. Replicate preparation retains every
observation, reports SD only with multiple replicates, and never fabricates
one-replicate uncertainty.

## Intentionally study-specific

Fine-mapping retains causal ranks and credible-set evaluation. Prediction
retains leakage controls, calibration and paired prior/representation contrasts.
Parameter estimation retains generating-versus-realized truth and posterior
intervals. Convergence retains chain traces, running means, diagnostic
trajectories, burn-in stability, recommendations and marker agreement. These
scientifically distinct outputs are not forced into one table or axis scale.
