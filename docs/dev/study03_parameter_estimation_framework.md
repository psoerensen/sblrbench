# Study 03 parameter-estimation framework

## Confirmed `sblr` parameter contract

The four supported public calls are `stblr_bed(method = "bayesc")`,
`stblr_bed(method = "bayesr")`, `stblr_csr(method = "sbayesc")`, and
`stblr_csr(method = "sbayesr")`. Their finalized single-trait fits retain trace
matrices named `vbs`, `vgs`, `ves`, `vle`, and `vld`; inclusion-trace
availability differs by finalizer as documented below. `fit$input`
records `nit`, `nburn`, `nthin`, `nchains`, marker count, model and scale
metadata. The public `summarise_posterior()` implementation confirms these
mean nonzero marker-effect variance, genetic variance, residual variance and
nonzero marker proportion, respectively. It derives `h2` as
`vgs / (vgs + ves)` and marker-effect variance as `vbs * pi_trace * m`, per
iteration.

Top-level trace matrices include burn-in. Study 03 reads `fit$input$nburn`,
removes it exactly once, checks a single trait, preserves joint row alignment,
and applies recorded thinning only when traces contain the unthinned scheduled
iterations. Development fits use `nthin = 1`, avoiding an ambiguous legacy
storage case.

BayesR `component_probabilities` and `dm_component_mean` are finalized
posterior summaries, not aligned draw-level class-proportion traces. They are
therefore not benchmarked as posterior parameters here. `vle` and `vld`
represent an LD decomposition (`vld = vgs - vle` in maintained tests), but a
defensible simulation truth on all supported BED/CSR operators is not yet
established. They are retained in native fits and excluded from primary
recovery. No scientific meaning was established for using these fields as
convergence diagnostics.

## Estimands and truths

The common primary estimands are realized effect variance (`vbs`), realized
genetic variance (`vgs`), realized residual variance (`ves`), and realized
heritability (`vgs / (vgs + ves)`). BayesC/SBayesC and the current CSR SBayesR
finalizer additionally retain
`pi_trace`, supporting method-specific causal proportion and total marker-effect
variance (`vbs * pi_trace * m`). The current BED BayesR finalizer does not retain an
aligned draw-level `pi_trace`; marker-level `component_probabilities` are
posterior summaries and are not substituted for it. Products and ratios are
calculated draw by draw before summarization.

The truth table explicitly distinguishes generating parameters from realized
quantities. Target heritability and the configured generating effect family are
generating parameters. Exact causal proportion, mean squared final causal
effect, sum of squared final effects, sample variance of genetic values, sample
variance of residuals, and their realized variance ratio are finite-sample
realized quantities. Sample variances use denominator `n - 1`, matching base R
`var()` and the analysis sample used by all fits.

## Design

Study 03 uses the pinned public qgdata genotype source, all 5,000 analysis
individuals, 37,991 canonical markers, 50 causal markers, and paired phenotypes
for the four single-trait methods. `sparse_homogeneous` is approximately
BayesC-like; `sparse_mixture` is approximately BayesR-like. Whole-vector
heritability rescaling means neither architecture is claimed to be an exact
draw from a fitted prior. Matching BayesC/SBayesC or BayesR/SBayesR results are
primary; cross-prior combinations are labelled misspecified robustness results.

The development profile uses one chain and short MCMC solely to validate the
contract and reporting path. The scientific profile is configuration only and
its lengths are provisional.

## Boundary with Study 04

Study 03 records posterior draws and computational metadata but does not
calculate R-hat, multichain effective sample size, between-chain diagnostics,
trace comparisons, convergence pass/fail rules, or select production MCMC
lengths. Those are reserved for Study 04.
