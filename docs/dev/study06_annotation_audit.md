# Study 06 annotation-informed models audit

## Audit boundary

This is a read-only scientific and implementation audit of the in-development
Study 06 material. No sampler, benchmark, or targets pipeline was run. No
capsule was promoted or regenerated. The pinned package inspected was `sblr`
0.2.0 at `02e8c74baa906e83c4a08d42a9cc6339b4e81072`; qgdata is pinned at
`6cca5819e711d326cfb2614d7e9d9f34942612cd`.

## Current state

The directory contains a complete old-style development graph rather than a
completed shared-framework study: `config.R`, deterministic annotation and
simulation code, method adapters, chain extraction, convergence diagnostics,
draft metrics, pilot/promotion code, a per-study targets graph, tests, and a
capsule-only convergence status page. The intended 40-fit benchmark never
started because both annotation-aware pilot fits failed the prespecified
maximum-history convergence gate.

The useful scientific content is real, but the orchestration is not an
authoritative long-term interface. The future migration should preserve the
scientific definitions in an ordinary specification and
`annotation-design.R`, while replacing the per-study target graph, CSV-only
checkpoint metadata, and duplicated adapters with the shared runner and
semantic checkpoints.

## Scientific question found

The apparent primary question is:

> When the marginal sparsity and effect-size mixture are held fixed, do
> correctly specified informative marker annotations improve marker
> prioritization and recovery under BayesRC/SBayesRC relative to homogeneous
> BayesR/SBayesR, and is that conclusion consistent across BED and CSR
> interfaces?

Secondary questions in the code are annotation-coefficient and implied-prior
recovery, effect and variance recovery, prediction, convergence, and runtime.
The code also drafts component classification, interaction contrasts, and CPO
summaries. Group priors, overlapping annotation models, learned-logistic
SBayesC, and multitrait annotation models are not part of the implemented
Study 06 design even though related package methods exist.

The current material does not clearly prioritize its many draft metrics. It
also mixes two ideas in the informative scenario--a binary enriched group and
a continuous signal--so the final report must state that the estimand is the
joint annotation design. Removing the continuous signal after observing that
its coefficient limited convergence would be a material post-pilot design
change and must not happen silently.

## Existing scientific design

### Data and preprocessing

- qgdata chromosome 1, all 37,991 markers passing the existing QC rules, and
  all 2,000 samples.
- Deterministic 70/30 train/test split with seed 3101 (1,400/600 samples).
- Marker and sample identities are explicitly aligned; training genotypes
  determine scaling and training-only summary statistics.
- CSR LD uses the existing 1,000-variant window, `r2_min = 0.001`, block size
  1,024, and one construction thread.

These mechanics should reuse `benchmark-data.R`; the annotation matrix must
remain marker-aligned after every filtering or reordering operation.

### Annotation design

The annotation matrix is deterministic at seed 5201 and contains, in order:

1. an explicit intercept;
2. `enriched_binary`, with exactly approximately 10% marker membership;
3. `continuous_signal`, a standardized Gaussian marker covariate;
4. `null_annotation`, an independently standardized Gaussian covariate.

The informative probit-stick coefficients use nonzero binary and continuous
signals and zero null coefficients. The intercept is numerically calibrated so
the expected number of nonnull markers is 50 with active-component weights
0.60/0.30/0.10 across mixture variances 0.01/0.1/1. The uninformative scenario
uses only an intercept calibrated to the same marginal mixture probabilities.
The expected informative nonnull share inside the enriched group is about
0.68668, versus about 0.10 under the uninformative design.

The current design has no categorical expansion, explicit overlap policy,
missing-value policy beyond rejecting non-finite values, or biological
annotation source. It is a controlled synthetic validation, not a real
functional-genomics benchmark. The construction is reproducible and compatible
with the pinned interfaces, but its exact marker IDs, column order,
standardization, and alpha matrices need frozen regression fixtures before a
final run.

### Simulation and coordinates

- Scenarios: `informative_annotations` and
  `uninformative_annotations`.
- Five replicates per scenario.
- Target heritability 0.30; expected 50 nonnull markers; mixture variances
  `0, 0.01, 0.1, 1`.
- Component allocation, Gaussian effect generation, training-scale
  heritability rescaling, and phenotype residual construction are deterministic.
- Simulation seeds use base 5301, scenario stride 100,000, replicate stride
  100, and offsets 11/23/37. Fit seeds use base 600,000 with scenario,
  replicate, and method offsets; four chain seeds are spaced by 101.
- Four methods yield 40 planned fits and 160 chains.

Generic seed construction, effect generation, phenotype construction, truth
tables, and oracle checks should reuse `benchmark-simulation.R`. Drawing
components from annotation-implied marker probabilities and recording the true
annotation probabilities remain annotation-specific extensions.

### Methods and controls

The intended comparison is exactly:

| Method | Interface | Role |
|---|---|---|
| `st_bed_bayesr` | `stblr_bed` | homogeneous BED baseline |
| `st_bed_bayesrc` | `stblr_bed` | annotation-aware BED model |
| `st_csr_sbayesr` | `stblr_csr` | homogeneous CSR baseline |
| `st_csr_sbayesrc` | `stblr_csr_annot` | annotation-aware CSR model |

Both annotation routes are implemented and tested in the pinned package. BED
uses full unscheduled sweeps; CSR uses the aligned sparse-LD input. The current
pilot requested four chains, true extended convergence traces, and the same
mixture multipliers. Other implemented CSR annotation policies are documented
in the support matrix but should not be added without changing the scientific
question.

## Fit-output and extraction support

Directly stored outputs include posterior mean marker effects, posterior mean
component probabilities, variance components, annotation coefficient means
and final states, annotation-scale means and final states, chain metadata,
runtime, and selected extended traces. BED also exposes an `annotation_prior`
field that is a final marker-prior state; CSR SBayesRC does not expose the same
field.

True retained alpha, sigma, and variance draws can be extracted when extended
traces are requested. Heritability and annotation-implied marker prior
probabilities can therefore be calculated draw by draw. Posterior intervals
must be based on those draws.

The current draft extraction is not yet scientifically uniform: it uses the
BED final marker-prior state when present, but transforms posterior-mean alpha
for CSR. Those are different quantities, and the nonlinear transform of a mean
is not the mean of the transformed posterior. Final Study 06 must use the same
draw-wise definition for both interfaces. It must not manufacture global
component-proportion traces or full marker-by-component posterior histories
that the fit does not store.

No package change is required to compute the narrow recommended estimands when
true alpha traces are retained. A nonbreaking package field such as
`annotation_prior_final`, consistently documented for BED and CSR, would make
final-state semantics clearer but is not a completion blocker.

## Convergence and feasibility

The stopped pilot used four chains and a 3,000-iteration maximum history, with
candidate burn-ins 250/500/750/1,000 and retained draws
500/1,000/1,500/2,000/2,500. Thresholds were rank-normalized R-hat <= 1.01,
bulk and tail ESS >= 400, and relative MCSE <= 0.05. Neither annotation method
passed. For both, the limiting quantity was
`alpha:continuous_signal:component_0_stick` at the selected 1,000 burn-in plus
2,000 retained draws.

The capsule records approximately 1,925 seconds for BED BayesRC and 586 seconds
for CSR SBayesRC at that maximum history. These values provide a planning
basis, not a runtime guarantee. The old local targets cache contains two native
fit objects with the same declared seeds but nonmatching trace values and
different runtimes; it is not safely reusable.

The ignored cache occupies about 214 MB. Its materialized scaled-genotype
object is about 131 MB and the two native fit objects are about 29--34 MB each.
At the recorded 3,000-iteration rates, simply summing ten fits of each of the
four methods gives roughly 12 serial fit-hours before qualification overhead;
the annotation fits will probably require longer histories, so that figure is
only a lower planning reference. Peak memory and parallel throughput were not
recorded well enough to promise a resource envelope.

Study 04 thresholds and scalar machinery are suitable, but Study 06 must add
all alpha and sigma coefficients plus draw-wise derived prior summaries to the
required convergence contract. The final benchmark must not start until a
prespecified annotation-specific qualification history passes.

## Results assessment

`results/reference/06_annotation_models/current-stop/` is a checksummed,
provenanced record of a valid stop decision at the pinned package SHA. It is
valid partial evidence, not a final method-performance capsule. It contains two
successful pilot fits' compact tables and traces, but no passing MCMC
recommendation and no 40-fit benchmark.

The ignored `results/local/study05_annotation_models/` tree is an old-number
development cache. Its native traces do not reproduce the stop capsule, so it
must remain retired rather than be translated into the semantic checkpoint
schema.

## Smallest authoritative metric set

Required:

- alpha recovery by annotation and stick: posterior mean, bias, RMSE, interval
  coverage, and interval width;
- draw-wise implied marker-prior recovery: RMSE/correlation, expected active
  count, and enriched-versus-unannotated contrast;
- marker nonnull PIP recovery: causal-marker PIP/rank, AUPRC, AUROC, and Brier
  score;
- posterior mean effect RMSE/correlation and genetic-value recovery;
- heritability, genetic-variance, and residual-variance recovery as model
  sanity checks;
- paired annotation-minus-baseline contrasts within interface and the
  informative-minus-uninformative interaction;
- convergence and runtime.

Optional: phenotype prediction summaries, component classification scores,
and sigmaSqAlpha behavior summaries (there is no generating truth for sigma).
BED-only CPO is not a fair cross-interface primary metric. Group-prior metrics,
credible sets, and multitrait metrics are outside this study.

## Primary gaps and recommendation

The blocking gap is convergence qualification of the continuous annotation
coefficient, followed by shared, semantically consistent annotation extraction.
The full grid, capsule, and report do not exist. The recommended design is to
retain the existing two scenarios, four methods, five replicates, annotations,
seeds, and controls; qualify the annotation samplers first; then run the 40-fit
grid only after the gate passes. This avoids redefining the study after seeing
the failed pilot.

If the longer prespecified qualification still fails, Study 06 remains blocked
and the next action belongs in `sblr`: diagnose and minimally correct annotation
coefficient mixing in both backends, then repeat package fixtures and the
qualification. A binary-only study is a possible redesigned alternative, but
it must be explicitly approved and registered as a new scientific design, not
substituted during implementation cleanup.

The ordered completion work and exact gates are in
`study06_annotation_completion_plan.md`; the proposed fit inventory is in
`study06_annotation_required_runs.csv`.
