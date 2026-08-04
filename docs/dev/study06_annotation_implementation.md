# Study 06 annotation-model implementation

## Status

The software implementation is complete and the scientific study remains
`qualification_pending`. No qualification or performance fit was run during
implementation. The unchanged `current-stop` capsule remains partial evidence
for the earlier stop decision, not a final Study 06 capsule.

## Final specification and design

`studies/06_annotation_models/spec.R` fixes qgdata SHA `6cca5819...`, 2,000
chromosome-1 samples, 37,991 ordered QC markers, the deterministic 70/30 split,
two scenarios, five replicates, four matched methods, target heritability 0.30,
50 expected active markers, mixture variances 0/0.01/0.1/1, annotation order,
priors, controls, and seeds. The benchmark grid has 40 fits and 160 chains.

`annotation-design.R` owns the genuinely study-specific construction:
intercept, deterministic 10% enriched binary group, sample-z-scored continuous
signal, sample-z-scored null annotation, informative coefficients, and a
marginally matched intercept-only uninformative control. It validates marker
identity/order, finite values, group size, rank, standardization, enrichment,
and expected active-marker range. Generic loading, alignment, scaling,
phenotype simulation, seeds, and truth storage remain shared.

## Method mapping and extraction

The method layer resolves BED BayesR (`stblr_bed/bayesr`), BED BayesRC
(`stblr_bed/bayesrc`), CSR SBayesR (`stblr_csr/sbayesr`), and CSR SBayesRC
(`stblr_csr_annot/sbayesrc`, `annotation_probit_stick`). Baselines use the
frozen Study 04 BayesR settings. Qualification requests true selected
annotation traces. CSR summary statistics and sparse LD use training data only.

BED's final marker-prior state and a transformation of CSR posterior-mean alpha
are not compared. `extract_annotation_coefficient_traces()` accepts only true,
identified retained alpha and sigmaSqAlpha traces. For each chain and retained
iteration, `summarise_drawwise_annotation_prior()` rebuilds the complete alpha
matrix and applies `sblr::sbayesrc_marker_pi()` to the same ordered annotation
matrix. Marker component priors are then summarized across draws. Posterior
means, final states, input priors, and posterior component allocations remain
separate. Missing traces produce `status = "unavailable"`; no substitution is
made.

## Metrics and output contract

`R/metrics-annotation.R` contains stable alpha, annotation-prior, marker PIP,
rank, top-k, precision/recall, PIP calibration, effect recovery, and matched
method/scenario contrasts. Shared prediction and parameter metrics are reused
where their definitions match. Qualification outputs are kept separate from
final truth, prior, marker, parameter, prediction, convergence, and runtime
tables. Empty placeholders are not written.

## Qualification and checkpoints

The four entries are informative replicate 1 and uninformative replicate 1 for
BED BayesRC and CSR SBayesRC. Maximum history is 9,000 iterations; candidate
burn-ins are 1,000/2,000/3,000 and retained windows 2,000/4,000/6,000. Every
required quantity needs four chains, R-hat <= 1.01, bulk and tail ESS >= 400,
and relative MCSE <= 0.05. Failure at maximum history blocks final execution.

The decision JSON records spec/provenance hashes, four entries, available
history, selected windows, diagnostic extrema, decisions, and semantic history
identities. Timestamp is metadata only. Final mode rejects a missing, failed,
or stale decision. Historical caches are not translated.

Checkpoint identity includes study/task/mode, scenario, replicate, method,
ordered train/test samples and markers, annotation matrix and column order,
phenotype/truth, priors, controls, MCMC history, seeds, LD settings, sblr SHA,
qgdata SHA, and `sblrbench-semantic-v2`. Paths, filenames, working directories,
reports, and timestamps are excluded. Qualification and final identities are
distinct; history reuse requires exact compatibility and required traces.

## Validation and remaining execution

Focused tests cover spec failures, 4/40 grids, 160 unique chain seeds,
deterministic annotation invariants, comparable draw-wise priors, unavailable
trace handling, marker metrics, safe default modes, decision identity, and
capsule-only reporting. Validation-only execution makes no fit call.

Remaining scientific work is separate: run four qualification entries,
review/freeze the decision, and only after a complete pass explicitly run the
40-fit benchmark. A package-side blocker exists if either route ceases to
expose true retained alpha and sigmaSqAlpha traces with chain identity;
final-state or posterior-mean substitution is prohibited.
