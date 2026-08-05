# Study 06 annotation-informed models completion plan

> **Execution update:** The shared-framework software described here has been
> implemented. The prespecified qualification ran on 2026-08-04 and failed;
> final fitting was not launched. See
> `study06_annotation_qualification_result.md`.
>
> **Versioning update:** the design below is retained as the historical v1
> sparse design and reclassified as `v1_sparse_stress`; it is not the current
> primary qualification. The current `v2_identifiable_qualification` contract
> is in `study06_v2_design.md`. Its four-entry qualification ran and failed;
> see `study06_v2_qualification_result.md`. It does not authorize the final
> benchmark.

## Recommended final design

### Required

- **Primary question:** whether correctly specified informative annotations
  improve marker prioritization and recovery under BayesRC/SBayesRC relative
  to marginally matched BayesR/SBayesR, consistently across BED and CSR.
- **Secondary questions:** recovery of annotation coefficients and implied
  marker priors; effect, genetic-value, and variance recovery; convergence and
  runtime.
- **Data:** the existing qgdata chromosome-1 design, 2,000 ordered samples,
  37,991 ordered QC markers, and deterministic 70/30 split.
- **Scenarios:** the existing informative and marginally matched uninformative
  annotation scenarios.
- **Annotations:** intercept, 10% enriched binary annotation, standardized
  continuous signal, and standardized null annotation in the existing order.
- **Replicates:** five per scenario. This supports descriptive replication,
  not precise frequentist coverage claims.
- **Methods:** BED BayesR, BED BayesRC, CSR SBayesR, and CSR SBayesRC only.
- **Simulation:** h2 0.30, expected 50 nonnull markers, mixture variances
  0/0.01/0.1/1, active weights 0.60/0.30/0.10, and the existing seeds.
- **Convergence:** four identifiable chains; rank-normalized R-hat <= 1.01;
  bulk/tail ESS >= 400; relative MCSE <= 0.05; true traces for all core
  variance, alpha, sigma, and derived annotation-prior summaries.
- **Metrics:** the required set enumerated in the audit, with paired
  annotation-minus-baseline and scenario-interaction contrasts.
- **Validation:** exact sample/marker/annotation alignment; deterministic
  simulation oracle; method-control equality; fit completeness; convergence
  gate; output-schema checks; pinned SHAs; capsule checksums.

### Useful optional extensions

- phenotype prediction correlation/NMSE/calibration after the core recovery
  claims are supported;
- component classification scores where the stored posterior probabilities
  make the estimand unambiguous;
- descriptive sigmaSqAlpha behavior, explicitly without a truth-recovery
  claim;
- a separately preregistered binary-only annotation diagnostic if the joint
  design requires scientific decomposition.

### Deferred

- fixed-marker, learned-logistic, and group-specific SBayesC policies;
- overlapping or real biological annotations;
- multitrait annotation models (Study 07 scope);
- broad prediction or model-selection claims based on CPO.

### Remove during migration

- per-study `targets.R` orchestration and CSV-only checkpoint metadata;
- duplicated data, simulation, method, convergence, and provenance helpers;
- the retired Study 02 targets-store dependency in `pilot.R`;
- hard-coded five-replicate assumptions outside the specification;
- inconsistent BED-final-versus-CSR-posterior-mean marker-prior extraction;
- promotion code until a final capsule actually passes every gate.

## Target source contract

### `spec.R`

An ordinary list containing the study ID/task, profiles, pinned data and package
SHAs, sample and marker rules, split and LD settings, two scenarios, five
replicates, exact annotation design and true alpha matrices, four method lists,
all controls and seed rules, convergence contract, estimands, metrics,
validation rules, and the future final capsule path. It must distinguish the
existing `current-stop` development evidence from a future final capsule.

### `annotation-design.R`

Only genuinely annotation-specific science: construction and validation of the
four ordered annotation columns; intercept calibration; conversion of alpha to
marker component probabilities; annotation-aware component sampling; true
annotation/prior tables; and small deterministic fixtures. Data loading,
generic simulation, checkpointing, and fit dispatch do not belong here.

### `analysis.R`

A readable exact workflow that prints provenance, data/scenario/annotation
tables, coordinates and seeds, methods and controls, estimands, convergence
requirements, calls `run_benchmark()`, verifies status/oracles/convergence,
summarizes extracted tables, creates named annotation-recovery and
marker-recovery plots, and prints the output inventory. It should not contain
sampler or extraction internals.

### `report.qmd`

A capsule-only report covering the question, synthetic annotation design,
methods, convergence gate, annotation/prior recovery, marker and effect
recovery, paired contrasts, computational results, validation boundaries, and
reproducibility. Until completion, the existing status page remains explicitly
development evidence.

## Minimal shared-framework extensions

| File | Responsibility | Why shared | Caller | Required tests |
|---|---|---|---|---|
| `benchmark-spec.R` | validate `task = "annotation_models"`, profiles, annotation columns, alpha dimensions, coordinates, and seeds | ordinary task validation belongs with other specs | Study 06 spec/CLI | valid/incomplete specs; profile counts; exact historical coordinates/seeds |
| `benchmark-methods.R` | translate the four plain method lists into existing BED/CSR calls | dispatch conventions and control capture are cross-study | shared runner | controls/priors equal old definitions; unsupported policy fails clearly |
| `benchmark-execution.R` | explicit annotation task branch with validate-only support and predictable tables | one runner should own execution/checkpoint/output mechanics | Study 06 analysis/CLI | validate-only makes zero fit calls; expected 40 coordinates |
| `benchmark-extraction.R` | extract true alpha/sigma traces and consistently label means, final states, and draw-derived marker priors | fit semantics must not differ by interface | Study 06 metrics | BED/CSR fixtures; no mean/final substitution; unavailable quantities explicit |
| `benchmark-reporting.R` | annotation/prior recovery and paired-contrast plots from tidy tables | reusable analysis/report presentation | exact analysis/report | functions return ggplot objects and preserve method order |
| `metrics-annotation.R` | alpha recovery, draw-wise prior recovery/enrichment, PIP recovery, and matched contrasts | a stable annotation-specific metric family exists | Study 06 runner/report | deterministic alpha/prior/PIP fixtures and denominator/formula regression |

`metrics-annotation.R` is justified; it should not become a universal metric
class. Generic prediction, parameter, convergence, and runtime functions should
be reused rather than copied.

## Package-side findings and possible changes

No package change is currently required to call the intended methods or obtain
true alpha/sigma traces. The following narrowly scoped issues must be handled:

1. **Annotation-coefficient mixing (conditional blocker).** The existing BED
   and CSR pilots both failed on the continuous-signal component-0 stick.
   First run the prespecified longer qualification. If it still fails, inspect
   the BayesRC annotation update in the BED core and the SBayesRC annotation
   update in the CSR native backend. The minimal package change must be based on
   that diagnostic--for example a corrected update or validated
   reparameterization--and must pass package fixed-alpha reduction,
   directional-enrichment, seed, threading, and extended-trace tests before the
   Study 06 gate is retried. This blocks completion only if the longer history
   fails.
2. **Final-state naming (nonblocking).** BED's `annotation_prior` is a final
   marker-prior state while CSR lacks a parallel field. A nonbreaking explicit
   `annotation_prior_final` field and documentation would remove ambiguity.
   Study 06 can instead compute the required posterior summaries draw-wise from
   alpha traces, so this is not required for completion.
3. **Unavailable full histories (nonblocking under the recommended metrics).**
   Full marker-by-component posterior histories and comparable global mixture
   traces are not stored. Do not add them unless a required estimand is approved;
   the recommended design does not require them.

## Ordered completion phases

### Phase A — freeze the scientific contract

1. Approve retention of the existing joint binary-plus-continuous informative
   design and the 40-coordinate grid.
2. Freeze annotation-column order, true alpha values, controls, seed tables,
   estimands, metric formulas, and convergence quantities as regression
   fixtures.
3. Record that five replicates support descriptive recovery comparisons rather
   than high-precision coverage claims.

### Phase B — migrate without fitting

1. Create `spec.R` and retain only focused annotation science in
   `annotation-design.R`.
2. Add the minimal shared extensions above and semantic-v2 checkpoints.
3. Replace the old targets graph and duplicated adapters; keep the stop capsule
   unchanged as historical development evidence.
4. Add deterministic tests for alignment, alpha-to-prior conversion, seed
   mapping, extraction semantics, metrics, and validate-only behavior.

### Phase C — qualify annotation MCMC

1. Pre-register one resumable 9,000-iteration maximum history with candidate
   burn-ins 1,000/2,000/3,000 and retained draws 2,000/4,000/6,000 for BED
   BayesRC and CSR SBayesRC at informative replicate 1.
2. Apply the unchanged Study 04 thresholds to all required alpha, sigma,
   variance, and derived prior quantities; choose the earliest passing window.
3. Confirm the selected window on uninformative replicate 1 for both
   annotation methods.
4. If any required quantity fails at the maximum history, stop. Diagnose the
   package implementation; do not run the benchmark or weaken the gate.

The 9,000-iteration ceiling is a recommended planning decision, not an existing
scientific result. It must be approved before execution. Qualification histories
may serve their matching final coordinates only if the selected window and all
semantic identity fields exactly match the frozen final specification.

### Phase D — run the minimum final grid

Run the two scenarios x five replicates x four methods (40 coordinates, four
chains each) using the qualified annotation controls and the existing Study 04
baseline controls. Resume only semantic checkpoints with exact scientific
identity. Stop on failed simulations, alignment, missing traces, fit failures,
or convergence failures.

### Phase E — capsule and reporting

1. Calculate only preregistered metrics and contrasts from extracted tables.
2. Promote one compact final capsule only after all implementation and
   scientific gates pass; exclude native fits and large traces unless required
   for an explicit audit table.
3. Render a capsule-only report with honest validation boundaries and update
   the website status from in development only after checksum and report audits.

## Completion gates

### Implementation

- annotation IDs/order and alpha dimensions validated;
- both annotation methods supported at the pinned package SHA;
- posterior means, final states, and true draws labelled distinctly;
- draw-wise marker-prior extraction implemented consistently for BED and CSR;
- required metrics and semantic-v2 checkpoints tested;
- validate-only execution resolves all 40 coordinates with zero sampler calls.

### Scientific validation

- deterministic simulation oracle and truth-table checks pass;
- annotation qualification passes unchanged thresholds;
- all required fits complete with four identifiable chains;
- null annotation remains null within the preregistered uncertainty criterion;
- uninformative scenario shows no unsupported annotation advantage;
- relevant annotation recovery and annotation-minus-baseline contrasts are
  evaluated without selectively dropping failed quantities.

### Reproducibility

- repository, sblr, qgdata, sample, marker, annotation, and LD identities pinned;
- every simulation, fit, and chain seed fixed;
- capsule manifest, inventory, and checksums pass;
- report reads frozen capsule tables only.

### Reporting

- exact analysis and final report complete;
- all claims correspond to preregistered estimands and passing coordinates;
- synthetic-design and five-replicate limitations stated;
- website status changes only after every prior gate passes.

## Decisions required before implementation

1. Approve retaining the current joint binary-plus-continuous design rather
   than retrospectively narrowing it after the convergence failure.
2. Approve the proposed 9,000-iteration qualification ceiling and candidate
   windows, or preregister a different finite design before fitting.
3. Confirm that the recommended core metric set is sufficient and that CPO,
   group-prior, and multitrait questions remain deferred.
