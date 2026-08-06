# Study 06 documentation inventory

## Scope and conventions

This inventory records the tracked Study 06 material reviewed during the
documentation consolidation at `sblrbench` commit
`d9ecfe20368daff2eec2a01a73c7a13b2dfa8f31` and read-only sibling `sblr`
commit `a165fb0635afcb8a712e8658175dfbb19896b3c3`. “README” and “report” links
refer to `studies/06_annotation_models/README.md` and `report.qmd`. Historical
records are retained even when their next-step wording has been superseded.

Evidence categories are the controlled vocabulary requested for this cleanup.
The compact “action” field records whether a file is linked, summarized,
retained as historical evidence, or left untouched.

## Benchmark documentation and navigation

| File | Title/purpose | Category | Status and authority | Linked README/report | Terminology, duplication, stale wording, decision | Cleanup action |
|---|---|---|---|---|---|---|
| `README.md` | Repository overview | navigation/index | Supporting | yes/indirect | Current but duplicated status prose | Short status and links |
| `studies/README.md` | Study-directory contract | navigation/index | Supporting | indirect/indirect | Previously too terse | Add Study 06 navigation link |
| `studies/index.qmd` | Website catalogue | navigation/index | Supporting | indirect/yes | Previously stopped at qualification | Summarize current diagnostic status |
| `_quarto.yml` | Site build/navigation | navigation/index | Authoritative build config | n/a | Current; no Study 06 defect | Leave untouched |
| `studies/06_annotation_models/README.md` | Study index and status | navigation/index | Authoritative navigation | self/yes | Needed hierarchy and reading paths | Rewrite as compact index |
| `studies/06_annotation_models/report.qmd` | Primary readable report | main report | Authoritative overview | yes/self | Mixed v1 detail with newer conclusions | Reorganize and consolidate |
| `docs/dev/study06_annotation_inference_evidence_synthesis.md` | Chronological technical synthesis | authoritative synthesis | Authoritative synthesis | yes/yes | Current ledger; missing consolidated implementation/runtime/roadmap sections | Extend without replacing ledger |
| `docs/dev/study06_documentation_inventory.md` | This file | navigation/index | Authoritative cleanup inventory | yes/yes | New | Link |
| `docs/dev/study06_documentation_cleanup_result.md` | Cleanup record | formal decision | Authoritative for this documentation task | yes/yes | New | Link |
| `docs/dev/study06_documentation_cleanup_decision.json` | Machine cleanup decision | formal decision | Authoritative for this documentation task | yes/yes | New | Link |

## Benchmark experiment records

| File | Title/purpose | Category | Status and authority | Linked README/report | Terminology, duplication, stale wording, decision | Cleanup action |
|---|---|---|---|---|---|---|
| `docs/dev/study06_annotation_qualification_result.md` | v1 qualification result | formal qualification | Immutable authoritative historical decision | yes/yes | Formal failure must remain | Link; untouched |
| `docs/dev/study06_v2_design.md` | Identifiable v2 protocol | design/protocol | Authoritative design | yes/yes | Current | Link; untouched |
| `docs/dev/study06_v2_qualification_result.md` | v2 qualification result | formal qualification | Immutable authoritative decision | yes/yes | Formal failure must remain | Link; untouched |
| `docs/dev/study06_v2_power_isolation_result.md` | Eight-fit paired diagnostic | diagnostic experiment | Authoritative diagnostic result | yes/yes | Later package audit supersedes its “next task,” not its decision | Link; untouched |
| `docs/dev/study06_v2_power_isolation_decision.json` | Paired diagnostic decision | formal decision | Machine-authoritative | yes/yes | Decision unchanged | Link; untouched |
| `docs/dev/study06_gctb_parity_design.md` | Official parity protocol | design/protocol | Authoritative design | yes/yes | Current | Link; untouched |
| `docs/dev/study06_gctb_parity_result.md` | Official export/API diagnostic; GCTB-P5 | diagnostic experiment | Authoritative blocker result | yes/yes | Multichain blocker must remain | Link; untouched |
| `docs/dev/study06_gctb_parity_decision.json` | GCTB-P5 decision | formal decision | Machine-authoritative | yes/yes | Decision unchanged | Link; untouched |
| `docs/dev/study06_gctb_single_trajectory_design.md` | Official descriptive protocol | design/protocol | Authoritative design | yes/yes | Current | Link; untouched |
| `docs/dev/study06_gctb_single_trajectory_result.md` | Official single-path result | descriptive comparison | Authoritative within single-trajectory scope | yes/yes | Must not be called convergence | Link; untouched |
| `docs/dev/study06_gctb_single_trajectory_decision.json` | GCTB-D decisions | formal decision | Machine-authoritative | yes/yes | Decisions unchanged | Link; untouched |
| `docs/dev/study06_annotation_audit.md` | Early implementation audit | implementation/correctness evidence | Supporting historical | yes/yes | Some pre-v2 framing | Retain and label supporting |
| `docs/dev/study06_annotation_implementation.md` | Initial implementation contract | design/protocol | Supporting historical | yes/yes | Detailed and partly superseded | Retain; link from implementation history |
| `docs/dev/study06_annotation_completion_plan.md` | Original completion plan | historical evidence | Historical, not current roadmap | yes/yes | Stale next-step wording by design | Retain historical; do not rewrite |
| `docs/dev/study06_migration.md` | CSR/operator migration history | historical evidence | Supporting historical | yes/yes | Approximate CSR no longer v2 gate | Retain historical |
| `docs/dev/sbayesr_gctb_diagnostic.md` | Earlier focused CSR/GCTB prior diagnostic | diagnostic experiment | Supporting historical | yes/yes | Different route and scope | Retain; link as precursor |
| `docs/dev/study06_annotation_model_support.csv` | Interface support matrix | design/protocol | Supporting machine-readable | yes/yes | Current only in its recorded context | Link; untouched |
| `docs/dev/study06_annotation_required_runs.csv` | Original run registry | design/protocol | Supporting historical | yes/yes | Superseded by versioned registries | Retain historical |
| `docs/dev/study06_annotation_results_inventory.csv` | Original results inventory | historical evidence | Supporting | yes/yes | Superseded by this broader inventory | Retain and cross-link |

## Study code, scripts, and tests

| File | Purpose | Category | Status and authority | Linked README/report | Notes | Cleanup action |
|---|---|---|---|---|---|---|
| `studies/06_annotation_models/spec.R` | Versioned v1/v2 specification | script/test | Scientific source of truth | yes/yes | Current | Link; untouched |
| `studies/06_annotation_models/analysis.R` | Study analysis/runner integration | script/test | Executable support | yes/yes | Current | Link; untouched |
| `studies/06_annotation_models/annotation-design.R` | Annotation/truth construction | script/test | Executable support | yes/yes | Current | Link; untouched |
| `studies/06_annotation_models/power-isolation.R` | Paired diagnostic helpers | script/test | Executable support | yes/yes | Current | Link; untouched |
| `studies/06_annotation_models/gctb-parity.R` | Official export/API helpers | script/test | Executable support | yes/yes | Current | Link; untouched |
| `studies/06_annotation_models/gctb-single-trajectory.R` | Official descriptive-analysis helpers | script/test | Executable support | yes/yes | Current | Link; untouched |
| `scripts/run_study06_power_isolation.R` | Paired diagnostic entry point | script/test | Reproducibility support | yes/yes | Must not be run in cleanup | Link; untouched |
| `scripts/run_study06_gctb_parity.R` | Official parity entry point | script/test | Reproducibility support | yes/yes | Must not be run in cleanup | Link; untouched |
| `scripts/run_study06_gctb_single_trajectory.R` | Official single-path entry point | script/test | Reproducibility support | yes/yes | Must not be run in cleanup | Link; untouched |
| `scripts/run_benchmark.R` | Shared benchmark entry point | script/test | Framework support | yes/yes | Study 06 final mode remains blocked | Link only |
| `scripts/run_benchmark.ps1` | PowerShell wrapper | script/test | Framework support | no/no | Not Study 06-specific | Leave untouched |
| `tests/testthat/test-study06-annotation-models.R` | Study 06 contracts | script/test | Verification | yes/yes | Add documentation/status contract | Extend minimally |

## Immutable v1 reference evidence

All 19 files below are **raw reference evidence**, historical failed-development
evidence, authoritative for v1 only, linked through the capsule README and
main report, current in their historical context, and left untouched.

| File | Purpose |
|---|---|
| `results/reference/06_annotation_models/current-stop/README.md` | Capsule scope and stop decision |
| `annotation_design_summary.csv` | Annotation design summary |
| `benchmark_manifest.json` | Identity and provenance manifest |
| `candidate_settings.csv` | Candidate convergence windows |
| `checksums.csv` | Capsule checksums |
| `computational_summary.csv` | Recorded runtime/resource summary |
| `config.R` | Frozen configuration |
| `contract_smoke_test.R` | Frozen contract smoke |
| `convergence_diagnostics.csv` | Quantity-level diagnostics |
| `example_data_manifest.csv` | Data identity |
| `fit_status.csv` | Fit completion status |
| `interface_audit_sources.csv` | Interface evidence sources |
| `method_recommendations.csv` | Historical stop recommendations |
| `reproduce.R` | Capsule validation entry point |
| `scalar_chain_draws.csv` | Compact retained scalar evidence |
| `seed_registry.csv` | Frozen seeds |
| `session_info.txt` | Software environment |
| `source_files.csv` | Source provenance |
| `true_alpha.csv` | v1 generating alpha |

## Read-only sibling `sblr` evidence

These files remain package-side evidence and are linked rather than copied.
Their low-level implementation detail remains authoritative in `sblr`.

| File in `psoerensen/sblr` | Purpose | Category | Authority/status | Study 06 action |
|---|---|---|---|---|
| `docs/dev/study06_alpha_hierarchy_joint_sampling_audit.md` | Conditional audit, trace correction, frozen/dynamic ablations | implementation/correctness evidence | Authoritative C decision | Summarize and link |
| `docs/dev/study06_alpha_hierarchy_decision.json` | C machine decision | formal decision | Authoritative | Link; unchanged |
| `docs/dev/study06_allocation_hierarchy_kernel_composition.md` | S1/H5/H20/A5/A20 audit | diagnostic experiment | Authoritative K4 | Summarize and link |
| `docs/dev/study06_allocation_hierarchy_kernel_decision.json` | K4 machine decision | formal decision | Authoritative | Link; unchanged |
| `docs/dev/study06_bed_coupling_tempering_screen.md` | Tiny exact validation and Study 06 ladder screen | diagnostic experiment | Authoritative T4 | Summarize and link |
| `docs/dev/study06_bed_coupling_tempering_decision.json` | T4 machine decision | formal decision | Authoritative | Link; unchanged |
| `docs/dev/study06_partial_exchange_feasibility.md` | Offline state-sufficiency audit | diagnostic experiment | Authoritative F6 | Summarize and link |
| `docs/dev/study06_partial_exchange_decision.json` | F6 machine decision | formal decision | Authoritative | Link; unchanged |
| `docs/dev/blr_block_eigen_contract.md` | Retained-factor operator contract | implementation/correctness evidence | Authoritative package contract | Link |
| `docs/methods/block_eigen_operator.qmd` | User-facing operator semantics | design/protocol | Authoritative method doc | Link |
| `docs/methods/annotation_priors.qmd` | Annotation preprocessing/prior semantics | design/protocol | Authoritative method doc | Link |
| `tools/study06_alpha_hierarchy_reference.R` | Independent conditional reference | script/test | Supporting exact validation | Retain in sibling |
| `tools/study06_alpha_hierarchy_audit.R` | Hierarchy audit driver | script/test | Supporting | Retain in sibling |
| `tools/study06_kernel_composition_audit.R` | Schedule audit driver | script/test | Supporting | Retain in sibling |
| `tools/study06_bed_coupling_tempering_screen.R` | Tempering screen driver | script/test | Supporting | Retain in sibling |
| `tools/coupling_tempering_tiny_reference.R` | Enumerable-state tempering reference | script/test | Supporting exact validation | Retain in sibling |
| `tools/study06_partial_exchange_feasibility.R` | Offline feasibility audit | script/test | Supporting | Retain in sibling |
| `tests/testthat/test-alpha-hierarchy-conditionals.R` | Conditional regression tests | script/test | Supporting correctness | Retain in sibling |
| `tests/testthat/test-bayesrc-coupling-tempering.R` | Tempering regression tests | script/test | Supporting correctness | Retain in sibling |
| `src/st_bayesrc_annotation_prior.h` | Shared scalar annotation kernel | implementation/correctness evidence | Authoritative source | Link conceptually; do not copy |

## Duplication and stale-language findings

- The repository README, study catalogue, Study 06 README, report, and synthesis
  repeated slightly different status summaries. They now use one status block.
- The report devoted most of its narrative to v1 and appended newer work. It is
  now organized by evidence level and experiment chronology; the v1 capsule
  remains rendered but is explicitly historical.
- Experiment-specific “next task” statements are historically correct at their
  commit. They remain untouched. The synthesis and README now own the current
  roadmap.
- “GCTB parity blocked” was ambiguous. Navigation now distinguishes blocked
  multichain convergence parity from the completed official single-trajectory
  descriptive comparison.
- “Active count” was used in contexts mixing expected non-null probability and
  sampled component occupancy. Principal documents now name these separately.
