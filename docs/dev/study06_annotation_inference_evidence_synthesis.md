# Study 06 annotation-inference evidence synthesis

## Scope and status

This is the authoritative chronological synthesis of the committed Study 06
annotation-inference diagnostics in `sblrbench` and the read-only sibling
`sblr` repository through scientific closure. It does not replace either failed
qualification decision or any other historical gate. The immutable v2 specification and informative truth
hashes are respectively
`241c15afab8fefc571e38e625130de6e4ab58b958c30cff369e26998ee30fa56`
and `169d52bef390022a9106d7e61b200493869b40cbb76f1c8d911eebbc80fea1eb`.

Current status:

```text
Study 06 status: CLOSED
Primary final decision: EST-R2
Official qualifier: EST-R5 (formal independent multichain replication unavailable)
Standard same-posterior continuous-alpha sampler development: closed at PMA-R3
Additional final benchmark / scientific fit: not required
Next methodology: separate Bayesian annotation-selection / annotation-PIP study
```

For the authoritative scientific endpoint, read the
[final conclusion](../studies/study06_final_conclusion.md) and
[machine-readable decision](../../results/reference/06_annotation_models/final_decision.json).
For the chronological account, continue here and use the
[main report](../../studies/06_annotation_models/report.qmd). For file-level
provenance and the authoritative/supporting distinction, use the
[documentation inventory](study06_documentation_inventory.md). Formal result
documents and decision JSON files remain authoritative for their individual
experiments; this synthesis does not rewrite those decisions.

## Evidence ledger

“Converged” below always means the recorded four-chain contract (R-hat <=
1.01, bulk/tail ESS >= 400, relative MCSE <= 0.05), except where an experiment
declared its own exact tiny-state criterion.

| Experiment | Question, route, data, and intervention | Result and classification | Ruled out; unresolved | Committed evidence |
|---|---|---|---|---|
| 1. v1 sparse qualification | Can learned BED BayesRC and approximate-CSR SBayesRC qualify on 37,991 markers with about 50 active markers? Informative and uninformative v1 truths; ordinary 9,000-iteration fits. | BED completed but failed 17/23 and 21/23 convergence quantities; CSR stopped with invalid projected residual scale. No scientific comparison was authorized. **Failed development qualification.** | Rules out treating v1 as a completed benchmark. Sparse late-stick identifiability and CSR operator fidelity were confounded. | `docs/dev/study06_annotation_qualification_result.md`; `results/reference/06_annotation_models/current-stop/` |
| 2. v2 identifiable design | Can a prespecified 1,500-marker, 171-active informative truth provide 171/87/37 stick continuations? Fifteen 100-marker blocks, matched annotations, h2 0.50. | All deterministic marker, block, marginal-matching, rank, cell, and truth gates passed. **Design accepted; not an empirical sampler result.** | Rules out accidental empty/rare truth cells as the explanation for v2. Finite-sample posterior identification remained empirical. | `docs/dev/study06_v2_design.md` |
| 3. v2 formal qualification | Do learned BED BayesRC and full-rank retained block-eigen SBayesRC converge and meet scientific/route gates for informative and uninformative truths? | All four fits completed; 74/92 convergence checks failed. Informative later-stick directions and both h2 route gates failed. Informative AUPRC was .5948 BED/.5503 block. **Failed—mixing plus scientific and route blockers.** | Rules out execution failure and numerical non-finiteness as the main v2 outcome. Did not identify whether failure was annotation-specific. | `docs/dev/study06_v2_qualification_result.md` |
| 4. paired power isolation | On the identical informative truth, what changes among no annotation, learned informative, learned shuffled, and fixed-alpha models on BED/block routes? | Baselines and fixed-alpha fits converged; all learned-alpha fits failed. Informative learned AUPRC exceeded matched baselines by .197 BED/.214 block and shuffled controls by .294/.267, but contrasts are **descriptive under non-convergence**. | Rules out a need for annotation learning to obtain the general BayesR baseline and isolates annotation feedback. Does not establish stable learned-model power. | `docs/dev/study06_v2_power_isolation_result.md`; decision JSON |
| 5. fixed-true-alpha comparison | Does heterogeneous marker-specific allocation remain intrinsically non-mixing when alpha is known? Public `alpha_init=true_alpha`, `updateAlpha=FALSE`; same truth and routes. | Both routes passed all available core convergence checks; expected-active counts were about 171 BED/180 block; AUPRC .6075/.5993. **Diagnostic upper bound passed.** | Strongly rules out heterogeneous marker priors alone as sufficient to reproduce the failure. Per-iteration full occupancy traces were unavailable in this earlier run. | paired isolation result and local fixed-alpha audit |
| 6. BayesR/SBayesR baselines | Does general mixture allocation fail without annotation learning? Same truth, BED BayesR and block SBayesR. | Both converged; AUPRC .3975/.3270, h2 .4036/.5050. **Baseline passed.** | Rules out a universal BayesR allocation failure. Leaves a separate likelihood-route calibration difference. | paired isolation result/decision |
| 7. shuffled-annotation controls | Can a non-oracle row permutation create a spurious learned advantage? One joint permutation (seed 6201) preserves annotation marginals/correlations. | Both learned shuffled fits failed convergence and ranked below informative fits (AUPRC .3013/.2742). **Negative-control result descriptive.** | Rules out annotation column marginals alone explaining the informative ranking gain. Non-convergence prevents a calibrated false-positive claim. | paired isolation result/decision |
| 8. mathematical and conditional audit | Are stick eligibility, Albert–Chib coefficient conditionals, variance conditional, intercept handling, or component-probability transforms wrong? Independent R reference versus native scalar owner. | Eligibility counts exact; conditional means/variances agreed to 1e-14; probability transform to 1e-15; 50,000-draw variance quantiles passed. **No conditional defect found.** | Rules out the audited conditional formulas, annotation row alignment, and probability transform. Joint posterior movement remained unresolved. | sibling `docs/dev/study06_alpha_hierarchy_joint_sampling_audit.md` |
| 9. frozen-allocation hierarchy | Does the alpha hierarchy converge when stick outcomes are fixed to truth? F1 current prior, F2 production prior, F3 fixed sigma, F4 fixed alpha; 12,000 iterations. | F1–F3 all passed (max R-hat <=1.0023; min bulk ESS >=2,760); F4 reproduced marker priors exactly. **Frozen hierarchy passed.** | Rules out an isolated alpha/variance conditional failure and a necessary hierarchy funnel when allocations are fixed. Allocation–hierarchy feedback remained. | sibling alpha-hierarchy audit/decision |
| 10. fixed-`sigmaSqAlpha` ablation | Is the centred variance dimension the primary dynamic bottleneck? Dynamic allocations, sigma fixed to 1, both routes. | BED failed 29/32; block failed 29/32. AUPRC .5913/.5448 remained high descriptively. **Decision C support.** | Rules out sampled annotation variance as the primary cause. Alpha/allocation feedback survives. | sibling alpha-hierarchy audit/decision |
| 11. production-equivalent variance-prior ablation | Does official-like `a=b=4` regularization fix dynamic mixing? Both routes, otherwise unchanged. | BED failed 29/35; block failed 32/35. AUPRC .5956/.5578. **Decision C; variance coupling secondary.** | Rules out the teaching prior alone as the cause. The block route still had broader variance/occupancy disagreement. | sibling alpha-hierarchy audit/decision |
| 12. component-trace defect/correction | Were missing occupancy histories caused by sampler failure? Audit selected-marker indexing against the fitted BED subset; add boundary validation only. | A native out-of-range trace index was found and corrected; full 9,000 x 4 x 1,500 traces then retained, with RNG-neutral regression tests. **Implementation defect corrected; sampler transition unchanged.** | Rules out final-state substitution and enables direct occupancy evidence. It does not fix mixing. | sibling alpha-hierarchy audit; component-trace tests/commit `8d0ad4c...` |
| 13a. S1 composition | How does legacy one allocation/one hierarchy update behave with complete traces? Both routes, exact v2 truth. | BED failed 31/35; block 34/35. AUPRC .5948/.5407. **K4 baseline.** | Establishes schedule baseline. | sibling kernel-composition note/decision |
| 13b. H5 composition | Does five hierarchy updates per allocation state improve communication? | BED failed 23/35; block 32/35. Alpha efficiency improved on BED, joint convergence did not. **K4.** | Shows conditional hierarchy equilibration can improve without joint convergence. | same |
| 13c. H20 composition | Does twenty hierarchy updates solve the problem? | BED alpha coefficients passed, but 14/35 joint quantities still failed; block failed 29/35. AUPRC .5945/.5419. **K4.** | Rules out simple under-updating of alpha as a complete solution. | same |
| 13d. A5 composition | Do five allocation sweeps per hierarchy state improve mixing? | BED failed 31/35; block 30/35; runtime-adjusted efficiency worsened. **K4.** | Evidence against merely repeating allocation sweeps. | same |
| 13e. A20 composition | Do twenty allocation sweeps improve mixing? | BED failed 32/35; block 35/35, with max R-hat 1.173/1.234. **K4; unfavorable.** | Stronger evidence that allocation-sweep repetition is not the missing transition. | same |
| 14. tiny exact coupling tempering | Is the three-level coupling/exchange algebra correct on an enumerable model? Eight markers, exact allocation enumeration and alpha quadrature. | PIP, active-count, alpha, detailed-balance, endpoint, and exchange-rate checks passed. **Tiny mechanism validation passed.** | Rules out the tested algebra and implementation on the tiny state space. Scaling to Study 06 remained open. | sibling tempering screen/decision |
| 15. Study 06 BED coupling-tempering screen | Does lambda 0/.5/1 replica exchange bridge Study 06 modes? Four 3,000-cycle ensembles. | 0/2,400 exchanges; median log ratios -348 and -455; target convergence failed, while AUPRC .5898 and prediction remained familiar. **T4: exchange mechanism fails.** | Rules out this sparse three-level ladder. The exact source of the extensive log-density gap was not decomposed. | sibling tempering screen/decision |
| 16. partial-exchange feasibility | Can retained tempering histories determine whether smaller exact exchange blocks would work? Offline only. | Lower-replica complete states were not retained; P1–P8 ratios unavailable, P9 invalid/unavailable. Aggregate count mismatch was not sufficient to explain penalties. **F6: retained state insufficient.** | Rules out retrospective exact partial-exchange claims. Does not rule in or out any proposed partial exchange. | sibling partial-exchange note/decision |
| 17. persistent route calibration | Is the BED/block h2 offset caused by annotation learning or eigen truncation? Baseline, learned, shuffled, fixed-alpha, hierarchy ablations and schedules; all blocks retain 100/100 modes. | Block h2 remains about .09–.10 above BED without annotations and with alpha fixed; predictions/effects remain highly correlated. **Separate route-calibration evidence.** | Rules out annotation learning and substantial eigenspace truncation as sole causes. Summary-likelihood/residual contracts remain unresolved. | paired isolation and sibling hierarchy/kernel notes |
| 18. official single trajectory | Does one default-stream official v0.2.6 path recover the same SNP signal? Exact validated export; D0/D1 matched four-component and D2 native five-component, 9,000/3,000. | D1 versus learned block `sblr`: PIP/effect/validation-g correlations .957/.998/.998; D1 AUPRC .557 versus D0 .320. D0/D2 occupancy and later-stick quantities visibly drift. **GCTB-D1/D3/D5/D6/D7; descriptive only.** | Establishes strong agreement for this one trajectory and reproduces annotation ranking benefit. Does not repair GCTB-P5, establish convergence, or validate latent architecture. | `docs/dev/study06_gctb_single_trajectory_result.md`; decision JSON |
| 19. GCTB-compatible block contract | Does new `gctb_block`/`allMixVe` reproduce pinned official D0/D1 residual and SNP behavior on the exact 1,500-marker export? New S0/S1 single trajectories, all modes. | S0 PIP/effect/g correlations .996/.998/.998; S1 .988/.997/.997. Mean block-Ve differences below .0001 and h2 differences below .026. **SMALL-G1.** | Validates the new contract on the small comparison; does not relabel historical fits or establish official multichain convergence. | `docs/dev/study06_gctb_block_contract_validation_result.md`; decision JSON |
| 20. large information-scale feasibility | Can ordinary BayesRC/SBayesRC converge and recover truth at n=5,000 and m=37,991? One frozen truth, six BED/block fits, all-positive modes. | All six fits completed. Baseline/fixed-alpha aggregate contracts pass; learned E1/B1 fail 23/31 core quantities. AUPRC is .312/.254 learned vs .151/.134 baseline. **LARGE-G2; G3/G4 secondary.** | Establishes that larger information does not remove learned alpha/allocation mixing in this replicate. SNP utility remains useful; route Vg/h2 and occupancy differences persist. | `docs/dev/study06_large_feasibility_result.md`; decision JSON; local compact manifest |

## Key quantitative comparison

Values marked `*` are descriptive under non-converged chains. “Active” is the
posterior mean expected-active count where available. Correlations are effect
with truth, validation genetic value with truth, and phenotype prediction.

| Fit | Conv. | AUPRC | AUROC | Active | h2 | Effect cor. | Val. g cor. | Pred. cor. |
|---|---|---:|---:|---:|---:|---:|---:|---:|
| BayesR BED baseline | pass | .3975 | .6943 | 167.8 | .4036 | .8922 | .8932 | .6314 |
| SBayesR block baseline | pass | .3270 | .6679 | 178.9 | .5050 | .8602 | .8698 | .6033 |
| learned informative BED | fail* | .5948 | .8528 | 58.7 | .4162 | .9138 | .9178 | .6545 |
| learned informative block | fail* | .5407 | .8300 | 131.7 | .5093 | .8847 | .8956 | .6316 |
| fixed-alpha BED | pass | .6075 | .8628 | 171.0 | .4273 | .9242 | .9267 | .6616 |
| fixed-alpha block | pass | .5993 | .8551 | 179.5 | .5134 | .8909 | .8994 | .6350 |
| shuffled BED | fail* | .3013 | .6404 | 56.1 | .3853 | .8781 | .8775 | .6185 |
| shuffled block | fail* | .2742 | .5835 | 141.3 | .4792 | .8582 | .8621 | .5985 |
| BED H20 | fail* | .5945 | .8520 | about 55.6 | .4164 | .9138 | .9178 | .6546 |
| block H20 | fail* | .5419 | .8309 | about 133.6 | .5090 | .8850 | .8958 | .6320 |
| short tempering target | fail* | .5898 | .8484 | about 53.6 | .416 | .913 | .917 | .654 |

Component occupancy is not interchangeable with expected-active count. The
kernel-composition evidence retains both; the paired diagnostic’s earlier
full component-trace request was unavailable until the trace correction.

### Large information-scale comparison

| Fit | Core convergence | AUPRC | AUROC | realized active | Vg | route h2 | effect cor. | in-sample g cor. |
|---|---|---:|---:|---:|---:|---:|---:|---:|
| E0 BED baseline | pass | .151 | .618 | 1,072 | .713 | .349 | .783 | .835 |
| B0 block baseline | pass | .134 | .601 | 1,143 | 1.054 | .515 | .769 | .838 |
| E2 fixed-alpha BED | pass | .304 | .872 | 1,088 | .761 | .372 | .829 | .866 |
| B2 fixed-alpha block | pass | .305 | .872 | 1,142 | 1.049 | .513 | .814 | .859 |
| E1 learned BED | fail* | .312 | .863 | 144 | .719 | .352 | .809 | .834 |
| B1 learned block | fail* | .254 | .840 | 1,093 | 1.048 | .512 | .807 | .855 |

`*` denotes descriptive under nonconvergence. Genetic-value recovery is on the
same 5,000 people used for fitting, not independent prediction. Block `ves` is
mean block residual variance and is not numerically interchangeable with BED
global residual variance.

## Established findings

- On the exact paired informative truth, annotations improve causal-marker
  ranking: learned and fixed-alpha AUPRCs materially exceed both no-annotation
  and shuffled controls on both routes. Learned comparisons remain descriptive
  because those chains did not converge; fixed-alpha comparisons provide the
  converged diagnostic upper bound.
- BayesR/SBayesR baselines converge, and true-alpha-fixed BayesRC/SBayesRC
  converge under the quantities the paired diagnostic could monitor.
- The audited probit eligibility, latent-variable, coefficient, variance, and
  component-probability conditional formulas are correct.
- The hierarchy converges rapidly when allocations are frozen.
- Fixing `sigmaSqAlpha` does not resolve dynamic allocation/alpha mixing.
- Replacing the teaching variance prior with the production-equivalent prior
  does not resolve dynamic mixing.
- Repeating hierarchy updates improves BED alpha convergence: all 12 BED H20
  coefficients passed, although occupancy and effect-variance quantities did
  not. Repeating allocation sweeps is unfavorable on both routes.
- Full joint convergence fails under the ordinary sampler in both the small
  Study 06 experiment and the registered 5,000-person/37,991-marker replicate;
  increased information alone did not resolve it in that replicate.
- The tested three-level complete-state coupling ladder has effectively zero
  overlap on Study 06 despite passing tiny exact validation.
- The retained tempering screen state is insufficient for exact partial-
  exchange or matched-endpoint analysis.
- The roughly .09–.10 BED/block h2 offset persists without annotation learning
  and with alpha fixed; it is separate from the alpha feedback failure.
- One matched official SBayesRC trajectory recovers nearly the same SNP effects,
  PIP ranking, and validation genetic values as learned block `sblr`, and it
  reproduces the annotation prioritization gain over official SBayesR. This is
  an established descriptive trajectory result, not a convergence result.
- The large baseline and fixed-alpha controls pass their aggregate convergence
  contracts. The learned failures therefore remain specifically associated
  with dynamic alpha/allocation feedback rather than general BayesR allocation
  or fixed heterogeneous priors.

## Strongly suggested, not established

- The later sticks have fewer eligible observations and may be only weakly
  informative for quantitative annotation architecture, although first-stick
  failures show this is not solely a late-stick problem.
- Stable PIPs and prediction may require less posterior exploration than alpha,
  variance, and occupancy: ranking/prediction are remarkably stable across
  non-converged schedules and the tempering target. This does not validate the
  latent posterior.
- Official SBayesRC may share the separation between stable SNP outputs and
  unstable latent architecture: the single path supports this interpretation,
  but genuinely independent official chains remain unavailable.
- GCTB’s practical success criterion may emphasize SNP weights/PIPs rather than
  complete latent-state convergence. That is an interpretation of exposed
  outputs, not evidence that unmonitored quantities converge.

## Unresolved questions

- Whether official alpha, occupancy, or SNP outputs converge across genuinely
  independent chains; v0.2.6 cannot answer this through its public seed API.
- Whether the native five-component architecture and its 0.001 component
  improve posterior geometry.
- Whether official residual-variance sampling/tuning changes the result.
- Whether official h2 agrees more closely with the BED or block-eigen `sblr`
  route.
- Whether a partially collapsed or other coordinated posterior-preserving
  transition is needed.
- Whether annotation coefficients can be interpreted quantitatively on this
  dataset even when SNP prioritization is useful.

## Evidence map

| Question | Best experiment | Answer | Remaining uncertainty |
|---|---|---|---|
| Is v2 accidentally unidentified by construction? | v2 truth/block audit | No; all prespecified truth gates pass. | Posterior information can still be weak. |
| Is the issue general BayesR allocation? | paired BayesR/SBayesR baselines | No; both baselines converge. | Route h2 calibration differs. |
| Are heterogeneous marker priors alone enough to fail? | fixed-true-alpha fits | No; both routes converge. | Learned feedback is still complex. |
| Are annotation conditionals wrong? | independent conditional audit | No audited formula defect. | Global scan geometry remains. |
| Is sampled `sigmaSqAlpha` the main funnel? | fixed-sigma dynamic fits | No. | Variance coupling can be secondary. |
| Is the variance prior too weak? | production-prior ablation | Not sufficient; failure persists. | Other prior parameterizations were not explored. |
| Are too few hierarchy updates the complete explanation? | H5/H20 | No; BED alpha improves, joint convergence does not pass. | A coordinated transition may help. |
| Are too few allocation sweeps the explanation? | A5/A20 | No; mixing worsens. | Different allocation kernels remain possible. |
| Can simple coupling tempering bridge modes? | tiny exact + Study 06 screen | Algebra works tiny; the three-level Study 06 ladder has zero exchanges. | Denser or alternative paths are untested. |
| Can existing screen data assess partial exchange? | partial-exchange audit | No; required lower-replica state is missing. | A new retention design would be required. |
| Is the h2 offset caused by annotation/eigen truncation? | baseline and fixed-alpha route pairs | No; it persists and all 100 modes are retained. | Summary-likelihood/residual calibration cause unknown. |
| Does one official path reproduce SNP signal? | official single-trajectory D0/D1/D2 | Yes descriptively: D1 strongly agrees with `sblr` SNP effects/ranking and improves over D0, while architecture differs. | No convergence or independent-chain inference. |
| Does official SBayesRC resolve posterior convergence? | official parity design/smoke | Not answerable: the public seed argument does not produce independent native chains. | Needs an upstream official seedable build or documented independent-chain mechanism. |

## Implementation and sampler lessons

The table below consolidates Study 06-relevant package lessons. Full source and
mathematical detail remain in the read-only sibling `sblr` repository.

| Question | Finding and implication | Status | Evidence |
|---|---|---|---|
| Shared annotation implementation | BED, CSR, and block-eigen scalar routes use the shared annotation-prior owner. Route consistency was audited directly. | correctness established | [`sblr` alpha-hierarchy audit](https://github.com/psoerensen/sblr/blob/master/docs/dev/study06_alpha_hierarchy_joint_sampling_audit.md) |
| Stick eligibility/orientation | Eligible/continue/stop truth counts are 1500/171/1329, 171/87/84, and 87/37/50; native and independent references agree. | correctness established | same |
| Albert–Chib and alpha | Truncated latent variables, residualized Gaussian updates, coefficient moments, and proper intercept handling agree to numerical precision. | correctness established | same |
| `sigmaSqAlpha` | The intercept is excluded and the three non-intercepts enter the exact inverse-chi-square conditional; current and production priors were validated. Conditional correctness does not imply joint convergence. | correctness established | same |
| Component probabilities | Sequential continuation-to-component probabilities, floor, normalization, and orientation pass numerical checks; the floor does not control ordinary Study 06 states. | correctness established | same |
| BED component traces | Full-Glist indices were incorrectly used for a fitted subset. The correction changed trace resolution/boundary validation only, not sampler transitions or RNG order. It did not explain mixing. | defect corrected | same |
| Frozen allocation | The hierarchy passes with current prior, production prior, fixed variance, and fixed alpha when allocations are truth-fixed. | candidate explanation ruled out | decision C |
| Fixed variance | Dynamic learned fits still fail with `sigmaSqAlpha = 1`. | candidate explanation ruled out | decision C |
| Production prior | The production-equivalent variance prior does not restore dynamic convergence. | candidate explanation ruled out | decision C |
| Kernel schedules | H20 improves BED alpha marginals but joint occupancy/variance still fails; A5/A20 worsen behavior. More conditional work is not the missing global transition. | candidate explanation ruled out | [`sblr` K4 audit](https://github.com/psoerensen/sblr/blob/master/docs/dev/study06_allocation_hierarchy_kernel_composition.md) |
| Coupling tempering | Tiny detailed balance and endpoints pass, but the Study 06 three-level ladder accepts 0/2,400 exchanges with gaps of hundreds of log-density units. Correct algebra is an ineffective bridge here. | candidate explanation ruled out | [`sblr` T4 screen](https://github.com/psoerensen/sblr/blob/master/docs/dev/study06_bed_coupling_tempering_screen.md) |
| Partial exchange | Exact formulas were derived, but lower-replica complete state was not retained. No proposal was validated or rejected; aggregate active-count mismatch is not the sole cause. | unresolved contract difference | [`sblr` F6 audit](https://github.com/psoerensen/sblr/blob/master/docs/dev/study06_partial_exchange_feasibility.md) |
| BED/block variance | The h2 offset persists without annotations, with fixed alpha, and with 100/100 modes. | unresolved contract difference | paired isolation; hierarchy/kernel audits |
| Block execution | No MCMC-time eigendecomposition or BED rereading; immutable operators are shared. Recorded block runtime is lower, while analytical and observed process-memory measures differ. | performance observation | [`sblr` block contract](https://github.com/psoerensen/sblr/blob/master/docs/dev/blr_block_eigen_contract.md) |
| Official/`sblr` contracts | Residual variance, `nbsq`/effect scale, p1 orientation, pi semantics, active/component counts, `sigma_anno`/`sigmaSqAlpha`, and wrapper labels remain partly different or unmapped. | unresolved contract difference | official result and source crosswalk |

These safeguards apply throughout: a corrected trace defect is not a mixing
cause when it did not alter transitions; exact conditional tests are not proof
of full posterior convergence; stable SNP outputs do not validate latent
architecture; and package experiments remain distinct from benchmark
experiments.

## Official SBayesRC: strongest delimited comparison

The external implementation is official SBayesRC R implementation v0.2.6 at
source SHA `b95d3fcbad8ff358290922a58fff893439296138`.

| Matched D1 comparison | PIP Pearson | PIP Spearman | Effect Pearson | Validation-g correlation | Top 50 | Top 100 | AUPRC official / `sblr` |
|---|---:|---:|---:|---:|---:|---:|---:|
| learned block `sblr` | .957 | .964 | .998 | .998 | 46 | 86 | .557 / .541 |
| learned BED `sblr` | .933 | .953 | .956 | .960 | 39 | 77 | .557 / .595 |

D1 resembled every saved learned chain: block PIP correlation .946–.965 and
effect correlation .99831–.99833; BED PIP correlation .925–.936 and effect
correlation .95483–.95664. Official D0 versus D1 AUPRC was .320 versus .557,
AUROC .664 versus .842, validation genetic-value correlation .869 versus .896,
and phenotype prediction .604 versus .634.

> The `sblr` SBayesRC implementation compares well with official SBayesRC for
> SNP effects, PIPs, causal prioritization, validation prediction, and practical
> runtime. Agreement is weaker or unresolved for latent architecture,
> later-stick alpha coefficients, component occupancy, and annotation-variance
> parameters.

This is single-trajectory descriptive evidence. Official v0.2.6 does not seed
its native RNG streams through the public seed formal, so formal independently
seeded multichain replication remains unavailable.

## Alpha and annotation-variance crosswalk

The mapped first-stick non-intercept truth is enriched 1.60, continuous .30,
and null 0. The calibrated truth intercept is available in the simulation
object but is not duplicated in the compact result documents.

| Stick-1 coefficient | Truth | Official D1 | learned BED | learned block | Evidence label |
|---|---:|---:|---:|---:|---|
| intercept | available in truth object | -2.406 | unavailable in compact table | unavailable in compact table | official single trajectory |
| enriched | 1.600 | 1.507 | 1.897 | 1.642 | truth / official single path / pooled nonconverged |
| continuous | .300 | .586 | .638 | .622 | same |
| null | 0 | -.156 | -.042 | -.371 | same |

First-stick enriched and continuous directions agree with truth and every
saved learned `sblr` chain; the official null coefficient is comparatively
small. Quantitative equality is not established. BED H20 brought all 12 alpha
marginals inside the diagnostic contract, but remained a nonconverged joint
state. Official p3 intercept and continuous drifted approximately 1.85 and
-1.06 posterior SD, and `sblr` later-stick quantities also mixed poorly.

Official final `sigma_anno` was 1.976/1.043/4.803 (D1) and
3.480/.384/.413/.605 (D2). Pooled nonconverged `sblr` informative
`sigmaSqAlpha` means were about 2.14/1.90/3.14 BED and 2.04/1.54/1.79 block.
Official retained variance histories are not exposed. Final `sigma_anno` values
cannot establish posterior parity, and `sigma_anno` and `sigmaSqAlpha` must not
be assumed identical without a source-contract audit.

## Runtime context

| Method | Trajectories/chains | Iterations | Wall time |
|---|---:|---:|---:|
| Official D0 SBayesR | 1 | 9,000 | 4.61 s |
| Official D1 matched SBayesRC | 1 | 9,000 | 39.28 s |
| Official D2 native SBayesRC | 1 | 9,000 plus tuning | 82.03 s |
| `sblr` informative block SBayesRC | 4 | 9,000 each | 54.74 s |
| `sblr` informative BED BayesRC | 4 | 9,000 each | 129.24 s |

Official runs used one OpenMP and one BLAS thread; `sblr` used four registered
chains and may execute them in parallel. Builds, output retention, diagnostics,
and startup differ. This is not a controlled performance benchmark. Descriptively,
four `sblr` block chains required about 1.4 times one official D1 trajectory.

## Current roadmap

1. Study 06 is closed at **EST-R2**.
2. Standard same-posterior continuous-alpha sampler development is closed at
   **PMA-R3**.
3. No additional Study 06 scientific fit is recommended.
4. The next methodology is a separate Bayesian annotation-selection /
   annotation-PIP study.

The historical `LARGE-F6`, `LARGE-G2/G3/G4`, qualification, official, and
sampler-audit records above remain unchanged evidence of the path to closure.
Official multichain work would require an upstream seedable release or another
documented independent-chain mechanism, but it is not a prerequisite for the
closed primary `sblr` conclusion.
