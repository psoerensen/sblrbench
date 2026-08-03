# Focused CSR SBayesR / GCTB prior diagnostic

## Scientific question

Studies 02, 03, and 06 show suspicious CSR SBayesR variance and prediction
behavior. This diagnostic was designed to isolate five explanations on one
prespecified Study 03 coordinate: the concentrated Dirichlet prior, `adjE`,
initial active probability, variance/residual reporting identities, or another
CSR SBayesR-specific issue.

The diagnostic is **blocked before fitting**. The required variance-prior
isolation cannot be expressed through the current public `sblr` CSR interface.
No result below is posterior evidence about which scientific explanation is
correct.

## Exact provenance

| Item | Value |
|---|---|
| Starting `sblrbench` HEAD | `8864c909a06fae0d44cd5674f8c575682ea2ced7` |
| Installed `sblr` | 0.2.0 |
| Installed `sblr` source SHA | `02e8c74baa906e83c4a08d42a9cc6339b4e81072` |
| qgdata revision | `6cca5819e711d326cfb2614d7e9d9f34942612cd` |
| R | 4.4.1 (2024-06-14 ucrt), x86_64-w64-mingw32 |
| Runtime | Windows 11 x64, build 26200 |
| C/C++ compiler configuration | `gcc` / `g++` |
| Input identity hash | `4ad97132d01634dadbf17abfce958081a0eb3536839eba9261dedf2c0faa7150` |

The diagnostic used the installed package and did not use `pkgload::load_all()`
for `sblr`.

## Simulation coordinate

| Field | Value |
|---|---:|
| Study | `03_parameter_estimation` |
| Architecture | `sparse_mixture` |
| Replicate | 1 |
| Trait | `trait1` |
| Samples | 5,000 |
| Markers after the frozen QC | 37,991 |
| Causal markers | 50 |
| Target / realized heritability | 0.30 / 0.30 |
| Data-selection seed | 3,301 |
| Architecture seed | 7,001 |
| Simulation seed | 7,002 |
| Method seed | 40,104 |
| Chain seeds | 140,104; 240,104; 340,104; 440,104 |
| Summary-data `yy` | 7,205.88881904012 |
| Phenotype variance, `yy / (n - 1)` | 1.44146605701943 |
| Oracle validation | Passed |

The script reused the current Study 03 QC, canonical marker order,
`qgg::getG(impute = TRUE, scale = TRUE)` genotype scaling, sparse-LD settings,
phenotype simulation, and `sblr::make_summary_stats(scale = TRUE)` construction.
It read existing input and LD caches directly and did not invoke a benchmark
target store.

## Fixed baseline prior quantities

The installed resolver was evaluated with the required baseline
`pi = (0.99, 0.01/3, 0.01/3, 0.01/3)`,
`alpha = pi * 500000`, `h2 = 0.30`, and mixture variances
`(0, 0.01, 0.1, 1)`.

| Quantity | Resolved value |
|---|---:|
| Initial mixture-multiplier weight | 0.0037000000000000002 |
| Prior-mean mixture-multiplier weight | 0.0037000000000000002 |
| `B` | 0.0030764029966260011 |
| `E` | 1.0090262399135990 |
| `ssb_prior` | 0.0015382014983130006 |
| `sse_prior` | 0.50451311995679948 |

These are the four values that had to be passed explicitly and confirmed in fit
metadata for every variant.

## Diagnostic variants and automatic confounding

| Variant | Initial active probability | `alpha` | `adjE` | `updatePi` | Initial mixture weight | Prior-mean mixture weight | Automatically resolved `B` | Automatically resolved `ssb_prior` | `ssb_prior` / baseline |
|---|---:|---|---:|:---:|---:|---:|---:|---:|---:|
| current_baseline | 0.0100000 | `pi * 500000` | 0.9 | yes | 0.0037000 | 0.0037000 | 0.0030764030 | 0.0015382015 | 1.0000000 |
| weak_dirichlet | 0.0100000 | `(1,1,1,1)` | 0.9 | yes | 0.0037000 | 0.2775000 | 0.0030764030 | 0.0000205094 | 0.0133333 |
| weak_dirichlet_adjE1 | 0.0100000 | `(1,1,1,1)` | 1.0 | yes | 0.0037000 | 0.2775000 | 0.0030764030 | 0.0000205094 | 0.0133333 |
| fixed_pi_adjE1 | 0.0100000 | `pi * 500000` | 1.0 | no | 0.0037000 | 0.0037000 | 0.0030764030 | 0.0015382015 | 1.0000000 |
| truth_matched_pi | 0.0013161 | `(1,1,1,1)` | 1.0 | yes | 0.0004870 | 0.2775000 | 0.0233751252 | 0.0000205094 | 0.0133333 |

`E = 1.0090262399135990` and `sse_prior = 0.50451311995679948`
resolve identically for these variants. In contrast, changing `alpha` from the
baseline to `(1,1,1,1)` changes the automatic `ssb_prior` by a factor of 75.
The truth-matched initial `pi` also changes automatic `B` by approximately
7.6-fold. Running those variants without explicit fixed priors would therefore
confound mixture-proportion behavior with effect-variance prior calibration.

## Public API isolation failure

The only exported CSR fitting entry point is `sblr::stblr_csr()`. It accepts
model-specific controls through `...`, but its installed BayesR dispatch has no
formal arguments named `B`, `E`, `ssb_prior`, or `sse_prior`. The BayesR backend
`sblr:::stblr_csr_bayesr` is unexported and was inspected only to verify the
limitation. It also has no supported fixed-prior arguments and reports:

```text
Unsupported argument(s) in ...: B, E, ssb_prior, sse_prior
```

Calling the unexported C++ sampler directly would bypass the supported public
interface and violate the protocol. No such call was made.

## Convergence results

| Variant | Expected chains | Completed chains | Fit attempted | Status |
|---|---:|---:|:---:|---|
| current_baseline | 4 | 0 | no | blocked before fitting |
| weak_dirichlet | 4 | 0 | no | blocked before fitting |
| weak_dirichlet_adjE1 | 4 | 0 | no | blocked before fitting |
| fixed_pi_adjE1 | 4 | 0 | no | blocked before fitting |
| truth_matched_pi | 4 | 0 | no | blocked before fitting |

R-hat, bulk ESS, tail ESS, MCSE, and relative MCSE are unavailable. It would be
incorrect to infer convergence from the prior-resolution audit.

## Mixture-proportion behavior

Initial component probabilities and Dirichlet prior means are recorded in
`tables/pi_summary.csv`. Posterior mean `pi`, final-chain `pi`, marker component
probabilities, sampled memberships, and component counts are unavailable
because no fit was attempted. No one of these objects has been relabelled as
another.

## Variance and residual identities

Direct final-state checks of `b' R b`, `b' X'X b / n`, `var(Z b)`, stored `vgs`,
direct residual SSE, the stored residual expression, `ves`, and residual drift
require fitted final effect states. They are deliberately recorded as
unavailable—not as passing or failing—in `direct_variance_checks.csv` and
`residual_checks.csv`.

The denominator contract is already fixed for the future rerun: `n = 5000`,
sample variances use `n - 1`, and `Z` is the exact qgg-standardized genotype
matrix used by Study 03. Posterior-mean-effect identities must remain distinct
from posterior means of draw-level variance quantities.

## Effect and genetic-value recovery

Effect RMSE/correlation, genetic-value correlation/RMSE, predicted genetic
variance, calibration, and PIP-threshold counts are unavailable for all five
variants. The optional BED anchor was not rerun: without any isolated CSR
variant, it cannot resolve the blocking question, and the frozen Study 03 BED
result remains available in the authoritative capsule.

## Variant comparison and interpretation

No variant comparison is scientifically estimable. In particular, the 75-fold
automatic change in `ssb_prior` is evidence of *confounding in the proposed
experiment*, not evidence that weak `alpha` fixes or worsens SBayesR.

The strongest supported interpretation is:

> The current public CSR SBayesR API cannot perform the requested causal
> prior-sensitivity diagnostic without simultaneously changing the variance
> prior. The cause of the suspicious SBayesR results remains unresolved.

This is an implementation-interface limitation, not Outcome 1–5 of the planned
posterior decision framework.

## Limitations

- No SBayesR variant was fitted, so there is no posterior, convergence,
  prediction, recovery, or direct-identity evidence.
- One prespecified simulation was reconstructed only to resolve the exact prior
  scale and prove the experimental confounding.
- The audit describes `sblr` SHA `02e8c74…`; a later interface may differ.
- The local coordinate checkpoint contains compact simulation and summary-data
  inputs, not a benchmark target store or a scientific fit.

## Recommended next action

Open a separate `sblr` implementation task to expose supported explicit
`B`, `E`, `ssb_prior`, and `sse_prior` arguments for public CSR BayesR fitting,
preserve current defaults when those arguments are omitted, and record the
resolved values in fit metadata. Add a small API contract test proving explicit
values are honored. Then rerun this unchanged five-variant design. Do not change
the defaults based on this pre-fit audit alone.

## Reproduction and local evidence

Run:

```powershell
Rscript studies/03_parameter_estimation/diagnostics/sbayesr-gctb-comparison.R
```

The final coordinate rebuild took approximately 116 seconds; the verified
checkpoint-reuse invocation took approximately 15 seconds. The second run
reported `input_checkpoint_reused = true` and zero sampler invocations.

Generated evidence is intentionally untracked:

- `results/local/sbayesr_gctb_diagnostic/manifest.json`
- `results/local/sbayesr_gctb_diagnostic/session_info.txt`
- `results/local/sbayesr_gctb_diagnostic/tables/design.csv`
- `results/local/sbayesr_gctb_diagnostic/tables/prior_resolution.csv`
- `results/local/sbayesr_gctb_diagnostic/tables/fit_status.csv`
- `results/local/sbayesr_gctb_diagnostic/tables/chain_summary.csv`
- `results/local/sbayesr_gctb_diagnostic/tables/pi_summary.csv`
- `results/local/sbayesr_gctb_diagnostic/tables/component_summary.csv`
- `results/local/sbayesr_gctb_diagnostic/tables/scalar_summary.csv`
- `results/local/sbayesr_gctb_diagnostic/tables/convergence_summary.csv`
- `results/local/sbayesr_gctb_diagnostic/tables/effect_recovery.csv`
- `results/local/sbayesr_gctb_diagnostic/tables/direct_variance_checks.csv`
- `results/local/sbayesr_gctb_diagnostic/tables/residual_checks.csv`
- `results/local/sbayesr_gctb_diagnostic/tables/variant_contrasts.csv`

No figures were generated because no fitted posterior quantities exist.
