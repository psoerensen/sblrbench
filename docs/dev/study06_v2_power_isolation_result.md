# Study 06 v2 paired power/isolation diagnostic

**Diagnostic outcome: completed; alpha hierarchy and feedback dominate. This
does not change the failed formal qualification and does not authorize the
final benchmark.**

Status remains:

```text
v1 sparse qualification: failed and preserved
v2 identifiable qualification: failed
paired power/isolation diagnostic: completed
final benchmark: not authorized
```

## Identity and execution contract

The diagnostic started from clean `master` at sblrbench
`954cb8f3d658e2861d07d68523a5d011aa4b9499` and clean sibling sblr
`f2e3647920ed7e8b1ea9d47a6571b3753285682a`. The isolated library was
`results/local/current_benchmark_refresh/rlib`; it loaded sblr 0.2.0 with
`RemoteSha` equal to that sibling SHA. The Study 06 specification hash was
`241c15afab8fefc571e38e625130de6e4ab58b958c30cff369e26998ee30fa56`.

The command was:

```text
Rscript scripts/run_study06_power_isolation.R
```

All eight fits used fit seed 701020, chain seeds
701121/701222/701323/701424, four chains, 9,000 maximum iterations, the
registered window selection and thinning, and the qualification priors and
thread controls. Block-eigen fits used `representation = "low_rank"`, the
cumulative positive-eigenvalue-mass policy, and `eigen_prop = 0.995`. No flat
intercept, sparse CSR, pair allocation, collapsed allocation, changed
initialization, or extended history was used.

No formal-qualification history was reused. The eight successful diagnostic
histories were newly produced in the diagnostic namespace. Subsequent
reporting-only passes reused all eight identity-matched diagnostic checkpoints.

## Shared truth and shuffled control

Every fit used the informative replicate-1 truth. Its aggregate hash was
`169d52bef390022a9106d7e61b200493869b40cbb76f1c8d911eebbc80fea1eb`.

| Item | Value/hash |
|---|---|
| individuals; train/test | 2,000; 1,400/600 |
| markers; blocks | 1,500; 15 blocks of 100 |
| component counts | 1,329 / 84 / 50 / 37 |
| realized heritability | 0.500000 |
| marker order | `135c3604e0c4395349475b8126e0957db265c4b348a2e339daa3e7ddf2316a29` |
| phenotype | `b3de7e139bc58bdd6ddb216e2596477b1c97e2bb95a33f120ed1728cc4f8acbb` |
| effects | `17c9f5e6ef73a626664830d4e7a17e7c255f3524b570a9fb59a194be35511ea9` |
| summary statistics | `2e5a915f7a70f6752b1a19812238e1034a46dc14782163ada1ee78f551cdc5f7` |
| blocks | `34a90410afa84165ac49961ed53b442d7c10c49929f7cb8384e0d9ad6c35a3d0` |

The single shuffled control used prespecified seed 6201. Its row-permutation
hash was `e26553b5e7d616d2075236fb201862f11358ce0eb147263aaaed4c68f7809ff5`
and matrix hash was
`37a74596185a708fe825e478dd251c090587e1e841160867fec828ae965d6e40`.
The intercept remained exact, all non-intercept marginals were identical, and
the maximum annotation-correlation difference was zero. Realized causal rates
were 0.128889 in shuffled-enriched and 0.111373 in shuffled-unenriched markers
(ratio 1.15728); the permutation was not selected using outcomes or truth.

## Registered fits and runtime

| Condition | BED route | Block-eigen route | Seconds BED/block |
|---|---|---|---:|
| no annotation | BayesR | SBayesR | 102.55 / 31.52 |
| learned informative | BayesRC | SBayesRC | 128.13 / 63.25 |
| learned shuffled | BayesRC | SBayesRC | 129.89 / 63.17 |
| fixed true alpha | BayesRC | SBayesRC | 100.55 / 31.91 |

All fits completed without retries. The fixed models used the public
`alpha_init = true_alpha`, `updateAlpha = FALSE` contract. Initialization and
retained traces equalled truth exactly. BED final alpha was exact; the
block-eigen final representation had maximum numerical round-trip error
`2.42e-13`. Marker probabilities recomputed with public
`sblr::sbayesrc_marker_pi()` differed from truth by 0 (BED) and `5.07e-14`
(block eigen). The block result does not expose BED's `annotation_prior`
convenience field, so that route-neutral public-contract audit is explicit.

## Convergence

Thresholds were unchanged: R-hat <= 1.01, bulk and tail ESS >= 400, and
relative MCSE <= 0.05.

| Fit | Burn + retained | Monitored | Failed | max R-hat | min bulk | min tail | max rel. MCSE | Pass |
|---|---:|---:|---:|---:|---:|---:|---:|---|
| baseline BED | 1,000 + 2,000 | 9 | 0 | 1.0047 | 665.4 | 600.5 | 0.0389 | yes |
| baseline block | 1,000 + 2,000 | 9 | 0 | 1.0026 | 1,029.2 | 1,795.7 | 0.0315 | yes |
| informative BED | 3,000 + 6,000 | 26 | 22 | 1.0844 | 34.6 | 33.0 | 0.1770 | no |
| informative block | 3,000 + 6,000 | 26 | 25 | 1.1413 | 23.9 | 39.0 | 0.2058 | no |
| shuffled BED | 3,000 + 6,000 | 26 | 18 | 1.1405 | 27.0 | 51.0 | 0.2124 | no |
| shuffled block | 3,000 + 6,000 | 26 | 19 | 1.0803 | 42.4 | 66.1 | 0.1525 | no |
| fixed-alpha BED | 1,000 + 2,000 | 4 | 0 | 1.0028 | 883.7 | 2,008.4 | 0.0332 | yes |
| fixed-alpha block | 1,000 + 2,000 | 4 | 0 | 1.0028 | 880.6 | 1,229.5 | 0.0357 | yes |

The informative BED failures were effect variance, all 12 alpha coefficients,
stick-1 and stick-2 annotation variances, and all seven marker-prior summaries.
Informative block failed every quantity except genetic variance: effect and
residual variance, heritability, all alpha coefficients, all three annotation
variances, and all prior summaries. Shuffled BED failed nine alpha
coefficients, stick-1/stick-2 annotation variance, and all prior summaries.
Shuffled block failed 11 alpha coefficients, stick-2 annotation variance, and
all prior summaries. The exact quantity-level table is
`results/local/06_annotation_models/v2_paired_power_isolation/selected_convergence.csv`.

This is chain-location disagreement, not merely broad chain-consistent
uncertainty: expected-active summaries differed materially by chain in learned
fits (informative BED 53.62--60.36, informative block 124.86--134.69; shuffled
BED 52.62--58.66, shuffled block 131.87--158.71), while baseline and fixed
chains were tight (baseline BED 167.57--167.94; baseline block 178.57--179.05;
fixed BED 170.93--171.15; fixed block 179.37--179.70).

The fixed-alpha table's four formal quantities are the available retained core
variance traces. Posterior allocation/PIP traces could not be added because of
the public trace abort below. Full-marker per-chain posterior mean PIPs and
their implied expected active counts were nevertheless nearly identical across
the four fixed chains. Accordingly, “fixed-alpha converged” here means the
available selected-window contract passes with strong chain-level allocation
agreement; it is not a claim that unavailable per-iteration occupancy traces
were tested.

### Occupancy limitation

The pinned public selected-component trace control aborted natively for all
1,500 markers, for a deterministic 496-marker set (all 171 causal plus 325
noncausal; estimated extended trace memory 0.998 GiB), and even in an untracked
20-iteration/four-chain/one-marker smoke. No such failed smoke produced
evidence. Complete per-iteration occupancy trajectories are therefore
unavailable. Final states are retained only as explicitly non-convergence
diagnostics; for example, learned informative BED final active counts were
45/46/100/54 and learned shuffled block counts were 82/146/303/115. They were
not substituted for traces.

## Variance and prediction

| Fit | effect variance | genetic variance | residual variance | h2 | validation-g truth cor. | phenotype cor. |
|---|---:|---:|---:|---:|---:|---:|
| baseline BED | 0.01710 | 0.80090 | 1.18136 | 0.40361 | 0.89321 | 0.63140 |
| baseline block | 0.02388 | 1.00166 | 0.98008 | 0.50499 | 0.86978 | 0.60327 |
| informative BED | 0.03063 | 0.82547 | 1.15550 | 0.41624 | 0.91784 | 0.65449 |
| informative block | 0.02862 | 1.01011 | 0.97151 | 0.50927 | 0.89556 | 0.63165 |
| shuffled BED | 0.03503 | 0.76372 | 1.21606 | 0.38532 | 0.87750 | 0.61854 |
| shuffled block | 0.02909 | 0.95059 | 1.03118 | 0.47918 | 0.86206 | 0.59855 |
| fixed BED | 0.01822 | 0.84843 | 1.13524 | 0.42727 | 0.92665 | 0.66164 |
| fixed block | 0.02443 | 1.01826 | 0.96311 | 0.51345 | 0.89943 | 0.63496 |

## Annotation learning and prior contrasts

For informative data, the first-stick learned means (95% intervals) were:

| Route | enriched | continuous | null | sigmaSqAlpha |
|---|---|---|---|---|
| BED | 1.897 [1.182, 2.851] | 0.638 [0.319, 1.045] | -0.042 [-0.514, 0.445] | 2.137 [0.408, 8.242] |
| block | 1.642 [0.826, 2.923] | 0.622 [0.261, 1.112] | -0.371 [-0.907, 0.239] | 1.885 [0.346, 7.406] |

Later-stick coefficients and annotation variances were broad and
chain-inconsistent. Shuffled non-intercept intervals included zero for every
stick. Full coefficient and variance summaries remain in local
`learned_alpha_summary.csv`.

| Fit | expected active | enriched contrast | continuous contrast | causal prior | noncausal prior |
|---|---:|---:|---:|---:|---:|
| informative BED | 57.19 | 0.1912 | 0.01261 | 0.1596 | 0.02250 |
| informative block | 132.42 | 0.2780 | 0.06155 | 0.2642 | 0.06565 |
| shuffled BED | 53.21 | -0.0120 | -0.01705 | 0.03657 | 0.03533 |
| shuffled block | 144.63 | -0.0128 | 0.00342 | 0.09622 | 0.09645 |

These learned-prior summaries are descriptive because their convergence gates
failed.

## Causal-marker identification

PIP AUPRC was primary. Bayesian FDR sorts PIPs descending and chooses the
largest prefix whose cumulative mean `1 - PIP` is at or below the requested
level.

| Fit | AUPRC | AUROC | R@10/25/50/100 | P@10/25/50/100 | median/mean causal rank |
|---|---:|---:|---|---|---|
| baseline BED | .3975 | .6943 | .058/.135/.205/.287 | 1/.92/.70/.49 | 323/492.3 |
| baseline block | .3270 | .6679 | .058/.105/.170/.234 | 1/.72/.58/.40 | 433/527.3 |
| informative BED | .5948 | .8528 | .058/.140/.246/.409 | 1/.96/.84/.70 | 136/281.6 |
| informative block | .5407 | .8300 | .058/.140/.228/.380 | 1/.96/.78/.65 | 168/311.9 |
| shuffled BED | .3013 | .6404 | .058/.117/.158/.205 | 1/.80/.54/.35 | 479/563.9 |
| shuffled block | .2742 | .5835 | .058/.111/.158/.199 | 1/.76/.54/.34 | 660/639.5 |
| fixed BED | .6075 | .8628 | .058/.140/.240/.415 | 1/.96/.82/.71 | 134/268.3 |
| fixed block | .5993 | .8551 | .058/.135/.246/.404 | 1/.92/.84/.69 | 137/278.6 |

| Fit | causal/noncausal mean PIP | ratio | FDR5 selected/true | FDR10 selected/true | effect-truth cor. |
|---|---|---:|---|---|---:|
| baseline BED | .2024/.1002 | 2.02 | 14/14 | 16/16 | .8922 |
| baseline block | .2154/.1069 | 2.02 | 16/16 | 19/18 | .8602 |
| informative BED | .2178/.0161 | 13.52 | 19/19 | 23/22 | .9138 |
| informative block | .3249/.0573 | 5.67 | 23/23 | 32/27 | .8847 |
| shuffled BED | .1196/.0268 | 4.47 | 13/13 | 15/15 | .8781 |
| shuffled block | .1831/.0828 | 2.21 | 16/16 | 19/16 | .8582 |
| fixed BED | .3981/.0775 | 5.14 | 22/22 | 29/27 | .9242 |
| fixed block | .4103/.0823 | 4.98 | 25/23 | 33/28 | .8909 |

By true component, AUPRC for component 1/2/3 was .077/.093/.598 (baseline
BED), .070/.069/.517 (baseline block), .162/.169/.594 (informative BED),
.140/.152/.558 (informative block), .066/.060/.508 (shuffled BED),
.055/.055/.483 (shuffled block), .167/.175/.601 (fixed BED), and
.168/.172/.564 (fixed block). The full component table also records mean PIP
and recall@50/@100.

## Paired power contrasts

Every required contrast includes at least one non-converged learned fit and is
therefore **descriptive under non-converged chains**, not credible.

| Comparison, route | delta AUPRC | AUROC | R@10/25/50/100 | FDR5/FDR10 selected | prediction |
|---|---:|---:|---|---|---:|
| informative - baseline, BED | .1973 | .1585 | 0/.006/.041/.123 | +5/+7 | .0231 |
| informative - baseline, block | .2138 | .1621 | 0/.035/.058/.146 | +7/+13 | .0284 |
| informative - shuffled, BED | .2935 | .2124 | 0/.023/.088/.205 | +6/+8 | .0359 |
| informative - shuffled, block | .2665 | .2465 | 0/.029/.070/.181 | +7/+13 | .0331 |
| fixed - learned, BED | .0127 | .0100 | 0/0/-.006/.006 | +3/+6 | .0071 |
| fixed - learned, block | .0586 | .0251 | 0/-.006/.018/.023 | +2/+1 | .0033 |

Thus annotation-informed fits show a large paired prioritization advantage over
both no-annotation and shuffled controls, and the fixed-alpha upper bound is
at least as strong. That benefit cannot be promoted to a stable method claim
while learned chains fail.

## BED versus block-eigen routes

| Condition | abs h2 delta | prediction-vector cor. | effect cor. | PIP cor. | abs AUPRC delta | active delta |
|---|---:|---:|---:|---:|---:|---:|
| baseline | .10138 | .9692 | .9668 | .8987 | .07050 | 11.06 |
| informative | .09303 | .9611 | .9569 | .8802 | .05408 | 73.06 |
| shuffled | .09386 | .9665 | .9634 | .8389 | .02712 | 85.27 |
| fixed alpha | .08618 | .9722 | .9701 | .9782 | .00816 | 8.52 |

The approximately 0.09--0.10 heritability offset persists without annotations
and with fixed alpha, so it is not caused solely by annotation learning. Both
routes nevertheless converge for baseline and fixed-alpha models. Since every
block retains all 100 modes, this diagnoses blockwise likelihood
factorization, summary-statistics semantics, residual/variance contracts, and
retained-factor execution—not substantial eigenspace truncation.

## Decision and limitations

The primary decision-matrix interpretation is **alpha hierarchy and feedback
dominate**: both BayesR/SBayesR baselines converge, both fixed-true-alpha
BayesRC/SBayesRC fits converge, every learned-alpha fit fails, and fixed alpha
shows clear descriptive annotation power. The persistent route variance offset
is a secondary, separable issue; it did not prevent baseline or fixed-alpha
convergence.

The next recommended task is a package-side, separately authorized audit of
the centred alpha--`sigmaSqAlpha` hierarchy and its update strategy, retaining
this exact truth and comparing proposed changes against the converged fixed
alpha upper bound. A separate blockwise summary-likelihood/residual-variance
calibration audit is also warranted because the no-annotation h2 offset
persists. Neither task should alter this diagnostic or the failed qualification.

Raw checkpoints and compact CSVs remain ignored under
`results/local/06_annotation_models/v2_paired_power_isolation/`. The tracked
machine decision is `docs/dev/study06_v2_power_isolation_decision.json`.
