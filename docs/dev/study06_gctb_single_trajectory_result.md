# Study 06 official SBayesRC single-trajectory result

## Decision and limits

The descriptive diagnostic completed. Supported decisions are **GCTB-D1,
GCTB-D3, GCTB-D5, GCTB-D6, and GCTB-D7**.

> One matched official SBayesRC trajectory recovered the same dominant
> SNP-level signal and annotation-informed prioritization benefit as the
> existing `sblr` learned fits, while active allocation, later-stick annotation
> quantities, and native five-component architecture remained trajectory- or
> implementation-sensitive. This is descriptive evidence only and does not
> establish convergence.

Both statements remain in force:

```text
Official multichain convergence parity:
blocked by the v0.2.6 seed contract.

Official single-trajectory descriptive parity:
completed and interpreted with explicit limitations.
```

No R-hat was calculated. Single-chain ESS and drift do not establish
convergence or exclude another native mode. The formal Study 06 v2 failure and
GCTB-P5 decision are unchanged; the final benchmark remains unauthorized.

## Provenance and identity

The run started from clean `sblrbench` `master` at
`36872d3ea7bc1be991963b614fc70cf05563ada4` and clean read-only sibling `sblr`
`master` at `a165fb0635afcb8a712e8658175dfbb19896b3c3` (version 0.2.0). Official
SBayesRC v0.2.6 was loaded from the existing isolated library at
`results/local/06_annotation_models/gctb_parity/rlib/SBayesRC`; its clean source
checkout remained `b95d3fcbad8ff358290922a58fff893439296138`.

The installed package-tree hash was
`137b3bca085651aa07908d5ba235f5b0a6f5be27fa48ede97efa1e00c31e23ee`
and its installed DESCRIPTION hash was
`b3c90d91696cb5077e715a773edbe124b20d3f9f7b5a037fe8a1c676e6b0b71b`.
The DESCRIPTION lacks `RemoteSha`; the clean checkout and earlier committed
export/smoke provenance supply the SHA pin.

R was 4.4.1 UCRT on x86_64 Windows. Rtools GCC 13.2 exposes OpenMP 4.5
(`_OPENMP=201511`). R used its Windows internal BLAS/LAPACK; every child set
`OMP_NUM_THREADS`, `OPENBLAS_NUM_THREADS`, `MKL_NUM_THREADS`, and
`BLAS_NUM_THREADS` to one.

All six immutable Study 06 hashes matched the design. Every one of the 20
exported files retained its manifest SHA-256. N was 1,400 for all 1,500
markers; marker order, A1, annotations, 15 x 100 block membership, and 100/100
positive-mode retention were exact. No export was recreated and no allele flip
was required.

## Execution record

All fits used 9,000 iterations, 3,000 burn-in, and 6,000 retained draws.
Nominal requested seeds are provenance only: official v0.2.6 does not connect
them to its native RNG streams.

| ID | Configuration | nominal seed | native wall time | warnings |
|---|---|---:|---:|---|
| D0 | matched four-component SBayesR; all modes | 711121 | 4.61 s from official log | none from official fit; post-fit SHA JSON coercion was repaired and the completed fit registered without rerun |
| D1 | matched four-component SBayesRC; informative annotations; all modes | 721121 | 39.28 s | benign wrapper label recycling for four-component NumSnp/Vg summaries; native dimensions used |
| D2 | native five-component SBayesRC plus tuning | 731121 | 82.03 s | none |

D0's first 40 native states were exactly equal to the previous D0 smoke for
`hsq_hist`, `ssq_hist`, `pi_hist`, `vg_comp_hist`, and `n_comp_hist` (560
values, maximum absolute difference zero). Post-burn-only annotation files and
sparse beta histories cannot overlap a smoke with a different burn-in. No
second full run was launched.

D2 tuning correlations were .55928, .57608, .56185, and .55755 for thresholds
.995, .99, .95, and .9. Relative correlations were 1, 1.03004, 1.00460, and
.99691. The native rule selected .995 because improvement was below 1.25; all
100 modes per block remained retained. The threshold matched the previous
smoke, although the tuning correlations did not repeat because the official
wrapper constructs its pseudo-summary tuning input with unseeded R `rnorm()`.
The ignored public seed therefore fixes neither that R-level tuning input nor
native streams. This does not alter the fixed D0 native-stream check.

## Exposed outputs

The comparison used posterior mean/SD/last SNP effects, SNP PIPs, sparse beta
history, `pi_hist`, component counts, component genetic-variance fractions,
`hsq_hist`, effect variance, block residual/effect/genetic-variance histories,
alpha histories, conditional and joint annotation probabilities, enrichment,
and official final summaries. Retained `sigmaSqAlpha` history is not exposed.
Four-component results use native array dimensions because wrapper labels
recycle NumSnp2:5/Vg2:5.

D1 alpha mapping is p1/p2/p3 x Intercept/enriched/continuous/null. D2 adds p4
with the same annotation order. Official p1 is the conditional probability of
leaving the null component; the marker-prior reconstruction follows that
source-verified orientation.

## Within-trajectory description

| ID | h2 mean (first/final third) | active mean (first/final third) | substantial drift flags | interpretation |
|---|---:|---:|---:|---|
| D0 | .465 (.455/.463) | 486 (229/509) | 6/17 | h2 visually stable; active allocation switched slowly between sparse regimes |
| D1 | .461 (.462/.457) | 78.9 (79.6/74.3) | 3/57 | SNP/variance summaries comparatively stable; p2-null and p3 intercept/continuous visibly drifted |
| D2 | .456 (.460/.452) | 749 (735/482; middle 1,030) | 21/72 | strong occupancy, mixture, p3/p4 intercept, and probability drift |

D0 active-count lag-1 autocorrelation was .9973 and single-chain bulk ESS was
3.9. D1 active-count standardized first-to-final drift was -.35 SD (ESS 42.6),
while p3 intercept drifted 1.85 SD. D2 active-count drift was -.80 SD (ESS 6.5)
and several annotation probabilities exceeded .75 SD drift. These are
**visible-drift descriptions**, never convergence diagnoses.

## SNP-level agreement

| Comparison | PIP Pearson/Spearman | effect Pearson/Spearman | validation-g agreement | top-50/top-100 overlap |
|---|---:|---:|---:|---:|
| C0 D0 vs `sblr` block SBayesR | .997/.841 | .999/.997 | .999 | 47/96 |
| C1 D1 vs learned block SBayesRC | .957/.964 | .998/.973 | .998 | 46/86 |
| C2 D1 vs learned BED BayesRC | .933/.953 | .956/.740 | .960 | 39/77 |
| C3 D1 vs fixed-alpha block | .884/.951 | .996/.968 | .997 | 43/84 |
| C3 D1 vs fixed-alpha BED | .844/.926 | .966/.763 | .970 | 37/76 |
| C4 D2 vs D1 | .542/.493 | .9995/.968 | .9996 | 35/63 |
| C5 D1 vs D0 | .778/.284 | .968/.895 | .968 | 30/47 |

For C1 PIP RMSE was .061; top-10/25/50/100/200 overlaps were
10/24/46/86/178. For C2 PIP RMSE was .053 and overlaps were
8/22/39/77/165. Dense near-zero PIPs account for part of rank-correlation
behavior; scatter and decile-bin plots are retained locally. The effect slope
was not post-hoc rescaled. C1's near-unity effect and genetic-value agreement
supports GCTB-D1, while the non-identical PIP scale supports cautious
calibration language rather than formal equivalence.

Against individual learned `sblr` chains, D1 PIP correlation ranged
.946-.965 for block chains and .925-.936 for BED chains. Effect correlation
ranged .99831-.99833 and .95483-.95664; validation-g agreement ranged
.99825-.99828 and .95928-.96079. D1 resembles every saved chain rather than a
single isolated `sblr` occupancy regime. The learned `sblr` fits themselves
failed multichain convergence, so these are per-chain descriptive agreements.

## Causal ranking, FDR, and prediction

| Method | AUPRC | AUROC | P/R @25 | P/R @50 | P/R @100 | P/R @200 | median causal rank |
|---|---:|---:|---:|---:|---:|---:|---:|
| D0 | .320 | .664 | .72/.105 | .58/.170 | .39/.228 | .275/.322 | 478 |
| D1 | .557 | .842 | .92/.135 | .78/.228 | .64/.374 | .50/.585 | 155 |
| D2 | .429 | .764 | .88/.129 | .70/.205 | .52/.304 | .37/.433 | 263 |
| `sblr` learned block* | .541 | .830 | .96/.140 | .78/.228 | .65/.380 | .48/.561 | see local table |
| `sblr` learned BED* | .595 | .853 | .96/.140 | .84/.246 | .70/.409 | .53/.620 | see local table |
| fixed-alpha block | .599 | .855 | -- | .84/.246 | .69/.404 | .52/.608 | see local table |
| fixed-alpha BED | .607 | .863 | -- | .82/.240 | .71/.415 | .53/.620 | see local table |

`*` Descriptive pooled summaries from non-converged learned chains.

At Bayesian FDR 5%, D0 selected 16/16 true, D1 selected 22 with 21 true, and
D2 selected 30 with 25 true. At 10%, those counts were 20/18, 28/25, and
45/32. D1 shared 21 of 22/23 5%-FDR discoveries with learned block `sblr` and
18 of 22/19 with learned BED; it shared all 16 D0 5%-FDR discoveries.

D1 improved over D0 by .237 AUPRC and .178 AUROC. Its validation genetic-value
correlation with truth increased from .869 to .896 and phenotype prediction
from .604 to .634; prediction RMSE decreased from .526 to .473. That supports
GCTB-D6. D1 mean causal/noncausal PIP was .251/.027 versus D0 .396/.315: D1's
ranking improvement arose largely from sharply suppressing noncausal PIPs.

The official effect correlations with truth were .859 (D0), .885 (D1), and
.886 (D2). Causal-only correlations were .879, .904, and .904. D1 validation
prediction had slope 1.055, RMSE .473, and predicted-genetic variance .807.

## Architecture and annotations

| Method | mean active | h2 | total block genetic variance | mean block residual variance |
|---|---:|---:|---:|---:|
| official D0 | 486 | .465 | .465 | 1.000 |
| official D1 | 78.9 | .461 | .461 | 1.000 |
| official D2 | 749 | .456 | .456 | 1.000 |
| `sblr` learned block* | about 132 | .509 | about 1.00 | about .98 |
| `sblr` learned BED* | about 59 | .416 | about .83 | about 1.18 |

Official D1's active count is closer to the BED regime, its h2 lies between BED
and block, and its SNP effects align most closely with the block route. It is
therefore not wholly on either `sblr` variance scale. SNP agreement coupled to
architecture differences supports GCTB-D3.

D1 p1 alpha means (SDs) were intercept -2.406 (.229), enriched 1.507 (.273),
continuous .586 (.160), and null -.156 (.185). The reconstructed expected
active count was 79.25; enriched and unenriched mean non-null priors were .223
and .0228. Enriched and continuous directions agree with truth and `sblr`; the
null coefficient is small relative to the enriched signal. D1 later sticks
were unstable: p2 intercept 10.56 and p3 intercept 8.95, with p3 intercept and
continuous drifting 1.85 and -1.06 SD. They are not quantitatively
interpretable from this trajectory.

D2 p1 enriched and continuous means were positive (.698 and .414), but its
null mean was also positive (.423), expected active count was 749, and several
p1 probabilities drifted substantially. D2 p2 showed the clearest informative
directions; p3/p4 were weak and drifting. D2 and D1 differ in PIP scale,
ranking, active allocation, and alpha geometry while changing components,
starting mixture, and tuning simultaneously. GCTB-D7 is supported, but no
single changed feature receives causal attribution.

## Evidence and recommendation

The ignored manifest hashes the package tree, all inputs by inherited export
identity, all three configurations, raw outputs, extracted tables/plots, and
both scripts. Raw scientific output remains under
`results/local/06_annotation_models/gctb_single_trajectory/`; it is not a
tracked posterior capsule.

The next focused task should audit official-versus-`sblr` latent architecture
contracts on the matched four-component summary route: trace the residual
variance, `nbsq` effect scale, p1 continuation orientation, pi update, and
active-count definitions on the same fixed early native states. Preserve SNP
outputs and do not change either sampler until a concrete contract difference
is identified. A separately released official seed contract is still required
before multichain parity can be attempted.
