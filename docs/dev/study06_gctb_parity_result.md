# Study 06 official SBayesRC parity result

## Decision

**GCTB-P5: blocked by the official v0.2.6 independent-chain contract.**

The exact Study 06 export and each G0/G1/G2 execution path passed a short
non-inferential smoke. However, two fresh G0 processes given distinct
preregistered requested seeds, 711121 and 711222, produced identical scientific
result signatures:
`087eca33bd0fa0faf6f3387ad9462f1d37422750010d9dece912554dbbc8c5b9`.
Pinned source inspection confirms that the wrapper declares `seed` but neither
uses it nor forwards it; native static/thread-local random engines are default
constructed and expose no public seed connection.

Four independent chains therefore cannot be produced without modifying the
official implementation, which is outside scope. The 12-coordinate registry
was not run. No convergence, causal-ranking, annotation-direction, h2, or
cross-implementation parity conclusion is drawn from the smoke states.

```text
v1 sparse qualification: failed and preserved
v2 identifiable qualification: failed
paired power isolation: completed
package-side hierarchy and transition audits: completed
official SBayesRC parity diagnostic: blocked
final benchmark: not authorized
```

## Provenance and installation

The diagnostic started from clean `sblrbench` `master` at
`de31f62e182d8540488d4135df4c58f052a515d9` and clean read-only sibling `sblr`
`master` at `a165fb0635afcb8a712e8658175dfbb19896b3c3`. The sibling version is
0.2.0. The official checkout is release v0.2.6 at
`b95d3fcbad8ff358290922a58fff893439296138`.

Official SBayesRC was installed only at
`results/local/06_annotation_models/gctb_parity/rlib/SBayesRC`. A first build
against BH 1.90 failed because that Boost requires C++14 while the package asks
R for C++11. The successful unchanged-source build used BH 1.87.0-1 and an
ignored Makevars mapping the requested standard to GNU++14. R was 4.4.1 UCRT
on x86_64 Windows; compilation used Rtools GCC 13.2 and `-fopenmp`. The
installed DESCRIPTION exposes version 0.2.6 and the official repository URL,
but no `RemoteSha`, so the clean source HEAD and manifest supply the SHA pin.
The successful installation was equivalent to:

```powershell
$env:R_MAKEVARS_USER = "results/local/06_annotation_models/gctb_parity/build-Makevars"
R CMD INSTALL --library=results/local/06_annotation_models/gctb_parity/rlib results/local/06_annotation_models/gctb_parity/official-source
```

## Export validation

The export reproduced specification
`241c15afab8fefc571e38e625130de6e4ab58b958c30cff369e26998ee30fa56`
and truth
`169d52bef390022a9106d7e61b200493869b40cbb76f1c8d911eebbc80fea1eb`.

| Identity | Hash/value |
|---|---|
| marker order | `135c3604e0c4395349475b8126e0957db265c4b348a2e339daa3e7ddf2316a29` |
| GWAS export | `1bd0abb220c4e7f9ca58ed11f2c2913e6e142868ba17820b988670cd01cd4610` |
| annotation matrix | `b2442b5c074b1cd6fa0eb047fcedae6d3296a15d2b5f985c2ac288b534a3b156` |
| LD audit/block identity | `369df6bbed513da7913f960b9f79e35d96cdeac321cf3105b303c932d8e2c0a4` |
| sample size | 1,400 for every marker |
| blocks/modes | 15 x 100; all 100 positive modes in every block |
| source LD reconstruction | maximum double-precision error `2.65e-14` |
| official-reader reconstruction | maximum float-format error `4.01e-08` |
| symmetry error | 0 in every block |
| diagonal range across blocks | 0.9204 to 1.1176 |

GWAS marker order, alleles, frequencies, effect direction, and cross-products
matched the Study 06 training objects. Official input alignment required zero
allele flips. The annotation reader received the required explicit intercept
slot and all three non-intercept columns; output headers verified the mapping.
G0/G1 used all positive modes. The G2 smoke selected threshold 0.995, which
still retained all 100 modes per block on this panel.

## Registry and smoke record

The unrun registry contains G0/G1/G2 x four intended chains, with requested
seeds formed from the registered Study 06 seeds plus fixed 10,000/20,000/30,000
condition offsets. It specifies 3,000 iterations and 1,000 burn-in exactly as
described in the design.

The first smoke-orchestrator attempt made each child regenerate the Study 06
truth. The child hash guard rejected those processes before the official fit
call. The runner was corrected so only the parent reconstructs, validates, and
exports truth; fresh children consume the hash-checked immutable export and
registry. This is an execution-scaffolding fix, not a simulation or model
change. No rejected child output was reused.
One direct G0-chain-1 smoke was then used to verify the child path; the final
orchestrator repeated that same 40-iteration non-inferential identity and
replaced only its ignored smoke directory. No 3,000-iteration coordinate was
run or retried.

| Smoke | Requested seed | Wall seconds | Threshold | Status |
|---|---:|---:|---:|---|
| G0 chain 1 | 711121 | 0.293 | 1.000 | completed |
| G0 chain 2 | 711222 | 0.348 | 1.000 | completed; identical to chain 1 |
| G1 chain 1 | 721121 | 0.453 | 1.000 | completed |
| G2 chain 1 | 731121 | 2.435 | 0.995 | tuning and fit completed |

The G2 tuning smoke evaluated the native grid and selected 0.995. These timings
exclude deterministic data reconstruction/export and are not estimates for a
four-chain 3,000-iteration registry.

| G2 threshold | Validation correlation | Relative correlation |
|---:|---:|---:|
| .995 | .65984 | 1.00000 |
| .990 | .65773 | .99681 |
| .950 | .66519 | 1.00811 |
| .900 | .68057 | 1.03142 |

The native rule retained .995 because the maximum relative improvement did not
exceed 1.25. Tuning was run only in the G2 smoke; there are no four-chain
threshold comparisons because the independent-chain stop rule fired first.

Official outputs exposed posterior mean/SD/last SNP effects, PIPs, `hsq_hist`,
`ssq_hist`, `pi_hist`, residual-variance history, block histories,
component-count and component-genetic-variance histories, final
`sigma_anno`, alpha text traces, conditional annotation probabilities, joint
component probabilities, per-annotation enrichment, and sparse beta histories.
A retained `sigmaSqAlpha` history was not exposed. Alpha files map stick p1,
p2, p3 (and p4 in G2) to columns Intercept, enriched binary, continuous, and
null in that order.

The four-component wrapper emitted benign table-recycling warnings because its
summary labels are hard-coded as NumSnp2:5/Vg2:5 even when only three active
columns exist. Native RDS dimensions, not those recycled labels, are the
authoritative output contract.

## Why scientific parity is unavailable

With default-constructed native streams, fresh processes repeat the same
trajectory. R `set.seed()` and the wrapper’s `seed=` argument cannot establish
independent chains in v0.2.6. Pooling identical histories would make R-hat/ESS,
alpha stability, occupancy stability, and chain-specific PIP agreement
meaningless. The smokes are therefore only API/format checks.

No G0 comparison with `sblr` SBayesR, G1 comparison with learned BED/block
BayesRC, or G2 native-workflow comparison is reported. The official parity
question, official alpha convergence, five-component geometry, residual-
variance effect, and official h2 calibration all remain unresolved.

## Recommendation

Request or identify an official SBayesRC commit/release with a documented
public-to-native seed contract and a regression test showing different seeds
produce different streams while the same seed reproduces. Then rerun the
unchanged 12-coordinate registry and existing export. Do not patch a private
fork and call it official parity; do not change Study 06 or `sblr` first.

Raw source, installed libraries, exports, manifests, and smoke outputs remain
ignored under `results/local/06_annotation_models/gctb_parity/`. No formal
qualification history was rerun and no final benchmark was launched.
