# Study 06 large information-scale feasibility result

## Decision

**LARGE-F6 — technical or design block.** No scientific fit was launched. The
single registered truth and input audits completed, followed only by short
non-inferential API smokes.

The formal small Study 06 qualification remains failed, historical evidence is
unchanged, and the final benchmark remains unauthorized.

## Continuation audit (2026-08-06)

The package-side continuation removed the first blocker: compact integer
component-count, realized-active-count, and stick-count histories now require
about 2.56 MiB per fit rather than a 13.59 GiB full component array. Exact
full-state oracles and tracing-on/off RNG-neutrality tests pass.

The continuation development installation is `sblr 0.2.0` under
`results/local/06_annotation_models/large_feasibility/continuation/devlib`,
with installed-tree SHA-256
`231de56787d101a15b87fd798bb4f247e3de7c9ee0967b04186498bb741fcdca`
and DLL SHA-256
`ae39c01607613db2af89ede5f31ea4063b7958d621ec5bf4e8cd196b3eba5825`.
Local installation does not provide `RemoteSha`; the package base HEAD and
unstaged source diff are recorded separately.

Package source tests passed (`4,688` expectations, no failures; one opt-in
skip). `R CMD check --no-manual --as-cran` completed compilation,
installation, documentation, examples, and compiled-code checks with zero
warnings, but ended with one test error because an existing test sources the
developer-only `tools/study06_partial_exchange_feasibility.R`, which is absent
from the built tarball. The check therefore does not meet the requested
zero-error contract; this packaging-test issue is separate from the B0
scientific execution blocker.

The B0 blocker remains and is now classified **B0-C4 (unresolved)**. The exact
frozen iteration-0 state was captured immediately before its residual draw.
Native and independent block-factor SSE calculations agree within `4.04e-10`:

| Residual audit term | Value |
|---|---:|
| `yy` | 10,223.891131 |
| transformed-score squared norm | 82,949.69 |
| reduced-residual squared norm | 68,502.41 |
| block fitted quadratic | 13,411.44 |
| block SSE | -4,223.39 |
| block residual scale | -4,221.35 |
| direct same-BED quadratic | 30,740.09 |
| direct same-BED SSE | 13,105.26 |
| omitted cross-block quadratic | 17,328.66 |

This is not floating-point cancellation. Full positive retention reproduces
every within-block cross-product but does not restore cross-block terms. A
0.995 diagnostic also failed materially, and `adjE = 0` did not change the
first failure. Clamping, applying a tolerance, changing rank, or adopting a
different residual-variance likelihood would violate the frozen execution
contract. Therefore no continuation smoke or scientific fit was promoted and
the decision remains `LARGE-F6`.

## Frozen identities and deterministic truth

The experiment used `sblrbench` commit
`6ef276eab99a37946901812d5466274c23f56f8a` and clean sibling `sblr` commit
`a165fb0635afcb8a712e8658175dfbb19896b3c3` (package 0.2.0). The exact sibling
source was rebuilt into the established isolated benchmark library; local-source
installation does not populate `RemoteSha`, so provenance is the clean source
HEAD, build command, and installed-tree hash.

| Identity | SHA-256 |
|---|---|
| specification | `b001bc36a5531e5e6b342286a253fc1fd34dad4265359d89d2feaa026d4533df` |
| truth | `e94a511540f600e61ef47b52947836f19a15388f5e8ce795c179929956817507` |
| marker order | `e394f4324f89ca7ad88691284a6216a24e26c691d7cb59d822d5d7f006b096f2` |
| sample order | `f3c14e98845fc565424a973a5a5850768c587f9d5a5306cc58ff36c5329a256b` |
| alleles | `996fb56147c14801f798a7bb21692a962aba18829c8a3a7abfa4dc768b06b082` |
| annotations | `5cc2a5dc64140b5cb2c3e77044d8a0b7bd698e2b33f1f8cb0f5b39ac278f7bd3` |
| blocks | `3b119c38fcababc5b70892f8fc29dcd773fb73e4e2a2b1cd6daf94556e16294d` |
| GWAS | `f00c01326d20d32c1b389f239ebe520327f4b8c859acf7320b10d531243cacfd` |
| alpha truth | `4766d00b77653825e9130a32ebcde1b16754ee99103f2ab4c4d3f1d715fbbf82` |
| selected panel | `0ae8cc37d0418d54cf52e4cf5271c5859d01506759a6d798f6c512a66b08438f` |

QC retained 37,991 of 50,000 markers for all 5,000 unique individuals. There
were no missing dosages among retained markers. MAF ranged from 0.0501 to 0.5.
The physical panel forms 76 blocks: 75 of size 500 and one of size 491.

The annotation prevalence is 0.199995, its binary/continuous correlation is
0.201909, both null correlations are below `1e-17`, matrix rank is four, and the
condition number is 2.44. Every block has enriched and unenriched markers.

Solved intercepts are `-2.471644`, `-0.693114`, and `-0.947510`. Expected
proportions match `0.970/0.015/0.010/0.005` to numerical precision, giving
expected counts `36,851.27 / 569.865 / 379.910 / 189.955`. The one registered
truth draw realized `36,791 / 618 / 392 / 190`, or 1,200 active markers. Stick
eligible/continuation counts were `37,991/1,200`, `1,200/582`, and `582/190`;
all eligible-subset rank and binary-outcome gates passed. Realized genetic and
residual variances are both 1.0 and realized heritability is exactly 0.50.

All 76 block matrices were symmetric, finite, and positive under the registered
tolerance. No positive mode was omitted. Maximum source-LD reconstruction error
was `4.99e-14`, transformed-score error `1.56e-13`, maximum sampled global
cross-block absolute correlation 0.0681, and maximum complete adjacent-block
absolute correlation 0.0791.

## Initial non-inferential smokes and blockers

| Fit | Status | Wall time | Scientific use |
|---|---|---:|---|
| E0 BED BayesR | passed, 12 iterations x 4 chains | 22.63 s | none |
| B0 block-eigen SBayesR | failed at iteration 0 | 83.5 s | none |
| E2 fixed-alpha BED BayesRC | passed, 12 iterations x 4 chains | 11.05 s | none |
| E1 learned BED BayesRC | passed, 12 iterations x 4 chains | 11.34 s | none |
| B2 fixed-alpha block SBayesRC | not launched after route stop | — | none |
| B1 learned block SBayesRC | not launched after route stop | — | none |

E2 held alpha exactly (`max absolute final error = 0`). E1 exposed all 12 alpha
and three `sigmaSqAlpha` histories. All successful smokes exposed the selected
300-marker effect, active, and component histories with finite outputs.

The original `LARGE-F6` phase recorded two independent blockers:

1. The pinned public trace API exposes component-probability histories and
   selected-marker component histories, but no compact genome-wide
   per-iteration component-count or realized-active-count history. Deriving it
   would require tracing all 37,991 component states. The materialized R trace is
   about 14.6 GB per fit, or 87.5 GB across six fits, before checkpoints and other
   histories; only 54.9 GB was free. Final component states cannot substitute for
   posterior traces. The continuation work has now removed this blocker through
   an RNG-neutral compact aggregate trace.
2. B0 stopped before a retained iteration with: `BayesR operator residual scale
   is invalid. trait=0, chain=0, iter=0`. This used the required same-sample,
   all-positive-mode, cumulative-mass operator. Changing eigen retention,
   operator route, seeds, or initialization after the failure was prohibited.

The original successful smokes and the continuation B0 reproducer cannot answer
the scientific question. There are no
convergence, recovery, PIP, ranking, effect, genetic-value, route-agreement, or
heritability-offset results. They must be reported as unavailable, not as
failures of learned-alpha sampling.

## Interpretation and next task

The correct conclusion is technical: the registered six-fit feasibility study
could not be executed under the current public trace and block-eigen contracts.
This is not evidence for or against information-scale-dependent convergence.

The compact-retention work is complete. The next task is a separately reviewed
block-eigen residual-likelihood contract decision. It must choose and validate
how residual variance is defined when block projections overlap in sample
space; this cannot be disguised as a numerical fix. Only after that decision
can this exact frozen experiment—not a new truth seed—be resumed.
