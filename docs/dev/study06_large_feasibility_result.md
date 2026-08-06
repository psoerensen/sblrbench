# Study 06 large information-scale feasibility result

## Decision

**LARGE-F6 — technical or design block.** No scientific fit was launched. The
single registered truth and input audits completed, followed only by short
non-inferential API smokes.

The formal small Study 06 qualification remains failed, historical evidence is
unchanged, and the final benchmark remains unauthorized.

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

## Non-inferential smokes and blockers

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

Two independent blockers prevent contract-compliant scientific execution:

1. The pinned public trace API exposes component-probability histories and
   selected-marker component histories, but no compact genome-wide
   per-iteration component-count or realized-active-count history. Deriving it
   would require tracing all 37,991 component states. The materialized R trace is
   about 14.6 GB per fit, or 87.5 GB across six fits, before checkpoints and other
   histories; only 54.9 GB was free. Modifying `sblr` to add an aggregate trace is
   outside scope. Final component states cannot substitute for posterior traces.
2. B0 stopped before a retained iteration with: `BayesR operator residual scale
   is invalid. trait=0, chain=0, iter=0`. This used the required same-sample,
   all-positive-mode, cumulative-mass operator. Changing eigen retention,
   operator route, seeds, or initialization after the failure was prohibited.

The successful smokes cannot answer the scientific question. There are no
convergence, recovery, PIP, ranking, effect, genetic-value, route-agreement, or
heritability-offset results. They must be reported as unavailable, not as
failures of learned-alpha sampling.

## Interpretation and next task

The correct conclusion is technical: the registered six-fit feasibility study
could not be executed under the current public trace and block-eigen contracts.
This is not evidence for or against information-scale-dependent convergence.

The next task should be a package-side, fixed-input audit with two narrowly
scoped goals: (1) expose compact per-iteration component-count and active-count
aggregates without retaining all marker states; and (2) reproduce and diagnose
the full-positive-mode B0 iteration-0 residual-scale rejection on the frozen
large-feasibility GWAS/block inputs. No sampler transition, prior, seed, or
scientific design change is justified by this blocked run. After those contracts
are corrected and independently tested, this exact frozen experiment—not a new
truth seed—can be considered for a separately authorized execution.

