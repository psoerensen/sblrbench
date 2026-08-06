# Study 06 documentation cleanup result

## Decision

**DOC-C1 — cleanup complete.** Study 06 now has one readable report, one
compact navigation page, and one authoritative chronological synthesis. Formal
decisions and experiment-specific result records remain separate and unchanged.
Every benchmark-side experiment and the relevant package-side correctness,
correction, and transition audits are discoverable through the
[inventory](study06_documentation_inventory.md).

## Provenance and scope

The documentation-only task started from clean `sblrbench` `master` at
`d9ecfe20368daff2eec2a01a73c7a13b2dfa8f31` and clean read-only sibling
`sblr` `master` at `a165fb0635afcb8a712e8658175dfbb19896b3c3`
(`sblr` 0.2.0). HEAD did not change. No MCMC, qualification, official
SBayesRC trajectory, benchmark, larger simulation, package build, install, or
compilation was run. No files were staged, committed, or pushed.

Immutable Study 06 identities:

```text
Specification:      241c15afab8fefc571e38e625130de6e4ab58b958c30cff369e26998ee30fa56
Informative truth:  169d52bef390022a9106d7e61b200493869b40cbb76f1c8d911eebbc80fea1eb
```

## Documentation hierarchy established

1. `studies/06_annotation_models/report.qmd` — primary user-facing report.
2. `studies/06_annotation_models/README.md` — compact status/index and reading
   paths.
3. `docs/dev/study06_annotation_inference_evidence_synthesis.md` — authoritative
   chronological technical synthesis.
4. Experiment result documents and decision JSON — authoritative within their
   own formal or descriptive scope.
5. Design/protocol and package mathematical/operator notes — supporting
   evidence.
6. V1 and superseded development records — explicitly historical and retained.

The principal documents now share this status:

```text
v1 sparse qualification: failed and preserved
v2 identifiable qualification: failed
paired power isolation: completed
package-side hierarchy and transition audits: completed
official SBayesRC multichain parity: blocked by the v0.2.6 seed contract
official SBayesRC single-trajectory descriptive comparison: completed
final benchmark: not authorized
larger n=5000, m≈38000 feasibility experiment: planned, not yet run
```

## Cleanup actions

- Reorganized the main report around purpose, designs, formal decisions,
  explicit evidence levels, conclusions, implementation lessons, official
  comparison, alpha/variance interpretation, runtime, timeline, and roadmap.
- Rebuilt the Study 06 README as the navigation authority with four recommended
  reading paths and reproducibility links.
- Extended the technical synthesis without replacing its ledger or evidence
  map. Added package lessons, official quantitative parity, alpha/variance
  crosswalk, runtime context, and the current roadmap.
- Updated the repository README, study README, and website catalogue with the
  multichain-versus-single-trajectory distinction.
- Added a complete tracked-file inventory and machine cleanup decision.
- Added focused documentation-contract tests. No simulation, sampler, fit,
  threshold, evidence, or result logic changed.

No experiment-specific result or design document was edited. In particular,
the v1/v2 qualification decisions, paired-isolation decision, GCTB-P5 result
and decision, and GCTB single-trajectory result and decision remain byte-for-
byte unchanged.

## Scientific consolidation

Established findings now clearly separate converged/exact evidence from the
delimited official descriptive comparison: baseline and fixed-alpha convergence,
conditional correctness, frozen-allocation convergence, failure of variance
and prior ablations to solve dynamic mixing, H20/A20 schedule lessons, failure
of the tested coupling ladder, persistent route calibration, annotation ranking
benefit, and strong matched official/`sblr` SNP-level agreement.

Descriptive findings are labelled as such: learned-alpha ranking/prediction,
pooled learned summaries, official alpha/architecture, within-path drift,
native five-component behavior, per-chain implementation agreement, and
runtime. Unresolved items include official multichain convergence,
alpha/annotation-variance and latent-architecture parity, residual/effect-scale
contracts, active/component-count semantics, BED/block calibration, and the
information-scale versus transition-design question.

## Implementation lessons represented

The report and synthesis cover the shared annotation kernel, stick orientation,
Albert–Chib and alpha conditionals, `sigmaSqAlpha`, component probabilities,
the corrected BED trace-index defect, frozen-allocation hierarchy, fixed-
variance and production-prior ablations, S1/H5/H20/A5/A20, tiny/exact and
Study 06 tempering, partial-exchange state limitations, block runtime/storage,
and unresolved official/`sblr` contracts. Package-side sources remain linked
and were not copied wholesale.

## Current roadmap

The next scientific phase is a planned, not-yet-run one-replicate feasibility
experiment with `n = 5000`, `m ≈ 38000`, four independent `sblr` chains, and
complete alpha, `sigmaSqAlpha`, occupancy, expected/realized active-count,
variance, heritability, PIP, and truth-recovery checks. Replication follows only
after success. A fixed-state official/`sblr` latent-contract audit is optional
and independent. Official multichain parity remains dependent on a seedable
official release or documented independent-chain mechanism.

## Preservation and validation

Before editing, the standalone Study 06 report and complete Quarto project both
rendered cleanly. After editing, changed R/JSON parsing, focused and complete
tests, both renders, link/heading checks, `git diff --check`, and immutable hash
comparisons were run. Exact counts and final status are reported in the task
handoff; all checks passed with zero failures or errors.

The 19 v1 `current-stop` files, v2 qualification evidence, paired-isolation
evidence, GCTB-P5 result/decision, single-trajectory result/decision and raw
manifest identity, and referenced sibling evidence retained their before/after
SHA-256 values. The sibling working tree remained clean and unchanged.
