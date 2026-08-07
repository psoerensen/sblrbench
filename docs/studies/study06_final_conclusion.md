# Study 06 final conclusion

## Primary decision

**EST-R2 — probability/ranking functions are stable but annotation contrasts
remain uncertain.** Official replication is secondarily **EST-R5-inconclusive**
because pinned SBayesRC v0.2.6 cannot supply independent native trajectories.

## Validated

- BayesR/SBayesR baseline controls and fixed-alpha BayesRC/SBayesRC controls.
- The retained block likelihood and residual-variance contract.
- Strong small-study official agreement for block residual variance, PIPs,
  effects, genetic values, and annotation benefit.
- Exactness of the audited hierarchy transitions, including PMA-R3 as an exact
  but computationally impractical global reference.
- Stable SNP-level annotation benefit and rankings on the frozen fixture.

## Difficult

- Unrestricted continuous-alpha joint posterior learning.
- Global alpha/allocation coupling and efficient joint movement.
- Quantitative recovery of the exact continuous coefficient decomposition and
  annotation-contrast magnitudes.

These observations support weak identification and poor practical exploration
language. They do not establish mathematical non-identifiability.

## Scientific stability

| Level | Study 06 conclusion |
|---|---|
| Raw alpha | Marginal R-hat up to 1.025/1.040, low ESS 70–81, relative MCSE up to 0.13, and material truth error; coefficient decomposition is not a reliable scientific endpoint. |
| A alpha | Representative scalar R-hat improves to at most 1.014, but whole-vector minimum correlations of 0.691/0.862 and large truth RMSE show residual instability. |
| Annotation probabilities | Rankings are stable: active-prior Spearman at least 0.986/0.950 and largest-component Spearman at least 0.999/0.995. Probability magnitudes and expected active count are less reliable. |
| Annotation contrasts | Directions are stable for enriched and continuous annotations and null intervals cover zero, but magnitudes are biased, route-sensitive, and incompletely covered. |
| SNP PIPs/effects | Highly stable: PIP Pearson at least 0.994/0.969 and signed beta Pearson at least 0.9999/0.9991. |
| Prediction/validation | Stable across chains: genetic-value agreement at least 0.99994/0.99925; validation correlations vary only in the third decimal place. |

Values before/after the slash are BED/block.

## Evidence types

Source-derived evidence comprises the frozen Study 06 identities, established
control and likelihood conclusions, official adapter provenance, official RNG
audit, and PMA-R3 endpoint. New evidence comprises the full 4 × 9,000-draw BED
and block transformations, counterfactual effects, posterior-direction audit,
rank comparisons, SNP stability metrics, and offline official D1 transformations.

Interpretation is deliberately narrower than the computation: stable rankings
support SNP-level utility; directional contrasts support qualitative biology;
neither licenses precise continuous annotation-effect claims. The official arm
cannot independently confirm multitrajectory convergence.

## Closure and next study

Study 06 is closed. Further unrestricted continuous-alpha sampler development
is not the next task. The recommended next methodology is a distinct model:
Bayesian annotation selection with annotation posterior inclusion probabilities,
targeting `P(delta_j = 1 | y)` rather than a continuous coefficient for every
annotation.

That model is motivated by the scientific estimand, not merely computation:
annotations improve SNP inference, the exact continuous decomposition is
difficult, and lower-dimensional relevance/directional summaries are more
stable and interpretable. The conceptual connection to Bayesian MAGMA and
Bayesian PoPS should inform the design, but neither model is implemented here.

## Reproduction

Run `Rscript scripts/run_study06_estimability.R` from the repository root with
the frozen ignored Study 06 artifacts present. The script consumes retained
chains and does not invoke a sampler. Compact outputs, hashes, probability
guards, and figures are written under
`results/local/06_annotation_models/estimability_and_contrasts/`.

STUDY 06 ESTIMABILITY ANALYSIS COMPLETE — REVIEW REQUIRED
