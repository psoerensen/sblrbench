# Study 06 GCTB-compatible block contract validation result

## Decision

**SMALL-G1 — contract validation passed.** Both preregistered SNP-inference and
residual-contract gates passed. This validates the new block contract on the
pinned small comparison; it does not alter the historical Study 06
qualification, establish official multichain convergence, or authorize the
final benchmark.

## Identity and execution

The comparison reused the validated 1,500-marker, 15-by-100-block official
export and official D0/D1 artifacts without modification. New S0/S1 artifacts
were written under the ignored
`results/local/06_annotation_models/gctb_block_contract_validation/` root.
The loaded package was `sblr` 0.2.0 from clean source SHA
`0c89234273389e14112ba0e08ef9d11d3e1819dc`; installed-tree SHA-256 was
`e723528e7d5d570a31b5b1d1c90551896ac48f86ab05261c181c8109af971fd0`.

| Comparison | PIP Pearson/Spearman | Effect Pearson | Validation-g Pearson | Top 50/100 | AUPRC new/official | AUROC new/official |
|---|---:|---:|---:|---:|---:|---:|
| S0-new vs D0 | .9955/.8043 | .9984 | .9985 | 49/97 | .3225/.3203 | .6513/.6641 |
| S1-new vs D1 | .9881/.9771 | .9974 | .9973 | 45/90 | .5669/.5572 | .8407/.8418 |

| Comparison | Mean block Ve new/official | Absolute difference | Summary h2 new/official | Absolute difference |
|---|---:|---:|---:|---:|
| S0-new vs D0 | 1.987706/1.987645 | .000060 | .48990/.46499 | .02491 |
| S1-new vs D1 | 1.987641/1.987544 | .000097 | .48651/.46085 | .02565 |

Runtime was 9.88 seconds for S0-new and 18.31 seconds for S1-new. Mean realized
active counts were 177.1 versus 486.4 for S0/D0 and 94.2 versus 78.9 for S1/D1;
these were explicitly nongating sampler/prior differences.

S1-new alpha means by stick (intercept/enriched/continuous/null) were
`-2.406/1.650/.590/-.230`, `.852/.499/.253/.266`, and
`-.007/-.867/.240/.986`. `sigmaSqAlpha` means were
`1.774/1.378/3.967`. These are descriptive single-trajectory summaries, not
official alpha parity or convergence evidence.

## Interpretation

On the exact small data, the new `gctb_block`/`allMixVe` implementation agrees
with the official D0/D1 artifacts on residual scale and on the primary SNP
outputs. This cleared the preregistered gate for resuming the separately frozen
large experiment. It does not relabel historical `global_projected` fits.
