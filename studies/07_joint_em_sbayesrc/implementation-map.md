# Study 07 current `sblr` 0.2.0 implementation map

Study 07 benchmarks clean sibling source SHA
`a85b749dfb903364c5cdb17ef8b8efc1583e347d`, installed as `sblr` 0.2.0.
The installed tree SHA-256 is
`a22927d8d1b6f7c8af3a5b18db5e7a827865832ea1f5bb47f6a4fd21ebb32b65`.

| Study method | Inference target | Route used | API status | B/E semantics |
|---|---|---|---|---|
| SBayesR | joint Bayesian | `stblr_block_eigen(method = "sbayesr")` | public | learned B/E; `gctb_block`, `allMixVe` |
| SBayesRC | joint Bayesian | `stblr_block_eigen(method = "sbayesrc")` | public | learned B/E; `gctb_block`, `allMixVe` |
| SBayesRC-EM | alpha MAP / empirical Bayes | `.stblr_mcem_sbayesrc_block_eigen()` | qualified internal | learned B/E across E-steps; separate final frozen-prior genomic block |
| SBayesRC-S | full joint selection posterior | `.st_sbayesrc_selection_csr()` | qualified internal CSR | learned B/E under the exact Study 06 block-diagonal CSR; global projected legacy residual contract |
| SBayesRC-S-EM | delta/alpha MAP plus conditional EB-PIP | `.stblr_mcem_sbayesrc_s_block_eigen()` | qualified internal | learned B/E; separate final frozen-prior genomic block |

The old repository implementation map described an installed 0.1.2 snapshot:
its public block route reconstructed dense storage and it did not contain the
Phase 5B/5C inference lines. It remains historical context. Current 0.2.0 adds
the retained-low-rank single-trait block backend and the two qualified internal
EM routes above. SBayesRC-EM fixes `sigmaSqAlpha`, has no independent global Pi
update, and reports `method = "SBayesRC-EM"`, `algorithm = "MCEM"`, and
`inference = "mcem"`. SBayesRC-S-EM fixes `pi_A` and `tau2` and reports
`algorithm = "MCEM-Laplace"`; its `annotation_pip_eb` is conditional on final
Rao-Blackwell responsibilities and is not the joint SBayesRC-S posterior PIP.

The narrow internal adapters are Study 07 orchestration only. They do not
redefine or modify any statistical transition in `sblr`.
