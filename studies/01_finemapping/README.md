# Fine-mapping pilot skeleton

Planned pilot: one chromosome, one trait; separated-causal and high-LD-causal architectures; 10 replicates; ST-BED BayesC, ST-BED BayesR, ST-CSR SBayesC, and ST-CSR SBayesR. The initial metric is PIP Brier score. Future metrics are average precision, causal-marker rank, credible-set coverage, and credible-set size.

Credible sets must use exported public `sblr` credible-set functions; `sblrbench` will not redefine their semantics. The first real study must verify marker scaling and ordering before scoring. This study may not make implementation changes in `sblr`.
