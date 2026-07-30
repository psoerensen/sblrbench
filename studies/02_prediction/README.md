# Study 02: paired ST/MT prediction

This development study asks how ST and MT BayesR predict genetic values in held-out individuals when effects are mostly shared or mostly trait-specific. It compares BED and training-summary/CSR representations with one deterministic 70/30 split.

The central leakage rule is strict: training individuals alone determine allele frequencies, genotype scaling, sparse LD and summary statistics. Test genotypes are transformed with the frozen training scale; test phenotypes appear only in evaluation.

Native fits remain in the ignored targets store. Depending on the public method contract they retain posterior/component quantities such as `vgs`, `vbs`, `ves`, `pis`, `vle`, and `vld`; Study 02 does not score or interpret them. Dedicated parameter-estimation and convergence studies will handle those quantities separately.

Run one, five or ten replicates with `SBLR_BENCH_STUDY=02_prediction` and `SBLR_BENCH_REPLICATES=1`, `5`, or `10`. Development controls use 500 iterations, 250 burn-in, one chain and one core. They do not establish convergence or method superiority.
