# Study 02: single-trait prediction development benchmark

Study 02 compares BayesC and BayesR through individual-level BED and
summary-statistics/sparse-LD CSR representations. Two architectures use the
same 50 causal markers per replicate and target h2 = 0.30: approximately
BayesC-like homogeneous effects and approximately BayesR-like variance-mixture
effects. Final heritability rescaling means neither is an exact sampler-prior
draw.

One deterministic split supplies 3,500 training and 1,500 test samples.
Training individuals alone determine allele frequencies, genotype scaling,
sparse LD and summary statistics. Test phenotypes appear only in evaluation.

Run `SBLR_BENCH_REPLICATES=1`, then expand to 5 or 10 without deleting the
targets store. The 500-iteration, 250-burn-in, one-chain settings validate
infrastructure only and do not establish convergence or method superiority.
