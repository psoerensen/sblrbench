# Study 05: integrated LD-operator validation capsule

This is the single authoritative compact capsule for Study 05. It combines the completed
main retained-low-rank/operator benchmark with the deterministic SBayesR LD-sensitivity
evidence. Files prefixed `sbayesr_` come from the former supplemental capsule.

The main component contains 60 validated four-chain fits across two architectures, five
paired replicates, and six operator configurations. The SBayesR component contains the
scheduler comparison, exact and hard-sparse CSR comparisons, full-rank and retained-0.995
block-eigen comparisons, corrected-score, quadratic/residual, spectral, and recovery audits.

No sampler was called during integration. Numerical tables were moved byte-for-byte except
for explicit study/path metadata. Native fits, genotype/LD/operator matrices, eigenvectors,
checkpoints, and logs are excluded. Historical source paths are retained only in the clearly
named provenance inventory `historical_operator_source_files.csv`.

Numerical fits retain their original source provenance. Compatibility with sblr SHA
`02e8c74baa906e83c4a08d42a9cc6339b4e81072` was validated without sampler reruns.
