# Study 06: current retained low-rank LD operator benchmark

This current capsule contains 60 validated four-chain fits across two
architectures, five paired replicates, and six operator configurations.
The numerical fits retain their original sblr source SHA
`bd8e2c8148a0d9540dc20716455706beeb16fa86`.

Compatibility with the current refresh SHA `02e8c74baa906e83c4a08d42a9cc6339b4e81072` was established without rerunning samplers: resolved BayesC and BayesR priors are identical, and all compiled scalar sampler and retained-low- rank transition sources are byte-identical. See `compatibility_manifest.json`, `prior_compatibility.csv`, and `source_compatibility.csv`.

Native fits, checkpoints, Q matrices, and target stores are excluded.
