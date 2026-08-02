# Current benchmark refresh

This run replaces versioned development-era benchmark presentations with one current,
validated benchmark per supported study. Git history remains the archive for superseded
capsules and reports.

## Frozen provenance

- Starting `sblrbench` commit: `0fffcb87d0c3e85d7685803caa2299cfcd8409fd`.
- Frozen `sblr` commit: `02e8c74baa906e83c4a08d42a9cc6339b4e81072` (`Align BLR prior variance calibration`).
- Installed package: `sblr` 0.2.0 in `results/local/current_benchmark_refresh/rlib`.
- qgdata commit: `6cca5819e711d326cfb2614d7e9d9f34942612cd`.

The refresh runner prepends the isolated package library to the ordinary R library search
path so installed dependencies remain available. Benchmark fits never use
`pkgload::load_all()` for `sblr`.

## Scientific policy

Study 04 is rerun first and produces the current method-specific MCMC recommendations.
Studies 02 and 03 use those frozen recommendations. Study 01 retains its established
ten-replicate fine-mapping design and its explicit development MCMC contract unless the
Study 04 result and a focused implementation review support a common replacement before
dispatch. Study 05 reruns only its prespecified maximum-history annotation convergence
pilot. Study 06 numerical evidence is reusable only if resolved prior inputs and all
sampler-transition source blobs are unchanged relative to its original package revision.

## Isolation and resume

All refresh outputs, target stores, logs, and the package library live below
`results/local/current_benchmark_refresh/`. A target store is unique to each phase.
The runner verifies the exact package SHA before each phase. Promotion occurs only after
coordinate, chain, diagnostic, provenance, and checksum validation.

The tracked cleanup inventory is generated locally at
`results/local/current_benchmark_refresh/audit/benchmark_inventory.csv`. No superseded
tracked capsule is removed until its current replacement and report validate.
