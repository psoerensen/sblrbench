# Contract smoke study

This deterministic infrastructure test creates a small named genotype matrix, calls installed `sblr::mtsim()`, validates `Z %*% B`, runs an oracle/mock method, and writes only compact metrics and provenance. It performs no MCMC and needs no external files or sibling source access.
