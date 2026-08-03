# Supplemental Study 06: SBayesR LD-operator robustness

This focused package-validation analysis supplements Study 06. It consolidates
the BED-scheduler and exact-versus-sparse diagnostics and extends the same
LD-rich window with full-rank and retained-low-rank block-eigen SBayesR fits.
It is not a general comparison of all `sblr` methods.

Run from the repository root:

```powershell
Rscript studies/06_ld_operator/sbayesr_ld_robustness/analysis.R
quarto render studies/06_ld_operator/sbayesr_ld_robustness/report.qmd
```

`analysis.R` validates and reuses every completed scheduler, BED, exact-CSR,
hard-sparse CSR, and block-eigen checkpoint. Large inputs and fits remain
ignored under `results/local/`; the report reads exclusively from the compact
capsule at `results/reference/06_ld_operator/sbayesr_ld_robustness/current/`.

The scripts retained under `scripts/` reproduce the two precursor diagnostics.
They continue to use their original local checkpoint directories so moving the
source does not invalidate completed fits.
