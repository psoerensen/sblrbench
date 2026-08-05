# sblrbench

`sblrbench` has two related purposes: it validates the scientific and software
behavior of [`sblr`](https://github.com/psoerensen/sblr), and it provides clear,
extensible analysis workflows built on the same validated framework.

## Validate `sblr`

Five completed scalar-validation studies are available:

1. Study 01 — Fine-mapping
2. Study 02 — Prediction
3. Study 03 — Parameter estimation
4. Study 04 — Convergence
5. Study 05 — LD-operator validation

Each completed report reads only a compact frozen capsule under
`results/reference/`. Capsules contain numerical tables, manifests, checksums,
source inventories, and provenance—not native fit objects or local checkpoints.
The reports document the limits of each design; these focused simulations do
not establish universal rankings or production defaults.

New validation studies should keep scientific choices in an ordinary-list
`spec.R`, use shared mechanics under `R/`, expose a readable `analysis.R`, and
publish only reviewed compact results through a frozen capsule.

Study 06 — Annotation-informed models is **in development**: its v1 sparse
qualification failed and is preserved, while its v2 identifiable qualification
completed all four fits but failed convergence, scientific-recovery, and route
gates; the final benchmark is unauthorized. Study 07 —
Multitrait validation is **in development**. Neither is a completed validation
benchmark.

## Practical analysis workflows

The exact study scripts under `studies/01_finemapping/` through
`studies/05_ld_operator/` show complete auditable workflows. Shorter reusable
examples live under `inst/templates/` for fine-mapping, prediction, parameter
estimation, convergence, and operator analysis.

The shared runner uses ordinary R functions and lists:

```r
library(sblrbench)

spec <- read_benchmark_spec("studies/02_prediction/spec.R")
result <- run_benchmark(
  spec,
  output_dir = "results/local/my-prediction-run",
  profile = "workshop",
  resume = TRUE
)
```

`workshop` profiles are for learning and workflow checks, not performance
claims. `benchmark` profiles retain the validated scientific coordinates. Local
outputs use predictable `checkpoints/`, `tables/`, `figures/`, manifest, and
session-information paths. Scenarios, supported methods, metrics, and output
locations can be adjusted through the ordinary specification and template
code; no plugin or workflow language is required.

The command-line entry point supports completed Studies 01–05:

```powershell
Rscript scripts/run_benchmark.R `
  --study 02_prediction `
  --profile workshop `
  --output-dir results/local/example `
  --resume true `
  --validate-only true
```

Validation-only mode resolves and checks the design without fitting.

## Repository roles

- `bayesian-linear-regression`: course and teaching material.
- `sblr`: software implementation and authoritative implemented methodology.
- `sblrbench`: validation benchmarks and extensible practical
  analysis/reporting workflows.

## Reproducibility and website

The benchmark data source, package revisions, seeds, controls, capsule
checksums, and report dependencies are recorded in the specifications and
reference capsules. Working outputs remain ignored under `results/local/`.

Render the capsule-only website with:

```powershell
quarto render
```

The published site is <https://psoerensen.github.io/sblrbench/>. See the
[benchmark catalogue](studies/index.qmd), [practical workflows](workflows.qmd),
[framework](framework.qmd), [metrics](metrics.qmd), and
[reproducibility guide](reproducibility.qmd).
