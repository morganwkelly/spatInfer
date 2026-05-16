# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## About

`spatInfer` is an R package implementing diagnostics and regression procedures from Conley and Kelly (2024) to estimate reliable results from regressions using spatial observations. It addresses inflated t-statistics caused by strong autocorrelation and directional trends in spatial data.

## Common Commands

Standard R package development workflow:

```r
# Install dependencies and load package
devtools::install_deps()
devtools::load_all()

# Build documentation from roxygen2 comments
devtools::document()

# Run R CMD CHECK
devtools::check()

# Run tests (if any exist)
devtools::test()

# Install from local source
devtools::install()

# Install from GitHub
pak::pak("morganwkelly/spatInfer")
```

Note: `README.md` is generated from `README.Rmd` — edit the `.Rmd` source, not the `.md`.

## Architecture

The package exposes a **sequential 4-step workflow**:

1. `optimal_basis()` — selects optimal spatial basis via BIC minimization over tensor spline dimensions (3×3 to N×N using `mgcv::bam`)
2. `placebo()` / `placebo_im()` — placebo tests using spatially-correlated noise simulations
3. `synth()` / `synth_im()` — synthetic outcome tests (null: outcome is trending spatial noise)
4. `basis_regression()` — final regression with spatial basis and large-cluster standard errors (returns a `fixest` object)

Two inference methods are available: **BCH** (bias-corrected heteroskedasticity-consistent) and **IM** (inverse measure). Functions ending in `_im` use IM inference; the base versions use BCH.

### Key Source Files

- **`R/spatial_helpers.R`** — core infrastructure: cluster generation (`generate_clusters`), principal components from tensor basis (`prin_comp`), Matérn spatial noise simulation (`Noise_Sim`), Moran test (`moran`), p-value computations
- **`R/bch_helpers.R`** — BCH inference: single simulation (`bch_sim`), parallel p-value computation (`bch_p_values`), result summary
- **`R/im_helpers.R`** — IM inference: analogous functions plus `coef_im()` for coefficient extraction
- **`R/optimal_basis.R`**, **`R/placebo.R`**, **`R/placebo_im.R`**, **`R/synth.R`**, **`R/synth_im.R`**, **`R/basis_regression.R`** — one file per exported function
- **`R/placebo_table.R`**, **`R/synth_table.R`** — format results to `tinytable` output

### Statistical Approach

- Spatial basis: tensor product splines (`bs="bs"`) via `mgcv`
- Correlation modeling: Matérn process with kriging (`fields`)
- Clustering: k-medoids (PAM) or Clara for large datasets (`cluster`, `spdep`)
- Parallelism: `foreach`/`doParallel` for simulation loops
- Final regression uses `fixest::feols` with cluster-robust SEs

### Data

`data/opportunity.rda` — Chetty et al. (2014) intergenerational mobility data; 693 US Census zones with variables: `mobility`, `single_mom`, `gini`, `dropout_rate`, `social_cap`, `X` (longitude), `Y` (latitude). Used in all package examples.

## Tutorial

Full step-by-step tutorial: https://morganwkelly.github.io/spatInfer_tutor/
