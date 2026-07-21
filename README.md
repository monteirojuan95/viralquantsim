# ViralQuantSim

[![R-CMD-check](https://github.com/monteirojuan95/viralquantsim/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/monteirojuan95/viralquantsim/actions/workflows/R-CMD-check.yaml)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE.md)

**An R framework for benchmarking internal-control normalisation in quantitative qPCR and sequencing workflows.**

Repository: <https://github.com/monteirojuan95/viralquantsim>

ViralQuantSim simulates a mechanistic assay in which the true target load is known. It separates shared technical variation from target- and internal-control-specific noise, generates qPCR measurements and sequencing read counts, and compares interpretable estimators of target concentration.

## Why this project exists

Paired sample-level datasets containing a reference pathogen load, pathogen Ct, spike-in internal-control Ct, and sequencing metrics are rarely released publicly. ViralQuantSim does not fabricate a supposedly real replacement dataset. Instead, it makes the data-generating assumptions explicit and asks a methodological question:

> Under which technical conditions does an internal control improve, not change, or worsen quantitative estimation?

No empirical dataset is included or reconstructed; every value is simulated from explicit, documented parameters.

## Current milestone

Version 0.1.0 contains:

- mechanistic simulation of true target load and fixed-dose internal-control input;
- shared and analyte-specific recovery and reverse-transcription effects;
- qPCR technical replicates, Ct censoring, and Delta Ct;
- overdispersed sequencing-derived target read counts and dropout;
- five candidate linear models;
- repeated V-fold cross-validation;
- two-parameter sensitivity analysis;
- an initial Shiny interface.

## Installation

Install the development release from GitHub:

```r
install.packages("remotes")
remotes::install_github("monteirojuan95/viralquantsim")
```

Optional interface packages:

```r
install.packages(c("shiny", "bslib", "ggplot2", "plotly", "DT"))
viralquantsim::run_app()
```

For local development:

```bash
git clone https://github.com/monteirojuan95/viralquantsim.git
cd viralquantsim
```

```r
install.packages(c("devtools", "testthat"))
devtools::load_all()
devtools::test()
devtools::check()
```

## Minimal example

```r
library(viralquantsim)

parameters <- default_parameters()
parameters$n_samples <- 1000
parameters$shared_recovery_sd <- 0.40
parameters$ic_specific_recovery_sd <- 0.10

simulated <- simulate_assay(parameters)
models <- fit_candidate_models(simulated)
compare_models(models, simulated)
```

## Candidate models

The initial benchmark compares:

1. true load ~ target Ct;
2. true load ~ Delta Ct;
3. true load ~ log10 target reads;
4. true load ~ log10 target reads + target Ct;
5. true load ~ log10 target reads + Delta Ct.

Delta Ct is defined as:

```text
Delta Ct = target Ct - internal-control Ct
```

## Interpretation of the sensitivity metric

`run_sensitivity_grid()` reports:

```text
delta_rmse = RMSE(reads + Delta Ct) - RMSE(reads only)
```

- Negative: the internal-control-adjusted model performed better.
- Near zero: the internal control added little predictive value.
- Positive: adding Delta Ct worsened prediction.

## Scientific boundaries

ViralQuantSim is a simulation and methodological benchmarking tool. It is not:

- a clinical viral-load calculator;
- a validated diagnostic device;
- a substitute for wet-laboratory validation;
- a reconstruction of a private clinical dataset.

## Roadmap

The next milestones are:

- migrate and modernise the original model-diagnostic logic;
- add influence, heteroscedasticity, calibration, and bootstrap modules;
- build a 3D Plotly regression surface;
- add schema-validated upload of user-owned data;
- expand the Shiny sensitivity explorer;
- create a Quarto documentation website.

## Citation

After installation, run:

```r
citation("viralquantsim")
```

GitHub also exposes citation metadata from the root `CITATION.cff` file.

## Licence

Code is released under the MIT License.
