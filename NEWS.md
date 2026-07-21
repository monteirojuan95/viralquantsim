# viralquantsim 0.1.0

First public release.

* Mechanistic simulation of a known true target load, with technical variation
  decomposed into shared and analyte-specific recovery, reverse-transcription,
  and inhibition effects (`simulate_assay()`).
* qPCR channel with technical replicates, Ct censoring at the detection limit,
  and Delta Ct (`simulate_qpcr()`).
* Sequencing channel with overdispersed target read counts and detection-limit
  dropout (`simulate_sequencing()`).
* Five interpretable candidate viral-load models with repeated V-fold
  cross-validation (`fit_candidate_models()`, `compare_models()`).
* Two-parameter sensitivity grid quantifying when internal-control
  normalisation helps (`run_sensitivity_grid()`).
* Shiny interface (`run_app()`).
* Sequencing-noise defaults set to a generic quantitative-sequencing regime:
  moderately strong read-based quantification (reads-only R-squared around 0.8),
  a detection limit near 100 copies/mL, and near-Poisson count sampling.
