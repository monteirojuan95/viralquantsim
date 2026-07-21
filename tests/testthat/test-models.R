test_that("candidate models fit and comparison is tidy", {
  parameters <- default_parameters()
  parameters$n_samples <- 150L
  parameters$seed <- 7L

  simulated <- simulate_assay(parameters)
  models <- fit_candidate_models(simulated)
  comparison <- compare_models(
    models,
    simulated,
    cross_validate = FALSE
  )

  expect_named(
    models,
    c("target_ct", "delta_ct", "reads", "reads_target_ct", "reads_delta_ct")
  )
  expect_equal(nrow(comparison), 5L)
  expect_true(all(c("model", "rmse", "mae", "bias") %in% names(comparison)))
})
