test_that("simulate_assay returns the requested number of samples", {
  parameters <- default_parameters()
  parameters$n_samples <- 100L
  parameters$seed <- 1L

  simulated <- simulate_assay(parameters)

  expect_s3_class(simulated, "viralquantsim_data")
  expect_equal(nrow(simulated), 100L)
  expect_true(all(grepl("^SYN", simulated$sample_id)))
})

test_that("simulation is reproducible for a fixed seed", {
  parameters <- default_parameters()
  parameters$n_samples <- 50L
  parameters$seed <- 99L

  first <- simulate_assay(parameters)
  second <- simulate_assay(parameters)

  expect_equal(first$true_log10_viral_load, second$true_log10_viral_load)
  expect_equal(first$target_ct, second$target_ct)
  expect_equal(first$target_reads, second$target_reads)
})

test_that("sequencing read counts respect total reads", {
  parameters <- default_parameters()
  parameters$n_samples <- 100L

  simulated <- simulate_assay(parameters)

  expect_true(all(simulated$target_reads >= 0))
  expect_true(all(simulated$target_reads <= simulated$total_reads))
})

test_that("Delta Ct follows the documented direction", {
  parameters <- default_parameters()
  parameters$n_samples <- 100L

  simulated <- simulate_assay(parameters)
  expected <- simulated$target_ct - simulated$internal_control_ct

  expect_equal(simulated$delta_ct, expected)
})
