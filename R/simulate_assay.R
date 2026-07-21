#' Simulate an integrated qPCR and sequencing assay
#'
#' Creates a mechanistic synthetic dataset in which the true target load is
#' known. Shared technical effects affect both the target and the spike-in
#' internal control, while analyte-specific effects prevent normalisation from
#' being unrealistically perfect.
#'
#' @param parameters A named list. Begin with [default_parameters()] and modify
#'   only the values required for a scenario.
#'
#' @return A data frame with one row per synthetic sample.
#' @examples
#' params <- default_parameters()
#' params$n_samples <- 200
#' simulated <- simulate_assay(params)
#' head(simulated[, c("true_log10_viral_load", "target_ct", "delta_ct")])
#' @export
simulate_assay <- function(parameters = default_parameters()) {
  .validate_parameters(parameters)
  set.seed(as.integer(parameters$seed))

  n <- as.integer(parameters$n_samples)

  data <- data.frame(
    sample_id = sprintf("SYN%05d", seq_len(n)),
    true_log10_viral_load = stats::runif(
      n,
      min = parameters$viral_load_min,
      max = parameters$viral_load_max
    ),
    stringsAsFactors = FALSE
  )

  data$true_viral_load_copies_ml <- 10^data$true_log10_viral_load
  data$nominal_ic_log10_input <- parameters$ic_log10_input
  data$actual_ic_log10_input <-
    data$nominal_ic_log10_input +
    stats::rnorm(n, mean = 0, sd = parameters$ic_pipetting_sd)

  data$shared_recovery_effect <- stats::rnorm(
    n,
    mean = 0,
    sd = parameters$shared_recovery_sd
  )

  data$target_recovery_effect <-
    data$shared_recovery_effect +
    stats::rnorm(
      n,
      mean = 0,
      sd = parameters$target_specific_recovery_sd
    )

  data$ic_recovery_effect <-
    data$shared_recovery_effect +
    stats::rnorm(
      n,
      mean = 0,
      sd = parameters$ic_specific_recovery_sd
    )

  data$shared_rt_effect <- stats::rnorm(
    n,
    mean = 0,
    sd = parameters$shared_rt_sd
  )

  data$target_rt_effect <-
    data$shared_rt_effect +
    stats::rnorm(
      n,
      mean = 0,
      sd = parameters$target_specific_rt_sd
    )

  data$ic_rt_effect <-
    data$shared_rt_effect +
    stats::rnorm(
      n,
      mean = 0,
      sd = parameters$ic_specific_rt_sd
    )

  data$target_log10_available <-
    data$true_log10_viral_load +
    data$target_recovery_effect +
    data$target_rt_effect

  data$ic_log10_available <-
    data$actual_ic_log10_input +
    data$ic_recovery_effect +
    data$ic_rt_effect

  data <- simulate_qpcr(data, parameters)
  data <- simulate_sequencing(data, parameters)

  class(data) <- c("viralquantsim_data", class(data))
  attr(data, "parameters") <- parameters
  data
}
