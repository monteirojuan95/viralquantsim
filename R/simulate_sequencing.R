#' Simulate sequencing-derived target reads
#'
#' Generates total reads, a sample-specific target fraction, target read
#' counts, and stochastic dropout. A beta-binomial construction is used to
#' represent overdispersion while ensuring target reads cannot exceed total
#' reads.
#'
#' @param data A data frame containing `target_log10_available`.
#' @param parameters A complete parameter list.
#'
#' @return The input data frame with sequencing columns.
#' @examples
#' input <- data.frame(target_log10_available = c(3, 4, 5))
#' seq <- simulate_sequencing(input)
#' seq[, c("total_reads", "target_reads", "log10_target_reads")]
#' @export
simulate_sequencing <- function(data, parameters = default_parameters()) {
  .validate_parameters(parameters)

  if (!"target_log10_available" %in% names(data)) {
    stop("`data` must contain `target_log10_available`.", call. = FALSE)
  }

  n <- nrow(data)

  shared_library_effect <- stats::rnorm(
    n,
    mean = 0,
    sd = parameters$library_shared_sd
  )

  target_library_effect <- shared_library_effect + stats::rnorm(
    n,
    mean = 0,
    sd = parameters$target_library_specific_sd
  )

  target_log10_library <- data$target_log10_available + target_library_effect

  total_reads <- pmax(
    1L,
    as.integer(round(stats::rlnorm(
      n,
      meanlog = log(parameters$mean_total_reads),
      sdlog = parameters$total_reads_sdlog
    )))
  )

  mean_target_fraction <- stats::plogis(
    parameters$target_fraction_intercept +
      parameters$target_fraction_slope * target_log10_library
  )

  epsilon <- .Machine$double.eps^0.5
  mean_target_fraction <- pmin(
    pmax(mean_target_fraction, epsilon),
    1 - epsilon
  )

  concentration <- parameters$sequencing_concentration
  alpha <- mean_target_fraction * concentration
  beta <- (1 - mean_target_fraction) * concentration
  realised_target_fraction <- stats::rbeta(n, alpha, beta)

  target_reads <- stats::rbinom(
    n,
    size = total_reads,
    prob = realised_target_fraction
  )

  dropout_probability <- stats::plogis(
    parameters$dropout_steepness *
      (parameters$dropout_midpoint_log10 - target_log10_library)
  )

  stochastic_dropout <- stats::rbinom(
    n,
    size = 1L,
    prob = dropout_probability
  ) == 1L

  target_reads[stochastic_dropout] <- 0L

  data$shared_library_effect <- shared_library_effect
  data$target_library_effect <- target_library_effect
  data$target_log10_library <- target_log10_library
  data$total_reads <- total_reads
  data$mean_target_fraction <- mean_target_fraction
  data$realised_target_fraction <- realised_target_fraction
  data$sequencing_dropout_probability <- dropout_probability
  data$sequencing_dropout <- stochastic_dropout
  data$target_reads <- target_reads
  data$log10_target_reads <- log10(target_reads + 1)

  data
}
