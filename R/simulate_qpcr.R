.row_mean_or_na <- function(x) {
  value <- rowMeans(x, na.rm = TRUE)
  value[is.nan(value)] <- NA_real_
  value
}

.row_sd_or_na <- function(x) {
  apply(x, 1L, function(values) {
    values <- values[!is.na(values)]
    if (length(values) < 2L) {
      return(NA_real_)
    }
    stats::sd(values)
  })
}

#' Simulate target and internal-control qPCR measurements
#'
#' Simulates technical qPCR replicates after shared and analyte-specific
#' recovery, reverse-transcription, and inhibition effects have been applied.
#' Delta Ct is defined as target Ct minus internal-control Ct.
#'
#' @param data A data frame containing `target_log10_available` and
#'   `ic_log10_available`.
#' @param parameters A complete parameter list, normally produced by
#'   [default_parameters()].
#'
#' @return The input data frame with qPCR summary columns and replicate
#'   matrices stored as list-columns.
#' @examples
#' input <- data.frame(
#'   target_log10_available = c(3, 4, 5),
#'   ic_log10_available = c(4, 4, 4)
#' )
#' qpcr <- simulate_qpcr(input)
#' qpcr[, c("target_ct", "internal_control_ct", "delta_ct")]
#' @export
simulate_qpcr <- function(data, parameters = default_parameters()) {
  .validate_parameters(parameters)

  required <- c("target_log10_available", "ic_log10_available")
  missing <- setdiff(required, names(data))
  if (length(missing) > 0L) {
    stop(
      "`data` is missing required columns: ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }

  n <- nrow(data)
  r <- as.integer(parameters$n_technical_replicates)

  shared_inhibition <- pmax(
    stats::rnorm(
      n,
      mean = parameters$shared_inhibition_mean_ct,
      sd = parameters$shared_inhibition_sd_ct
    ),
    0
  )

  target_inhibition <- shared_inhibition + stats::rnorm(
    n,
    mean = 0,
    sd = parameters$target_specific_inhibition_sd_ct
  )

  ic_inhibition <- shared_inhibition + stats::rnorm(
    n,
    mean = 0,
    sd = parameters$ic_specific_inhibition_sd_ct
  )

  target_ct_expected <-
    parameters$target_qpcr_intercept -
    parameters$target_qpcr_slope * data$target_log10_available +
    target_inhibition

  ic_ct_expected <-
    parameters$ic_qpcr_intercept -
    parameters$ic_qpcr_slope * data$ic_log10_available +
    ic_inhibition

  target_replicates <- matrix(
    stats::rnorm(
      n * r,
      mean = rep(target_ct_expected, times = r),
      sd = parameters$qpcr_measurement_sd_ct
    ),
    nrow = n,
    ncol = r
  )

  ic_replicates <- matrix(
    stats::rnorm(
      n * r,
      mean = rep(ic_ct_expected, times = r),
      sd = parameters$qpcr_measurement_sd_ct
    ),
    nrow = n,
    ncol = r
  )

  target_replicates[target_replicates > parameters$ct_limit] <- NA_real_
  ic_replicates[ic_replicates > parameters$ct_limit] <- NA_real_

  target_ct <- .row_mean_or_na(target_replicates)
  ic_ct <- .row_mean_or_na(ic_replicates)

  data$shared_inhibition_ct <- shared_inhibition
  data$target_ct_expected <- target_ct_expected
  data$internal_control_ct_expected <- ic_ct_expected
  data$target_ct <- target_ct
  data$internal_control_ct <- ic_ct
  data$delta_ct <- target_ct - ic_ct
  data$target_ct_sd <- .row_sd_or_na(target_replicates)
  data$internal_control_ct_sd <- .row_sd_or_na(ic_replicates)
  data$target_detected <- !is.na(target_ct)
  data$internal_control_detected <- !is.na(ic_ct)
  data$assay_valid <- data$internal_control_detected
  data$target_ct_replicates <- I(split(target_replicates, row(target_replicates)))
  data$internal_control_ct_replicates <- I(split(ic_replicates, row(ic_replicates)))

  data
}
