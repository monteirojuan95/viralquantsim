#' Default simulation parameters
#'
#' Returns a named list containing the parameters used by [simulate_assay()].
#' The defaults describe a generic RNA-pathogen workflow and are not calibrated
#' to any specific clinical assay or private dataset.
#'
#' @return A named list.
#' @examples
#' params <- default_parameters()
#' params$n_samples <- 200
#' @export
default_parameters <- function() {
  list(
    n_samples = 500L,
    seed = 42L,
    viral_load_min = 1.5,
    viral_load_max = 6.0,

    ic_log10_input = 4.0,
    ic_pipetting_sd = 0.08,

    shared_recovery_sd = 0.30,
    target_specific_recovery_sd = 0.15,
    ic_specific_recovery_sd = 0.15,

    shared_rt_sd = 0.15,
    target_specific_rt_sd = 0.08,
    ic_specific_rt_sd = 0.08,

    shared_inhibition_mean_ct = 0.15,
    shared_inhibition_sd_ct = 0.12,
    target_specific_inhibition_sd_ct = 0.08,
    ic_specific_inhibition_sd_ct = 0.08,

    target_qpcr_intercept = 40.0,
    target_qpcr_slope = 3.32,
    ic_qpcr_intercept = 40.0,
    ic_qpcr_slope = 3.32,
    qpcr_measurement_sd_ct = 0.25,
    n_technical_replicates = 2L,
    ct_limit = 40.0,

    library_shared_sd = 0.10,
    target_library_specific_sd = 0.10,
    mean_total_reads = 2e6,
    total_reads_sdlog = 0.35,
    target_fraction_intercept = -13.0,
    target_fraction_slope = 2.0,
    sequencing_concentration = 400.0,
    dropout_midpoint_log10 = 2.0,
    dropout_steepness = 4.0
  )
}

.validate_parameters <- function(parameters) {
  required <- names(default_parameters())
  missing <- setdiff(required, names(parameters))

  if (length(missing) > 0L) {
    stop(
      "Missing simulation parameters: ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }

  if (parameters$n_samples < 3L) {
    stop("`n_samples` must be at least 3.", call. = FALSE)
  }

  if (parameters$viral_load_min >= parameters$viral_load_max) {
    stop("`viral_load_min` must be lower than `viral_load_max`.", call. = FALSE)
  }

  if (parameters$n_technical_replicates < 1L) {
    stop("`n_technical_replicates` must be at least 1.", call. = FALSE)
  }

  non_negative <- c(
    "ic_pipetting_sd",
    "shared_recovery_sd",
    "target_specific_recovery_sd",
    "ic_specific_recovery_sd",
    "shared_rt_sd",
    "target_specific_rt_sd",
    "ic_specific_rt_sd",
    "shared_inhibition_mean_ct",
    "shared_inhibition_sd_ct",
    "target_specific_inhibition_sd_ct",
    "ic_specific_inhibition_sd_ct",
    "qpcr_measurement_sd_ct",
    "library_shared_sd",
    "target_library_specific_sd",
    "total_reads_sdlog",
    "sequencing_concentration",
    "dropout_steepness"
  )

  invalid <- non_negative[vapply(
    parameters[non_negative],
    function(x) !is.numeric(x) || length(x) != 1L || is.na(x) || x < 0,
    logical(1)
  )]

  if (length(invalid) > 0L) {
    stop(
      "These parameters must be non-negative numeric scalars: ",
      paste(invalid, collapse = ", "),
      call. = FALSE
    )
  }

  if (parameters$mean_total_reads <= 0) {
    stop("`mean_total_reads` must be greater than zero.", call. = FALSE)
  }

  if (parameters$sequencing_concentration <= 0) {
    stop("`sequencing_concentration` must be greater than zero.", call. = FALSE)
  }

  invisible(parameters)
}
