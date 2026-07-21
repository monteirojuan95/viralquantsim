#' Run a two-parameter sensitivity grid
#'
#' Evaluates how the performance difference between a reads-only model and a
#' reads-plus-Delta-Ct model changes across combinations of two simulation
#' parameters. Negative `delta_rmse` values mean the internal-control-adjusted
#' model performed better.
#'
#' @param parameter_x Name of the first simulation parameter.
#' @param values_x Numeric values for the first parameter.
#' @param parameter_y Name of the second simulation parameter.
#' @param values_y Numeric values for the second parameter.
#' @param base_parameters A complete parameter list.
#' @param repetitions Number of independent simulations per grid cell.
#' @param cross_validate Whether to use cross-validated RMSE.
#'
#' @return A data frame with one row per parameter combination and repetition.
#' @examples
#' \donttest{
#' grid <- run_sensitivity_grid(
#'   "shared_recovery_sd", c(0.1, 0.4),
#'   "ic_specific_recovery_sd", c(0.05, 0.3),
#'   repetitions = 2L
#' )
#' }
#' @export
run_sensitivity_grid <- function(
    parameter_x,
    values_x,
    parameter_y,
    values_y,
    base_parameters = default_parameters(),
    repetitions = 10L,
    cross_validate = TRUE) {
  .validate_parameters(base_parameters)

  valid_names <- names(base_parameters)
  if (!parameter_x %in% valid_names || !parameter_y %in% valid_names) {
    stop("Both parameter names must exist in `base_parameters`.", call. = FALSE)
  }

  grid <- expand.grid(
    value_x = values_x,
    value_y = values_y,
    repetition = seq_len(as.integer(repetitions)),
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )

  rows <- lapply(seq_len(nrow(grid)), function(index) {
    scenario <- grid[index, , drop = FALSE]
    parameters <- base_parameters
    parameters[[parameter_x]] <- scenario$value_x
    parameters[[parameter_y]] <- scenario$value_y
    parameters$seed <- as.integer(base_parameters$seed + index - 1L)

    simulated <- simulate_assay(parameters)
    models <- fit_candidate_models(simulated)
    comparison <- compare_models(
      models,
      simulated,
      cross_validate = cross_validate,
      v = 5L,
      repeats = 1L,
      seed = parameters$seed
    )

    metric <- if (isTRUE(cross_validate)) "cv_rmse" else "rmse"
    reads_rmse <- comparison[comparison$model == "reads", metric]
    adjusted_rmse <- comparison[comparison$model == "reads_delta_ct", metric]

    data.frame(
      parameter_x = parameter_x,
      value_x = scenario$value_x,
      parameter_y = parameter_y,
      value_y = scenario$value_y,
      repetition = scenario$repetition,
      reads_rmse = reads_rmse,
      reads_delta_ct_rmse = adjusted_rmse,
      delta_rmse = adjusted_rmse - reads_rmse,
      stringsAsFactors = FALSE
    )
  })

  do.call(rbind, rows)
}
