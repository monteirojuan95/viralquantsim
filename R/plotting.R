#' Plot observed versus predicted viral load
#'
#' @param model An `lm` object.
#' @param title Optional plot title.
#'
#' @return A `ggplot` object.
#' @examples
#' \dontrun{
#' simulated <- simulate_assay(default_parameters())
#' models <- fit_candidate_models(simulated)
#' plot_model_predictions(models$reads_delta_ct)
#' }
#' @export
plot_model_predictions <- function(model, title = NULL) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package `ggplot2` is required for plotting.", call. = FALSE)
  }

  frame <- stats::model.frame(model)
  plot_data <- data.frame(
    observed = stats::model.response(frame),
    predicted = stats::fitted(model)
  )

  ggplot2::ggplot(
    plot_data,
    ggplot2::aes(x = observed, y = predicted)
  ) +
    ggplot2::geom_point(alpha = 0.65) +
    ggplot2::geom_abline(intercept = 0, slope = 1, linetype = 2) +
    ggplot2::coord_equal() +
    ggplot2::labs(
      title = title,
      x = "True viral load (log10 copies/mL)",
      y = "Predicted viral load (log10 copies/mL)"
    ) +
    ggplot2::theme_minimal()
}

#' Plot the benefit of internal-control adjustment
#'
#' @param sensitivity_data Data returned by [run_sensitivity_grid()].
#'
#' @return A `ggplot` heatmap. Negative values favour the adjusted model.
#' @examples
#' \dontrun{
#' grid <- run_sensitivity_grid(
#'   "shared_recovery_sd", c(0.1, 0.4),
#'   "ic_specific_recovery_sd", c(0.05, 0.3), repetitions = 2L
#' )
#' plot_internal_control_benefit(grid)
#' }
#' @export
plot_internal_control_benefit <- function(sensitivity_data) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package `ggplot2` is required for plotting.", call. = FALSE)
  }

  required <- c("value_x", "value_y", "delta_rmse")
  missing <- setdiff(required, names(sensitivity_data))
  if (length(missing) > 0L) {
    stop("Missing required sensitivity columns.", call. = FALSE)
  }

  aggregated <- stats::aggregate(
    delta_rmse ~ value_x + value_y,
    data = sensitivity_data,
    FUN = mean,
    na.rm = TRUE
  )

  ggplot2::ggplot(
    aggregated,
    ggplot2::aes(x = value_x, y = value_y, fill = delta_rmse)
  ) +
    ggplot2::geom_tile() +
    ggplot2::labs(
      x = unique(sensitivity_data$parameter_x),
      y = unique(sensitivity_data$parameter_y),
      fill = "Delta RMSE",
      subtitle = "Negative values favour reads + Delta Ct"
    ) +
    ggplot2::theme_minimal()
}
