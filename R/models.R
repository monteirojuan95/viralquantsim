.model_formulas <- function() {
  list(
    target_ct = true_log10_viral_load ~ target_ct,
    delta_ct = true_log10_viral_load ~ delta_ct,
    reads = true_log10_viral_load ~ log10_target_reads,
    reads_target_ct = true_log10_viral_load ~ log10_target_reads + target_ct,
    reads_delta_ct = true_log10_viral_load ~ log10_target_reads + delta_ct
  )
}

#' Fit candidate viral-load models
#'
#' Fits five interpretable linear models reflecting the main comparisons from
#' the original analytical rationale: target Ct, Delta Ct, reads, reads plus
#' target Ct, and reads plus Delta Ct.
#'
#' @param data A data frame returned by [simulate_assay()] or a compatible
#'   user-supplied dataset.
#'
#' @return A named list of `lm` objects.
#' @examples
#' simulated <- simulate_assay(default_parameters())
#' models <- fit_candidate_models(simulated)
#' names(models)
#' @export
fit_candidate_models <- function(data) {
  formulas <- .model_formulas()

  models <- lapply(formulas, function(formula) {
    stats::lm(formula, data = data, na.action = stats::na.omit)
  })

  class(models) <- c("viralquantsim_models", class(models))
  models
}

#' Calculate regression performance metrics
#'
#' @param observed Numeric observed values.
#' @param predicted Numeric predicted values.
#'
#' @return A one-row data frame.
#' @examples
#' model_metrics(observed = c(2, 3, 4), predicted = c(2.1, 2.9, 4.2))
#' @export
model_metrics <- function(observed, predicted) {
  keep <- stats::complete.cases(observed, predicted)
  observed <- observed[keep]
  predicted <- predicted[keep]

  if (length(observed) == 0L) {
    return(data.frame(
      n = 0L,
      rmse = NA_real_,
      mae = NA_real_,
      bias = NA_real_,
      within_0_25_log10 = NA_real_,
      within_0_50_log10 = NA_real_
    ))
  }

  error <- predicted - observed

  data.frame(
    n = length(observed),
    rmse = sqrt(mean(error^2)),
    mae = mean(abs(error)),
    bias = mean(error),
    within_0_25_log10 = mean(abs(error) <= 0.25),
    within_0_50_log10 = mean(abs(error) <= 0.50)
  )
}

.cross_validate_formula <- function(formula, data, v = 5L, repeats = 3L, seed = 42L) {
  model_data <- stats::model.frame(
    formula,
    data = data,
    na.action = stats::na.omit
  )

  n <- nrow(model_data)
  if (n < 3L) {
    return(data.frame(
      cv_n = n,
      cv_rmse = NA_real_,
      cv_mae = NA_real_,
      cv_bias = NA_real_,
      cv_within_0_50_log10 = NA_real_
    ))
  }

  v <- max(2L, min(as.integer(v), n))
  repeats <- max(1L, as.integer(repeats))
  set.seed(as.integer(seed))

  all_observed <- numeric(0)
  all_predicted <- numeric(0)

  for (repeat_id in seq_len(repeats)) {
    fold_id <- sample(rep(seq_len(v), length.out = n))

    for (fold in seq_len(v)) {
      train <- model_data[fold_id != fold, , drop = FALSE]
      test <- model_data[fold_id == fold, , drop = FALSE]

      if (nrow(train) < 2L || nrow(test) == 0L) {
        next
      }

      fit <- stats::lm(formula, data = train)
      prediction <- stats::predict(fit, newdata = test)
      observed <- stats::model.response(stats::model.frame(formula, data = test))

      all_observed <- c(all_observed, observed)
      all_predicted <- c(all_predicted, prediction)
    }
  }

  metrics <- model_metrics(all_observed, all_predicted)
  data.frame(
    cv_n = metrics$n,
    cv_rmse = metrics$rmse,
    cv_mae = metrics$mae,
    cv_bias = metrics$bias,
    cv_within_0_50_log10 = metrics$within_0_50_log10
  )
}

#' Compare candidate models
#'
#' @param models A named list returned by [fit_candidate_models()].
#' @param data The dataset used to fit the models.
#' @param cross_validate Whether to calculate repeated V-fold cross-validation.
#' @param v Number of folds.
#' @param repeats Number of repetitions.
#' @param seed Random seed for resampling.
#'
#' @return A data frame with one row per model.
#' @examples
#' simulated <- simulate_assay(default_parameters())
#' models <- fit_candidate_models(simulated)
#' compare_models(models, simulated, cross_validate = FALSE)
#' @export
compare_models <- function(
    models,
    data,
    cross_validate = TRUE,
    v = 5L,
    repeats = 3L,
    seed = 42L) {
  if (is.null(names(models)) || any(names(models) == "")) {
    stop("`models` must be a named list.", call. = FALSE)
  }

  rows <- lapply(names(models), function(model_name) {
    model <- models[[model_name]]
    frame <- stats::model.frame(model)
    observed <- stats::model.response(frame)
    predicted <- stats::fitted(model)
    apparent <- model_metrics(observed, predicted)

    row <- data.frame(
      model = model_name,
      formula = paste(deparse(stats::formula(model)), collapse = " "),
      n = apparent$n,
      r_squared = summary(model)$r.squared,
      adjusted_r_squared = summary(model)$adj.r.squared,
      rmse = apparent$rmse,
      mae = apparent$mae,
      bias = apparent$bias,
      within_0_25_log10 = apparent$within_0_25_log10,
      within_0_50_log10 = apparent$within_0_50_log10,
      aic = stats::AIC(model),
      bic = stats::BIC(model),
      stringsAsFactors = FALSE
    )

    if (isTRUE(cross_validate)) {
      cv <- .cross_validate_formula(
        stats::formula(model),
        data = data,
        v = v,
        repeats = repeats,
        seed = seed
      )
      row <- cbind(row, cv)
    }

    row
  })

  result <- do.call(rbind, rows)
  rownames(result) <- NULL
  ordering_metric <- if ("cv_rmse" %in% names(result)) result$cv_rmse else result$rmse
  result[order(ordering_metric, na.last = TRUE), , drop = FALSE]
}
