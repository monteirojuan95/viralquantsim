#' Run the ViralQuantSim Shiny application
#'
#' @param ... Arguments passed to [shiny::runApp()].
#'
#' @return No return value; launches a Shiny application.
#' @examples
#' \dontrun{
#' run_app()
#' }
#' @export
run_app <- function(...) {
  required <- c("shiny", "bslib", "ggplot2", "plotly", "DT")
  missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]

  if (length(missing) > 0L) {
    stop(
      "Install the following packages before running the app: ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }

  app_directory <- system.file("shiny", package = "viralquantsim")
  shiny::runApp(app_directory, ...)
}
