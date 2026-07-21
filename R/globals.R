# Column names referenced inside ggplot2 aes() calls are resolved from the
# plotting data frame at run time, so R CMD check's static analysis cannot see
# them. Declaring them here removes the spurious "no visible binding" NOTE.
utils::globalVariables(c(
  "value_x",
  "value_y",
  "delta_rmse",
  "observed",
  "predicted"
))
