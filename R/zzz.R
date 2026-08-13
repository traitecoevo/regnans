##' @import plant
##' @importFrom Rcpp evalCpp
##' @importFrom stats rpois uniroot
##' @importFrom utils capture.output modifyList
##' @importFrom tibble tibble is_tibble
##' @importFrom purrr map_dbl
##' @importFrom dplyr select mutate relocate any_of
##' @importFrom ggplot2 ggplot aes aes_string geom_point geom_line geom_abline geom_text scale_x_log10 scale_y_log10 xlab ylab theme_classic theme element_text
##' @importFrom ggplot2 .data
##' @useDynLib regnans, .registration = TRUE
NULL

## Column names used in non-standard evaluation (dplyr/ggplot2), declared here
## so R CMD check does not flag them as undefined global variables.
utils::globalVariables(c(
  "batch_nr", "births", "data", "fitness", "invader",
  "resident", "strategy_id", "y"
))
