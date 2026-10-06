##' @import plant
##' @importFrom Rcpp evalCpp
##' @importFrom stats optim optimise rpois uniroot
##' @importFrom utils capture.output modifyList
##' @importFrom dplyr mutate relocate any_of
##' @importFrom ggplot2 .data autoplot
##' @useDynLib regnans, .registration = TRUE
NULL

## Column names used in non-standard evaluation (dplyr/ggplot2), declared here
## so R CMD check does not flag them as undefined global variables.
utils::globalVariables(c("data", "resident", "strategy_id"))
