## The message is logged as written: braces in it (a vector printed as
## "{1.66, 1.41}") are text, not a glue expression.
plant_log_info <- function(message, routine=NULL, ...) {
  if (!is.null(routine)) message <- paste0(crayon::yellow(routine), "> ", message)
  logger::log_info(logger::skip_formatter(as.character(message)), namespace="plant", .topenv=parent.frame())
}

plant_log_debug <- function(message, routine=NULL, ...) {
  if (!is.null(routine)) message <- paste0(crayon::yellow(routine), "> ", message)
  logger::log_debug(logger::skip_formatter(as.character(message)), namespace="plant", .topenv=parent.frame())
}

plant_log_max_fitness <- function(...) {
  plant_log_info(..., routine="max_fitness")
}

plant_log_assembler <- function(...) {
  plant_log_info(..., routine="assembler")
}

plant_log_births <- function(...) {
  plant_log_info(..., routine="births")
}

plant_log_state <- function(...) {
  plant_log_info(..., routine = "state")
}

plant_log_deaths <- function(...) {
  plant_log_info(..., routine="deaths")
}
