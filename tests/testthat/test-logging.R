test_that("a log message with braces is logged as written", {
  # The equilibrium runner prints a community's birth rates as "{1.66, 1.41}";
  # read as a glue expression that is a syntax error.
  msg <- "eq> 1: {1.661716, 1.411285} -> {1.64714, 0.004249168}"
  out <- capture.output(plant_log_info(msg, routine = "equilibrium"), type = "message")
  expect_true(any(grepl(msg, out, fixed = TRUE)))
})
