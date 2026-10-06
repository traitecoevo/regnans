## Cost of the canonical equation under each odelia stepper, in equilibrium
## solves (the only cost that matters), on the reference models.
## Run from the package root: Rscript scripts/canonical-stepper-benchmark.R
suppressMessages(devtools::load_all(quiet = TRUE))
eval(parse(text = readLines("tests/testthat/test-canonical.R")[1:28]))

cases <- list(
  jj12_css = list(comm = jj12(), x0 = -1.5, extra = list()),
  dd99_split = list(comm = dd99(sigma_C = 0.4, x0 = 0), x0 = 0.8, extra = list(branch = "immediate", max_residents = 2)),
  dd99_split_unpolished = list(comm = dd99(sigma_C = 0.4, x0 = 0), x0 = 0.8,
                               extra = list(branch = "immediate", max_residents = 2, polish = FALSE, max_steps = 400)),
  gk98_split = list(comm = gk98(d = 1.5), x0 = 1.2, extra = list(branch = "immediate", max_residents = 2))
)
steppers <- c("rodas", "rkck", "dopri")

rows <- list()
for (nm in names(cases)) {
  cs <- cases[[nm]]
  for (stepper in steppers) {
    ctrl <- canonical_control(c(list(stepper = stepper), cs$extra))
    el <- system.time(ce <- suppressWarnings(community_canonical_equation(cs$comm, x0 = cs$x0, control = ctrl)))[["elapsed"]]
    fin <- final(ce)
    rows[[length(rows) + 1L]] <- data.frame(
      case = nm, stepper = stepper,
      outcome = ce$outcome, steps = ce$steps, rejections = ce$rejections, solves = ce$evaluations,
      t_end = signif(max(ce$trajectory$time), 4), max_grad = signif(max(abs(fin[[grep("^gradient_", names(fin))[1]]])), 2),
      seconds = round(el, 1))
  }
}
out <- do.call(rbind, rows)
print(out, row.names = FALSE)
