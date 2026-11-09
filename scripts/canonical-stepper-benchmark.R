## Cost of the canonical equation under each odelia stepper, on the reference
## models, with the model's own equilibrium and with one iterated on the
## one-generation map (as plant's is). RODAS is run with its Jacobian both ways
## (derivative_control(resident_jacobian = )): differenced across equilibrium
## solves, and by implicit differentiation at a fixed environment. `runs`
## counts what a model run costs on plant -- one call of the one-generation map
## inside an equilibrium solve, or one evaluation at given traits and
## densities -- so the two routes can be compared in one unit.
## Run from the package root: Rscript scripts/canonical-stepper-benchmark.R
suppressMessages(devtools::load_all(quiet = TRUE))

runner_calls <- 0L
trace("explicit_community_demography_runner_cleanup", where = asNamespace("regnans"), print = FALSE,
      exit = quote({
        e <- environment(runner)
        if (is.function(e$runner_full)) e <- environment(e$runner_full)
        runner_calls <<- runner_calls + length(e$history)
      }))
eval(parse(text = readLines("tests/testthat/test-canonical.R")[1:28]))

cases <- list(
  jj12_css = list(comm = jj12(), x0 = -1.5, extra = list()),
  dd99_split = list(comm = dd99(sigma_C = 0.4, x0 = 0), x0 = 0.8, extra = list(branch = "immediate", max_residents = 2)),
  dd99_split_unpolished = list(comm = dd99(sigma_C = 0.4, x0 = 0), x0 = 0.8,
                               extra = list(branch = "immediate", max_residents = 2, polish = FALSE, max_steps = 400)),
  gk98_split = list(comm = gk98(d = 1.5), x0 = 1.2, extra = list(branch = "immediate", max_residents = 2))
)
routes <- list(c("rodas", "stencil"), c("rodas", "implicit"), c("rkck", "stencil"), c("dopri", "stencil"))
equilibria <- c("model", "equilibrium_iteration")

rows <- list()
for (equilibrium in equilibria) {
  for (nm in names(cases)) {
    cs <- cases[[nm]]
    comm <- cs$comm
    comm$demography_control$equilibrium_solver_name <- equilibrium
    for (route in routes) {
      comm$derivative_control <- derivative_control(list(resident_jacobian = route[2]))
      ctrl <- canonical_control(c(list(stepper = route[1]), cs$extra))
      runner_calls <- 0L
      el <- system.time(ce <- suppressWarnings(community_canonical_equation(comm, x0 = cs$x0, control = ctrl)))[["elapsed"]]
      fin <- final(ce)
      rows[[length(rows) + 1L]] <- data.frame(
        equilibrium = sub("equilibrium_", "", equilibrium), case = nm, stepper = route[1],
        jacobian = if (route[1] == "rodas") route[2] else "-",
        outcome = ce$outcome, steps = ce$steps, rejections = ce$rejections, solves = ce$evaluations,
        environments = ce$environments,
        runs = if (equilibrium == "model") NA_integer_ else runner_calls + ce$environments,
        t_end = signif(max(ce$trajectory$time), 4), max_grad = signif(max(abs(fin[[grep("^gradient_", names(fin))[1]]])), 2),
        seconds = round(el, 1))
    }
  }
}
out <- do.call(rbind, rows)
print(out, row.names = FALSE)
