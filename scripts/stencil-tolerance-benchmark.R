# How accurate is a finite difference across equilibrium solves, and what
# does it cost?
#
# The resident Jacobian and the parameter Jacobian difference the selection
# gradient across equilibrium solves. On the reference models each is compared
# with its closed form (tests/testthat/helper-coalition.R), with the
# equilibrium iterated rather than taken from the model, at the default
# equilibrium_eps; on plant, where there is no closed form, the cost in
# runner calls and seconds is reported. Runner calls are counted for the
# derivative alone, after the community has been solved.
#
#   Rscript scripts/stencil-tolerance-benchmark.R

suppressMessages(devtools::load_all(".", quiet = TRUE))
source("tests/testthat/helper-coalition.R")
source("tests/testthat/helper-assembly.R")

## Count the calls of every demography runner, through the connector (a
## parameter map rebuilds the harness, so a count on the harness would miss
## them); the inner runner is named runner_full so that the cleanup connector
## reads the model's state from it.
calls <- new.env()
calls$n <- 0L
make_runner <- community_make_demography_runner
assignInNamespace("community_make_demography_runner", function(community) {
  runner_full <- make_runner(community)
  function(birth_rates) {
    calls$n <- calls$n + 1L
    runner_full(birth_rates)
  }
}, "regnans")
measure <- function(f) {
  calls$n <- 0L
  seconds <- system.time(value <- f())[["elapsed"]]
  list(value = value, runner_calls = calls$n, seconds = seconds)
}
max_rel <- function(a, b) max(abs(as.numeric(a) - as.numeric(b))) / max(abs(as.numeric(b)))

rows <- list()
add <- function(case, solver, r, m, error = NA_real_) {
  rows[[length(rows) + 1L]] <<- data.frame(
    case = case, solver = sub("equilibrium_(solve_)?", "", solver), r = r,
    error = signif(error, 2), runner_calls = m$runner_calls, seconds = round(m$seconds, 2))
}

## ---- reference models, against the closed form ---------------------------------
reference <- function(harness, x, solver, r, bounds, birth_rate) {
  community_start(bounds, trait_scale = "linear", birth_rate_initial = birth_rate,
                  demography_control = demographic_step_control(list(equilibrium_solver_name = solver)),
                  derivative_control = derivative_control(list(r_jacobian = r, r_parameter = r)),
                  harness = harness) |>
    community_add(trait_matrix(x, "x"), birth_rate = rep(birth_rate, length(x))) |>
    community_demography()
}
solvers <- c("equilibrium_solve_newton", "equilibrium_iteration")

x_dd <- c(-0.3, 0.6)
J_dd <- oracle_jacobian(dd99_coalition_gradient, x_dd)
G_dd <- oracle_jacobian(function(p) dd99_coalition_gradient(x_dd, x0 = p[1], sigma_K = p[2]), c(0, 1))
x_gm <- c(0.25, 0.65)
J_gm <- oracle_jacobian(function(x) gm99_coalition_gradient(x, alpha = 7, beta = 15), x_gm, h = 1e-5)

for (solver in solvers) {
  for (r in 1:2) {
    dd <- reference(harness_dd99(sigma_C = 0.4), x_dd, solver, r, bounds(x = c(-2, 2)), 300)
    m <- measure(function() community_selection_gradient_jacobian(dd))
    add("DD99 pair, resident J", solver, r, m, max_rel(m$value, J_dd))
    m <- measure(function()
      community_selection_gradient_parameter_jacobian(dd, community_parameter_map(dd, c("x0", "sigma_K"))))
    add("DD99 pair, parameter G", solver, r, m, max_rel(m$value, G_dd))
    gm <- reference(harness_gm99(alpha = 7, beta = 15), x_gm, solver, r, bounds(x = c(0.1, 0.95)), 1)
    m <- measure(function() community_selection_gradient_jacobian(gm))
    add("GM99 pair, resident J", solver, r, m, max_rel(m$value, J_gm))
  }
}

## ---- plant, cost only ----------------------------------------------------------
## the steps of test-plant-smoke-singularity.R (a 1e-3 relative step sits in
## the noise of the SCM's cohort schedule), and the default for one resident.
## The pair is not a coalition: its second resident is being excluded, at a
## density of ~1e-5 that halves every generation.
plant_start <- function(d = 1e-2) {
  community_start(bounds(lma = c(0.02, 0.6)), model_support = assembly_model_support(),
                  derivative_control = list(d_second = d, eps_second = d))
}
one <- plant_start() |>
  community_add(trait_matrix(0.0825, "lma"), birth_rate = 200) |>
  community_demography()
pair <- plant_start() |>
  community_add(trait_matrix(c(0.06, 0.3), "lma"), birth_rate = c(200, 200)) |>
  community_demography()
one_fine <- plant_start(1e-3) |>
  community_add(trait_matrix(0.0825, "lma"), birth_rate = 200) |>
  community_demography()

for (case in list(list("plant one resident, d = 1e-2", one), list("plant one resident, d = 1e-3", one_fine),
                  list("plant pair, one being excluded, d = 1e-2", pair))) {
  m <- measure(function() community_selection_gradient_jacobian(case[[2]]))
  add(case[[1]], "equilibrium_iteration", 1L, m)
  cat(case[[1]], ": J =", format(as.numeric(m$value), digits = 6), "\n")
}

print(do.call(rbind, rows), row.names = FALSE)
