# The "model" equilibrium solver and the iterated demography of the reference
# models. The reference harnesses supply their own equilibrium (the "model"
# solver, their default) and a one-generation demography runner that every
# other solver iterates, so the solvers can be tested against known answers.

cases <- list(
  list(harness_dd99(), bounds(x = c(-2, 2)), "linear", 0.4, 100),
  list(harness_gk98(d = 1.5), bounds(x = c(-3, 3)), "linear", 0.3, 1),
  list(harness_jj12(), bounds(x = c(-3, 3)), "linear", 0.2, 1),
  list(harness_gm99(alpha = 7, beta = 15), bounds(x = c(0.08, 0.9)), "log", 0.3, 1)
)

solve_with <- function(cs, solver = NULL, eps = 1e-8) {
  ctrl <- if (is.null(solver)) demographic_step_control() else
    demographic_step_control(list(equilibrium_solver_name = solver, equilibrium_eps = eps,
                                  equilibrium_nsteps = 1000))
  community_start(cs[[2]], trait_scale = cs[[3]], harness = cs[[1]], demography_control = ctrl) |>
    community_add(trait_matrix(cs[[4]], "x"), birth_rate = cs[[5]]) |>
    community_demography()
}

test_that("a control that leaves the solver open takes the harness default", {
  expect_null(demographic_step_control()$equilibrium_solver_name)
  expect_equal(community_start(bounds(x = c(-2, 2)), harness = harness_dd99())$demography_control$equilibrium_solver_name,
               "model")
  expect_equal(community_start(bounds(lma = c(0.05, 2)), harness = harness_plant())$demography_control$equilibrium_solver_name,
               "equilibrium_iteration")
  expect_equal(community_start(bounds(x = c(-2, 2)), harness = harness_dd99(),
                               demography_control = demographic_step_control(
                                 list(equilibrium_solver_name = "equilibrium_hybrid")))$demography_control$equilibrium_solver_name,
               "equilibrium_hybrid")
  expect_error(demographic_step_control(list(equilibrium_solver_name = "newton")), "should be one of")
  expect_error(community_start(bounds(lma = c(0.05, 2)), harness = harness_plant(),
                               demography_control = demographic_step_control(list(equilibrium_solver_name = "model"))),
               "needs a harness that supplies its own equilibrium")
  expect_true(harness_provides(harness_gk98(), "equilibrium"))
  expect_false(harness_provides(harness_plant(), "equilibrium"))
})

test_that("the model solver returns the model's own equilibrium in one evaluation", {
  for (cs in cases) {
    comm <- solve_with(cs)
    expect_true(attr(comm, "converged"), info = cs[[1]]$label)
    expect_equal(as.numeric(comm$birth_rate), as.numeric(cs[[1]]$equilibrium(cs[[4]])), info = cs[[1]]$label)
    expect_equal(NROW(attr(comm, "progress")), 1L, info = cs[[1]]$label)
    expect_equal(as.numeric(comm$resident_fitness), 0, tolerance = 1e-8, info = cs[[1]]$label)
  }
})

test_that("iterating the demography reaches the same equilibrium for every model", {
  for (cs in cases) {
    closed <- solve_with(cs)
    iterated <- solve_with(cs, "equilibrium_iteration")
    expect_true(attr(iterated, "converged"), info = cs[[1]]$label)
    expect_equal(as.numeric(iterated$birth_rate), as.numeric(closed$birth_rate),
                 tolerance = 1e-6, info = cs[[1]]$label)
    # it iterated: GK98's single-resident recursion lands in one step, the others take several
    expect_gte(NROW(attr(iterated, "progress")), 2L, label = cs[[1]]$label)
    if (cs[[1]]$label != "gk98") expect_gt(NROW(attr(iterated, "progress")), 4L, label = cs[[1]]$label)
  }
})

test_that("every equilibrium solver reaches the model's own equilibrium through the iterated demography", {
  cs <- cases[[4]]
  closed <- solve_with(cs)
  for (solver in c("equilibrium_iteration", "equilibrium_solve_newton", "equilibrium_solve_nleqslv",
                   "equilibrium_solve_dfsane", "equilibrium_hybrid")) {
    comm <- solve_with(cs, solver)
    expect_true(attr(comm, "converged"), info = solver)
    expect_equal(as.numeric(comm$birth_rate), as.numeric(closed$birth_rate),
                 tolerance = 1e-6, info = solver)
    expect_gt(NROW(attr(comm, "progress")), 1L, label = solver)
  }
})

test_that("near GM99's viability edge Newton needs far fewer evaluations than fixed-point iteration", {
  # at x = 0.9 the one-generation map's multiplier is ~0.9, so the iteration
  # crawls; Newton does not care
  cs <- list(harness_gm99(alpha = 7, beta = 15), bounds(x = c(0.08, 0.95)), "log", 0.9, 1)
  iterated <- solve_with(cs, "equilibrium_iteration")
  newton <- solve_with(cs, "equilibrium_solve_newton")
  expect_equal(as.numeric(newton$birth_rate), as.numeric(iterated$birth_rate), tolerance = 1e-5)
  expect_gt(NROW(attr(iterated, "progress")), 40L)
  expect_lt(NROW(attr(newton, "progress")), NROW(attr(iterated, "progress")) / 3)
})

test_that("a start near the equilibrium converges in fewer evaluations than a cold one", {
  cs <- cases[[4]]
  ctrl <- demographic_step_control(list(equilibrium_solver_name = "equilibrium_iteration"))
  base <- community_start(cs[[2]], trait_scale = cs[[3]], harness = cs[[1]], demography_control = ctrl)
  cold <- base |> community_add(trait_matrix(0.3, "x")) |> community_demography()
  warm <- base |> community_add(trait_matrix(0.3, "x"), birth_rate = 0.98 * cold$birth_rate) |>
    community_demography()
  expect_lt(NROW(attr(warm, "progress")), NROW(attr(cold, "progress")))
  expect_equal(as.numeric(warm$birth_rate), as.numeric(cold$birth_rate), tolerance = 1e-4)
})

test_that("a dimorphic community reaches the same equilibrium either way", {
  h <- harness_gk98(d = 1.5)
  x <- trait_matrix(c(-0.8, 0.8), "x")
  closed <- community_start(bounds(x = c(-3, 3)), trait_scale = "linear", harness = h) |>
    community_add(x, birth_rate = 1) |> community_demography()
  iterated <- community_start(bounds(x = c(-3, 3)), trait_scale = "linear", harness = h,
                              demography_control = demographic_step_control(
                                list(equilibrium_solver_name = "equilibrium_iteration", equilibrium_eps = 1e-8))) |>
    community_add(x, birth_rate = 1) |> community_demography()
  expect_equal(as.numeric(iterated$birth_rate), as.numeric(closed$birth_rate), tolerance = 1e-5)
})

test_that("print.harness says how the equilibrium is found", {
  out <- paste(utils::capture.output(print(harness_gk98())), collapse = "\n")
  expect_match(out, "supplied by the model")
  out <- paste(utils::capture.output(print(harness_plant())), collapse = "\n")
  expect_match(out, "iterating the demography runner")
})

test_that("Newton reaches GM99's equilibrium from the cold default start at every resident", {
  h <- harness_gm99(alpha = 7, beta = 15)
  ctrl <- demographic_step_control(list(equilibrium_solver_name = "equilibrium_solve_newton",
                                        equilibrium_eps = 1e-8, equilibrium_nsteps = 1000))
  for (x in exp(seq(log(0.08), log(0.9), length.out = 9))) {
    closed <- community_start(bounds(x = c(0.08, 0.9)), trait_scale = "log", harness = h) |>
      community_add(trait_matrix(x, "x")) |> community_demography()
    newton <- expect_no_warning(
      community_start(bounds(x = c(0.08, 0.9)), trait_scale = "log", harness = h, demography_control = ctrl) |>
        community_add(trait_matrix(x, "x")) |> community_demography())
    expect_true(attr(newton, "converged"), info = x)
    expect_equal(as.numeric(newton$birth_rate), as.numeric(closed$birth_rate), tolerance = 1e-6, info = x)
  }
})
