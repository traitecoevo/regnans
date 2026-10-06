# harness_iterate_demography(): the reference models reaching their equilibrium by
# iterating their own demography, so that the package's solvers actually solve.

test_that("the iterated dynamics reach the analytic equilibrium of every model", {
  cases <- list(
    list(harness_dd99(), bounds(x = c(-2, 2)), "linear", 0.4, 100),
    list(harness_gk98(d = 1.5), bounds(x = c(-3, 3)), "linear", 0.3, 1),
    list(harness_jj12(), bounds(x = c(-3, 3)), "linear", 0.2, 1),
    list(harness_gm99(alpha = 7, beta = 15), bounds(x = c(0.08, 0.9)), "log", 0.3, 1)
  )
  for (cs in cases) {
    h <- cs[[1]]
    analytic <- community_start(cs[[2]], trait_scale = cs[[3]], harness = h) |>
      community_add(trait_matrix(cs[[4]], "x"), birth_rate = cs[[5]]) |>
      community_demography()
    numerical <- community_start(cs[[2]], trait_scale = cs[[3]], harness = harness_iterate_demography(h)) |>
      community_add(trait_matrix(cs[[4]], "x"), birth_rate = cs[[5]]) |>
      community_demography()
    expect_true(attr(numerical, "converged"), info = h$label)
    expect_equal(as.numeric(numerical$birth_rate), as.numeric(analytic$birth_rate),
                 tolerance = 1e-4, info = h$label)
    expect_equal(as.numeric(numerical$resident_fitness), 0, tolerance = 1e-4, info = h$label)
    # it iterated: GK98's single-resident recursion lands in one step, the
    # others take several
    expect_gte(NROW(attr(numerical, "progress")), 2L, label = h$label)
    if (h$label != "gk98") expect_gt(NROW(attr(numerical, "progress")), 4L, label = h$label)
  }
})

test_that("a start near the equilibrium converges in fewer evaluations than a cold one", {
  h <- harness_iterate_demography(harness_gm99(alpha = 7, beta = 15))
  base <- community_start(bounds(x = c(0.08, 0.9)), trait_scale = "log", harness = h)
  cold <- base |> community_add(trait_matrix(0.3, "x")) |> community_demography()
  warm <- base |> community_add(trait_matrix(0.3, "x"), birth_rate = 0.98 * cold$birth_rate) |>
    community_demography()
  expect_lt(NROW(attr(warm, "progress")), NROW(attr(cold, "progress")))
  expect_equal(as.numeric(warm$birth_rate), as.numeric(cold$birth_rate), tolerance = 1e-4)
})

test_that("harness_iterate_demography applies to explicit harnesses only and prints its mode", {
  expect_error(harness_iterate_demography(harness_plant()), "explicit")
  out <- paste(utils::capture.output(print(harness_iterate_demography(harness_gk98()))), collapse = "\n")
  expect_match(out, "iterated from the model's own demography")
})

test_that("a dimorphic community reaches the same equilibrium either way", {
  h <- harness_gk98(d = 1.5)
  x <- trait_matrix(c(-0.8, 0.8), "x")
  analytic <- community_start(bounds(x = c(-3, 3)), trait_scale = "linear", harness = h) |>
    community_add(x, birth_rate = 1) |> community_demography()
  numerical <- community_start(bounds(x = c(-3, 3)), trait_scale = "linear",
                               harness = harness_iterate_demography(h)) |>
    community_add(x, birth_rate = 1) |> community_demography()
  expect_equal(as.numeric(numerical$birth_rate), as.numeric(analytic$birth_rate), tolerance = 1e-4)
})

test_that("every equilibrium solver reaches the closed-form equilibrium through the iterated demography", {
  h <- harness_gm99(alpha = 7, beta = 15)
  closed <- community_start(bounds(x = c(0.08, 0.9)), trait_scale = "log", harness = h) |>
    community_add(trait_matrix(0.3, "x"), birth_rate = 1) |>
    community_demography()
  for (solver in c("equilibrium_iteration", "equilibrium_solve_nleqslv",
                   "equilibrium_solve_dfsane", "equilibrium_hybrid")) {
    comm <- community_start(bounds(x = c(0.08, 0.9)), trait_scale = "log",
                            harness = harness_iterate_demography(h),
                            demography_control = demographic_step_control(
                              list(equilibrium_solver_name = solver, equilibrium_eps = 1e-8))) |>
      community_add(trait_matrix(0.3, "x"), birth_rate = 1) |>
      community_demography()
    expect_true(attr(comm, "converged"), info = solver)
    expect_equal(as.numeric(comm$birth_rate), as.numeric(closed$birth_rate),
                 tolerance = 1e-6, info = solver)
    expect_gt(NROW(attr(comm, "progress")), 1L, label = solver)
  }
})
