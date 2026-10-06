# harness_numerical(): the reference models reaching their equilibrium by
# iterating their own dynamics, so that the package's solvers actually solve.

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
    numerical <- community_start(cs[[2]], trait_scale = cs[[3]], harness = harness_numerical(h)) |>
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
  h <- harness_numerical(harness_gm99(alpha = 7, beta = 15))
  base <- community_start(bounds(x = c(0.08, 0.9)), trait_scale = "log", harness = h)
  cold <- base |> community_add(trait_matrix(0.3, "x")) |> community_demography()
  warm <- base |> community_add(trait_matrix(0.3, "x"), birth_rate = 0.98 * cold$birth_rate) |>
    community_demography()
  expect_lt(NROW(attr(warm, "progress")), NROW(attr(cold, "progress")))
  expect_equal(as.numeric(warm$birth_rate), as.numeric(cold$birth_rate), tolerance = 1e-4)
})

test_that("harness_numerical applies to explicit harnesses only and prints its mode", {
  expect_error(harness_numerical(harness_plant()), "explicit")
  out <- paste(utils::capture.output(print(harness_numerical(harness_gk98()))), collapse = "\n")
  expect_match(out, "iterated from the model's own dynamics")
})

test_that("a dimorphic community reaches the same equilibrium either way", {
  h <- harness_gk98(d = 1.5)
  x <- trait_matrix(c(-0.8, 0.8), "x")
  analytic <- community_start(bounds(x = c(-3, 3)), trait_scale = "linear", harness = h) |>
    community_add(x, birth_rate = 1) |> community_demography()
  numerical <- community_start(bounds(x = c(-3, 3)), trait_scale = "linear",
                               harness = harness_numerical(h)) |>
    community_add(x, birth_rate = 1) |> community_demography()
  expect_equal(as.numeric(numerical$birth_rate), as.numeric(analytic$birth_rate), tolerance = 1e-4)
})
