# The parallel map and the model evaluations routed through it.

# DD99 solved by the Newton equilibrium solver rather than the model's own
# equilibrium, so that where each solve starts from matters to its answer,
# and started near its densities (K0 = 500): from the default 1e-3 the capped
# Newton steps stall.
dd99_newton <- function(r_jacobian = 1L) {
  community_start(bounds(x = c(-2, 2)), trait_scale = "linear", birth_rate_initial = 300,
                  demography_control = demographic_step_control(list(equilibrium_solver_name = "equilibrium_solve_newton")),
                  derivative_control = derivative_control(list(r_jacobian = r_jacobian)),
                  harness = harness_dd99(sigma_C = 0.4))
}
dd99_newton_off <- function(r_jacobian = 1L) {
  dd99_newton(r_jacobian) |>
    community_add(trait_matrix(c(-0.3, 0.6), "x"), birth_rate = c(300, 300)) |>
    community_demography()
}

# ---- the routed paths are pinned ----------------------------------------------
# Where each stencil solve starts from, and the tolerance it is solved to, move
# these values far beyond 1e-12 but within the oracle tolerances of the other
# tests, so they are pinned here: every point of a resident Jacobian starts
# from the centre's densities and solver state, which is also what makes the
# answer the same under any plan.

test_that("the resident Jacobian, singularity solve and classification are pinned", {
  J <- community_selection_gradient_jacobian(dd99_newton_off(r_jacobian = 2L))
  expect_equal(as.numeric(J), c(-2.3204972883327541, 2.1524888129590094, 1.6746455137339302, -2.6885672276952843),
               tolerance = 1e-12)
  expect_identical(attr(J, "evaluations"), 8L)

  pair <- community_solve_singularity(dd99_newton(), x0 = trait_matrix(c(-0.3, 0.6), "x"),
                                      birth_rate = c(300, 300))
  expect_true(attr(pair, "converged"))
  expect_equal(as.numeric(pair$traits), c(-0.44202687274337538, 0.44202688729386219), tolerance = 1e-12)
  expect_identical(attr(pair, "evaluations"), 10L)

  cl <- community_classify_singularity(pair)
  expect_equal(as.numeric(cl$jacobian),
               c(-2.5329623198889788, 1.9609547177981892, 1.9609546953174599, -2.5329623368624841),
               tolerance = 1e-12)
  expect_equal(as.numeric(cl$invasion_fitness), c(0.91304347863477864, 0.91304347751622095), tolerance = 1e-12)
  expect_identical(cl$evaluations, 6L)
  expect_equal(cl$classification, "branching point")
})

test_that("the parameter Jacobian and the canonical branching rate are pinned", {
  off <- dd99_newton_off()
  G <- community_selection_gradient_parameter_jacobian(off, community_parameter_map(off, c("x0", "sigma_K")))
  expect_equal(as.numeric(G), c(0.64585172568761351, 0.53607835036222762, -0.70624586393160793, 1.0608255475494446),
               tolerance = 1e-12)
  expect_identical(attr(G, "evaluations"), 4L)

  one <- dd99_newton() |>
    community_add(trait_matrix(0, "x"), birth_rate = 300) |>
    community_demography()
  br <- canonical_branch_rate(community_clear_residents(one), one, community_trait_transform(one),
                              canonical_control(), matrix(0, 1, 1), 1L, 1, 0.05, 1)
  expect_equal(br$rate, 0.01299719785847951, tolerance = 1e-12)
  expect_identical(br$evaluations, 12L)
  expect_identical(br$failed, 0L)
})

# ---- the map ------------------------------------------------------------------

test_that("residents are split into contiguous chunks, and without a plan the map is lapply", {
  expect_equal(regnans_chunks(0), list())
  expect_equal(regnans_chunks(5, parallel = FALSE), list(1:5))
  expect_equal(regnans_map(1:3, function(i) i * 2, parallel = FALSE), list(2, 4, 6))
  expect_identical(regnans_workers(parallel = FALSE), 1L)
  expect_identical(regnans_pointers(list(1, "a", function(x) x)), 0L)
  expect_identical(regnans_pointers(list(new("externalptr"), list(new("externalptr")))), 2L)
})

test_that("a plan's two workers are separate processes running this source", {
  for (type in c("multicore", "multisession")) {
    with_two_workers(type, {
      expect_identical(regnans_workers(), 2L)
      expect_equal(sort(lengths(regnans_chunks(11))), c(5L, 6L))
      expect_equal(unlist(regnans_chunks(11)), 1:11)
      w <- worker_identity()
      pid <- vapply(w, `[[`, integer(1), "pid")
      expect_length(unique(pid), 2L)
      expect_false(Sys.getpid() %in% pid)
      expect_equal(vapply(w, `[[`, character(1), "path"),
                   rep(getNamespaceInfo("regnans", "path"), 2))
    })
  }
})

test_that("a plan leaves the caller's random-number stream alone, and draws on workers are reproducible", {
  set.seed(1)
  before <- runif(1)
  draws <- list()
  for (type in c("multicore", "multisession")) {
    with_two_workers(type, {
      set.seed(1)
      draws[[type]] <- regnans_map(1:4, function(i) runif(1))
      expect_identical(runif(1), before)
      set.seed(1)
      expect_identical(regnans_map(1:4, function(i) runif(1)), draws[[type]])
    })
  }
  expect_identical(draws$multicore, draws$multisession)
})

test_that("work holding an external pointer is refused unless the workers are forks", {
  ptr <- new("externalptr")
  with_two_workers("multicore", {
    expect_equal(regnans_map(1:2, function(i, p) typeof(p), p = ptr), list("externalptr", "externalptr"))
  })
  with_two_workers("multisession", {
    expect_error(regnans_map(1:2, function(i, p) i, p = ptr), "external pointer.*plan\\(future::multicore\\)")
    # a closure carries its enclosing frame, and here that frame holds ptr:
    # what the map exists to catch
    expect_error(regnans_map(1:2, function(i) i), "external pointer")
    expect_equal(regnans_map(1:2, identity), list(1L, 2L))
  })
})

# ---- every routed path, sequentially and under each plan ------------------------
# No solve depends on the plan: the answer and its cost are the same bit for
# bit, with the model's own equilibrium and with an iterated one.
dd99_model <- function() {
  community_start(bounds(x = c(-2, 2)), trait_scale = "linear", harness = harness_dd99(sigma_C = 0.4))
}
routed <- function(start) {
  off <- start() |>
    community_add(trait_matrix(c(-0.3, 0.6), "x"), birth_rate = c(300, 300)) |>
    community_demography()
  pair <- community_solve_singularity(start(), x0 = trait_matrix(c(-0.3, 0.6), "x"), birth_rate = c(300, 300))
  one <- start() |> community_add(trait_matrix(0, "x"), birth_rate = 300) |> community_demography()
  list(J = community_selection_gradient_jacobian(off),
       pair = pair,
       classification = community_classify_singularity(pair),
       G = community_selection_gradient_parameter_jacobian(off, community_parameter_map(off, c("x0", "sigma_K"))),
       branch = canonical_branch_rate(community_clear_residents(one), one, community_trait_transform(one),
                                      canonical_control(), matrix(0, 1, 1), 1L, 1, 0.05, 1),
       canonical = community_canonical_equation(start(), x0 = 0.8,
         control = canonical_control(list(branch = "stochastic", max_residents = 2, mutation_sd = 0.1, seed = 7))))
}
expect_routed_identical <- function(a, b) {
  expect_identical(a$J, b$J)
  expect_identical(a$pair$traits, b$pair$traits)
  expect_identical(attr(a$pair, "evaluations"), attr(b$pair, "evaluations"))
  expect_identical(a$classification, b$classification)
  expect_identical(a$G, b$G)
  expect_identical(a$branch, b$branch)
  expect_identical(a$canonical$trajectory, b$canonical$trajectory)
  expect_identical(a$canonical$events, b$canonical$events)
  expect_identical(a$canonical$evaluations, b$canonical$evaluations)
}

test_that("every routed path gives the sequential answer and cost under each plan", {
  for (start in list(dd99_model, dd99_newton)) {
    sequential <- routed(start)
    expect_true(attr(sequential$pair, "converged"))
    expect_true(any(sequential$canonical$events$event == "branch"))
    for (type in c("multicore", "multisession")) {
      expect_routed_identical(with_two_workers(type, routed(start)), sequential)
    }
  }
})

test_that("a stencil point that merges two residents is refused under a plan as it is sequentially", {
  # 8e-4 apart at x ~ 0.5, inside two finite-difference steps (1e-3 |x| each)
  close <- dd99_model() |>
    community_add(trait_matrix(c(0.5, 0.5008), "x"), birth_rate = c(100, 100)) |>
    community_demography()
  expect_error(community_selection_gradient_jacobian(close),
               "residents 1 and 2 are within the finite-difference step")
  for (type in c("multicore", "multisession")) {
    with_two_workers(type, {
      expect_error(community_selection_gradient_jacobian(close),
                   "residents 1 and 2 are within the finite-difference step")
    })
  }
})

test_that("a stencil solve that fails still counts the solves that finished", {
  off <- dd99_model() |>
    community_add(trait_matrix(c(-0.3, 0.6), "x"), birth_rate = c(300, 300)) |>
    community_demography()
  gradient <- singularity_gradient_fn(off, 2L)
  attr(gradient, "prime")(off)
  solve <- singularity_solve_point
  local_mocked_bindings(singularity_solve_point = function(x, ...) {
    if (x[2] > 0.6) stop("the solve failed")
    solve(x, ...)
  })
  expect_error(attr(gradient, "points")(list(c(-0.3, 0.59), c(-0.3, 0.61), c(-0.31, 0.6)), 1e-3), "the solve failed")
  expect_identical(attr(gradient, "evaluations")(), 2L)
})
