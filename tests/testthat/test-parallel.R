# The parallel map and the model evaluations routed through it.

# DD99 solved by the Newton equilibrium solver rather than the model's own
# equilibrium, so that where each solve starts from matters to its answer.
dd99_newton <- function(r_jacobian = 1L, equilibrium_eps = 1e-5) {
  community_start(bounds(x = c(-2, 2)), trait_scale = "linear",
                  demography_control = demographic_step_control(list(equilibrium_solver_name = "equilibrium_solve_newton",
                                                                     equilibrium_eps = equilibrium_eps)),
                  derivative_control = derivative_control(list(r_jacobian = r_jacobian)),
                  harness = harness_dd99(sigma_C = 0.4))
}
dd99_newton_off <- function(r_jacobian = 1L, equilibrium_eps = 1e-5) {
  dd99_newton(r_jacobian, equilibrium_eps) |>
    community_add(trait_matrix(c(-0.3, 0.6), "x"), birth_rate = c(300, 300)) |>
    community_demography()
}

# ---- the sequential path is pinned -------------------------------------------
# With no plan, every routed path must be the arithmetic it was before routing,
# the warm-start chain from one stencil point to the next included. Seeding
# each point from the centre instead moves these values far beyond 1e-12 but
# within the oracle tolerances of the other tests, so they are pinned here.

test_that("the resident Jacobian, singularity solve and classification are pinned sequentially", {
  J <- community_selection_gradient_jacobian(dd99_newton_off(r_jacobian = 2L))
  expect_equal(as.numeric(J), c(-2.3204965543239702, 2.1524899760059748, 1.67464366695251, -2.6885701560685891),
               tolerance = 1e-12)
  expect_identical(attr(J, "evaluations"), 8L)

  pair <- community_solve_singularity(dd99_newton(), x0 = trait_matrix(c(-0.3, 0.6), "x"),
                                      birth_rate = c(300, 300))
  expect_true(attr(pair, "converged"))
  expect_equal(as.numeric(pair$traits), c(-0.4420268727422344, 0.44202688729495382), tolerance = 1e-12)
  expect_identical(attr(pair, "evaluations"), 10L)

  cl <- community_classify_singularity(pair)
  expect_equal(as.numeric(cl$jacobian),
               c(-2.5329623352992989, 1.9609547023982508, 1.9609546877876345, -2.5329623443894227),
               tolerance = 1e-12)
  expect_equal(as.numeric(cl$invasion_fitness), c(0.9130434786348407, 0.91304347751611148), tolerance = 1e-12)
  expect_identical(cl$evaluations, 6L)
  expect_equal(cl$classification, "branching point")
})

test_that("the parameter Jacobian and the canonical branching rate are pinned sequentially", {
  off <- dd99_newton_off()
  G <- community_selection_gradient_parameter_jacobian(off, community_parameter_map(off, c("x0", "sigma_K")))
  expect_equal(as.numeric(G), c(0.64585173892039993, 0.53607837136718239, -0.70624586226180475, 1.0608255501970212),
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

test_that("a plan leaves the caller's random-number stream alone", {
  set.seed(1)
  before <- runif(1)
  for (type in c("multicore", "multisession")) {
    with_two_workers(type, {
      set.seed(1)
      regnans_map(1:4, function(i) i)
      expect_identical(runif(1), before)
    })
  }
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

# With the model's own equilibrium no solve depends on where it starts, so
# solving the points apart gives the sequential answer bit for bit.
dd99_model <- function() {
  community_start(bounds(x = c(-2, 2)), trait_scale = "linear", harness = harness_dd99(sigma_C = 0.4))
}
routed_model <- function() {
  off <- dd99_model() |>
    community_add(trait_matrix(c(-0.3, 0.6), "x"), birth_rate = c(300, 300)) |>
    community_demography()
  pair <- community_solve_singularity(dd99_model(), x0 = trait_matrix(c(-0.3, 0.6), "x"))
  one <- dd99_model() |> community_add(trait_matrix(0, "x")) |> community_demography()
  list(J = community_selection_gradient_jacobian(off),
       pair = pair,
       classification = community_classify_singularity(pair),
       G = community_selection_gradient_parameter_jacobian(off, community_parameter_map(off, c("x0", "sigma_K"))),
       branch = canonical_branch_rate(community_clear_residents(one), one, community_trait_transform(one),
                                      canonical_control(), matrix(0, 1, 1), 1L, 1, 0.05, 1),
       canonical = community_canonical_equation(dd99_model(), x0 = 0.8,
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

test_that("every routed path gives the sequential answer and cost under each plan, with the model's equilibrium", {
  sequential <- routed_model()
  expect_true(any(sequential$canonical$events$event == "branch"))
  for (type in c("multicore", "multisession")) {
    expect_routed_identical(with_two_workers(type, routed_model()), sequential)
  }
})

# Where the equilibrium is iterated, a stencil point solved apart starts from
# the centre's densities and solver state rather than the previous point's.
# Either way a finite difference of a solve carries the solve's error divided
# by the step, so the two agree only as well as the equilibrium is solved: at
# a tight equilibrium tolerance both reach the closed-form Jacobian to the
# stencil's truncation error, and each other well inside it. The cost is the
# same.
routed_newton <- function() {
  pair <- community_solve_singularity(dd99_newton(equilibrium_eps = 1e-9), x0 = trait_matrix(c(-0.3, 0.6), "x"),
                                      birth_rate = c(300, 300))
  list(J = community_selection_gradient_jacobian(dd99_newton_off(r_jacobian = 2L, equilibrium_eps = 1e-9)),
       pair = pair, classification = community_classify_singularity(pair))
}

test_that("with an iterated equilibrium a plan gives the same derivatives to the solve's accuracy, at the same cost", {
  sequential <- routed_newton()
  oracle <- oracle_jacobian(dd99_coalition_gradient, c(-0.3, 0.6))
  root <- c(-1, 1) * dd99_pair_root()
  expect_equal(as.numeric(sequential$J), as.numeric(oracle), tolerance = 1e-5)
  for (type in c("multicore", "multisession")) {
    parallel <- with_two_workers(type, routed_newton())
    expect_equal(as.numeric(parallel$J), as.numeric(oracle), tolerance = 1e-5)
    expect_equal(parallel$J, sequential$J, tolerance = 1e-6)
    expect_identical(attr(parallel$J, "evaluations"), attr(sequential$J, "evaluations"))
    expect_true(attr(parallel$pair, "converged"))
    expect_equal(as.numeric(parallel$pair$traits), root, tolerance = 1e-7)
    expect_identical(attr(parallel$pair, "evaluations"), attr(sequential$pair, "evaluations"))
    expect_equal(parallel$classification$jacobian, sequential$classification$jacobian, tolerance = 1e-6)
    expect_equal(parallel$classification$invasion_fitness, sequential$classification$invasion_fitness,
                 tolerance = 1e-8)
    expect_identical(parallel$classification$evaluations, sequential$classification$evaluations)
    expect_equal(parallel$classification$classification, "branching point")
  }
})

test_that("derivative_control(parallel = FALSE) keeps the sequential numbers under a plan", {
  expect_error(derivative_control(list(parallel = NA)), "TRUE or FALSE")
  jacobian <- function(parallel) {
    comm <- dd99_newton_off(r_jacobian = 2L)
    comm$derivative_control$parallel <- parallel
    community_selection_gradient_jacobian(comm)
  }
  sequential <- jacobian(TRUE)
  with_two_workers("multicore", {
    expect_identical(jacobian(FALSE), sequential)
    # at the default equilibrium tolerance, points solved apart do move it
    expect_false(identical(jacobian(TRUE), sequential))
  })
})

test_that("a stencil point that merges two residents is refused under a plan as it is sequentially", {
  # 8e-4 apart at x ~ 0.5, inside two finite-difference steps (1e-3 |x| each)
  close <- dd99_model() |>
    community_add(trait_matrix(c(0.5, 0.5008), "x"), birth_rate = c(100, 100)) |>
    community_demography()
  for (type in c("multicore", "multisession")) {
    with_two_workers(type, {
      expect_error(community_selection_gradient_jacobian(close),
                   "residents 1 and 2 are within the finite-difference step")
    })
  }
})
