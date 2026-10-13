# Coexistence labels and their scouting in parameter space (R/boundary.R).
# Oracles in helper-boundary.R: GK98's 1 | 2 | 3 regions with the pair and
# triple in closed form, DD99's single boundary at sigma_C = sigma_K.

# ---- boundary_control --------------------------------------------------------

test_that("boundary_control has defaults and validates", {
  ctrl <- boundary_control()
  expect_equal(ctrl$evolve, "newton")
  expect_equal(ctrl$max_residents, 8L)
  expect_identical(boundary_control(ctrl), ctrl)
  expect_equal(boundary_control(list(max_residents = 3))$max_residents, 3L)
  expect_error(boundary_control(list(nope = 1)), "Unknown control parameters nope")
  expect_error(boundary_control(list(evolve = "none")), "should be one of")
  expect_error(boundary_control(list(max_residents = 0)), "positive whole number")
  expect_error(boundary_control(list(edge_tol = -1)), "single positive number")
  expect_error(boundary_control(list(fitness_tol = NA)), "single finite number")
})

test_that("assembler_control validates its evolution settings", {
  expect_equal(assembler_control()$evolve, "none")
  expect_equal(assembler_control(list(evolve = "newton"))$evolve, "newton")
  expect_error(assembler_control(list(evolve = "newton", birth_type = "stochastic", vcv = diag(1))),
               "evolve needs birth_type")
  expect_error(assembler_control(list(evolve = "newton", branch_offset = 0)), "single positive number")
  expect_error(assembler_control(list(max_residents = 0)), "at least 1")
  expect_error(assembler_control(list(edge_tol = -1)), "non-negative")
})

# ---- community_label ----------------------------------------------------------

test_that("GK98 supports one, two or three species as its patches move apart", {
  # d = 1 < sqrt(3/2): the generalist at 0 is a CSS; 1.5: a branching point,
  # the pair +-a*; 2.5 > 1.929: the pair's middle patch invadable, one species
  # per patch. Residents to 1e-6 of the closed forms (evolve_tol 1e-6 on the
  # gradient, whose slope is order one), densities exact.
  one <- community_label(gk98_kappa_start(1))
  l <- attr(one, "label")
  expect_equal(l[c("n", "status")], list(n = 1L, status = "stable"))
  expect_null(l$reason)
  expect_equal(as.numeric(one$traits), 0, tolerance = 1e-6)
  expect_equal(as.numeric(one$birth_rate), 3)
  expect_true(isTRUE(attr(one, "converged")))

  pair <- community_label(gk98_kappa_start(1.5))
  expect_equal(attr(pair, "label")[c("n", "status")], list(n = 2L, status = "stable"))
  a <- gk98_kappa_pair_root(1.5)
  expect_equal(sort(as.numeric(pair$traits)), c(-a, a), tolerance = 1e-6)
  expect_equal(as.numeric(pair$birth_rate), c(1.5, 1.5), tolerance = 1e-6)

  triple <- community_label(gk98_kappa_start(2.5))
  expect_equal(attr(triple, "label")[c("n", "status")], list(n = 3L, status = "stable"))
  oracle <- gk98_kappa_triple(2.5)
  ord <- order(as.numeric(triple$traits))
  expect_equal(as.numeric(triple$traits)[ord], oracle$traits, tolerance = 1e-6)
  expect_equal(as.numeric(triple$birth_rate)[ord], oracle$density, tolerance = 1e-6)
  # one round per species: the start evolved, then an immigrant admitted and
  # the coalition re-solved, until the maximum of fitness is at a resident
  expect_equal(vapply(list(one, pair, triple), function(x) attr(x, "label")$rounds, 0L), 1:3)
  # equilibrium solves: the start, each round's solve, classification and
  # immigrant; pinned so a change in cost is seen
  expect_identical(vapply(list(one, pair, triple), function(x) attr(x, "label")$evaluations, 0L),
                   c(7L, 27L, 54L))
})

test_that("evolving by the canonical equation reaches the same community", {
  triple <- community_label(gk98_kappa_start(2.5), boundary_control(list(evolve = "canonical")))
  l <- attr(triple, "label")
  expect_equal(l[c("n", "status")], list(n = 3L, status = "stable"))
  oracle <- gk98_kappa_triple(2.5)
  # the canonical route polishes to evolve_tol on the gradient, so the middle
  # resident settles within ~1e-4 of 0
  expect_equal(sort(as.numeric(triple$traits)), oracle$traits, tolerance = 1e-3)
  expect_gt(l$evaluations, 100L)
})

test_that("DD99 supports one species when competition is wider than the resource, else fills the box", {
  one <- community_label(dd99_start(sigma_C = 1.5))
  expect_equal(attr(one, "label")[c("n", "status")], list(n = 1L, status = "stable"))
  expect_equal(as.numeric(one$traits), 0, tolerance = 1e-8)
  # below the diagonal the branching point's landscape rises towards the
  # box's edges, where the immigrant is refused: the box sets the count
  edge <- community_label(dd99_start(sigma_C = 0.7))
  l <- attr(edge, "label")
  expect_equal(l[c("n", "status")], list(n = 1L, status = "edge"))
  expect_match(l$reason, "edge of the bounds")
  expect_equal(as.numeric(edge$traits), 0, tolerance = 1e-8)
  # a first immigrant into an empty community may land on the edge: the
  # empty DD99 community's fitness is flat, and evolution moves it to x0
  from_empty <- community_label(community_start(bounds(x = c(-2, 2)), trait_scale = "linear",
                                                harness = harness_dd99(sigma_C = 1.5)))
  expect_equal(attr(from_empty, "label")[c("n", "status")], list(n = 1L, status = "stable"))
  expect_equal(as.numeric(from_empty$traits), 0, tolerance = 1e-8)
})

test_that("a label is censored at the cap, and unresolved where no attractor is reached", {
  capped <- community_label(gk98_kappa_start(2.5), boundary_control(list(max_residents = 2)))
  l <- attr(capped, "label")
  expect_equal(l[c("n", "status")], list(n = 2L, status = "capped"))
  expect_match(l$reason, "beyond 2 residents")
  # the pair the immigrant settled into is a singular coalition (here an
  # asymmetric one: a specialist on an outer patch and a resident straddling
  # the other two), not a count the model decided
  expect_lt(max(abs(community_selection_gradient(capped)$selection_gradient)), 1e-6)
  expect_gt(abs(sum(capped$traits)), 0.5)

  rounds <- community_label(gk98_kappa_start(2.5), boundary_control(list(max_rounds = 1)))
  expect_equal(attr(rounds, "label")$status, "unresolved")
  expect_match(attr(rounds, "label")$reason, "still assembling after 1 round")

  # a singular strategy that is an ESS the dynamics flees: Newton finds it,
  # the classifier rejects it, and the dynamics run from the start instead,
  # here to the edge of the bounds
  goe <- community_label(garden_of_eden_start())
  l <- attr(goe, "label")
  expect_equal(l[c("n", "status")], list(n = 1L, status = "unresolved"))
  expect_match(l$reason, "canonical equation ended")
  expect_error(community_label(list()), "needs a community")
})

# ---- boundary_design and boundary_scout -----------------------------------------

test_that("boundary_design fills a box deterministically", {
  box <- bounds(d = c(1, 3), kappa = c(0.5, 3))
  des <- boundary_design(box, 7)
  expect_equal(dim(des), c(7L, 2L))
  expect_equal(colnames(des), c("d", "kappa"))
  expect_true(all(des[, "d"] >= 1 & des[, "d"] <= 3 & des[, "kappa"] >= 0.5 & des[, "kappa"] <= 3))
  expect_identical(boundary_design(box, 7), des)
  expect_identical(boundary_design(box, 3), des[1:3, ])
  # the sequence spreads: no two points share a coordinate
  expect_equal(length(unique(round(des[, 1], 10))), 7L)
  expect_error(boundary_design(box, 0), "positive whole number")
  expect_error(boundary_design(bounds(d = c(1, Inf)), 2), "finite")
})

test_that("scouting labels the DD99 plane by the diagonal", {
  by <- community_parameter_map(dd99_start(), c("sigma_C", "sigma_K"))
  des <- boundary_design(bounds(sigma_C = c(0.5, 1.5), sigma_K = c(0.8, 1.2)), 6)
  sc <- boundary_scout(by, des)
  expect_s3_class(sc, "boundary_labels")
  expect_equal(names(sc)[1:2], c("sigma_C", "sigma_K"))
  expect_equal(sc$sigma_C, des[, "sigma_C"])
  expect_equal(sc$n, rep(1L, 6))
  expect_equal(sc$status, ifelse(des[, "sigma_C"] > des[, "sigma_K"], "stable", "edge"))
  expect_true(all(vapply(sc$traits, function(x) abs(x) < 1e-6, logical(1))))
  expect_equal(unlist(sc$birth_rate), rep(500, 6))
  expect_equal(attr(sc, "evaluations"), sum(sc$evaluations))
  out <- paste(utils::capture.output(print(sc)), collapse = "\n")
  expect_match(out, "6 points over sigma_C, sigma_K")
  expect_match(out, "1 edge: 5 points")
  expect_match(out, "1 stable: 1 point\n")
  # one parameter as a vector, named by the map
  sc1 <- boundary_scout(community_parameter_map(dd99_start(), "sigma_C"), c(0.7, 1.3))
  expect_equal(names(sc1)[1], "sigma_C")
  expect_equal(sc1$status, c("edge", "stable"))
})

test_that("scouting validates its inputs and names a point that fails", {
  by <- community_parameter_map(dd99_start(), "sigma_C")
  expect_error(boundary_scout("sigma_C", 1), "must be a function")
  expect_error(boundary_scout(by, c(1, NA)), "finite parameter values")
  expect_error(boundary_scout(function(p) p, 1), "must return a community")
  broken <- function(p) if (p > 1) stop("no model here") else by(p)
  expect_error(boundary_scout(broken, c(0.7, 1.5)), "failed at p = 1.5: no model here")
})

# The map runs on the workers, which load the package and not the test
# helpers, so it is written with package functions only.
gk98_by_d_kappa <- function(p) {
  community_start(bounds(x = c(-4, 4)), trait_scale = "linear",
                  harness = harness_gk98(sigma = 1, mu = c(-p[1], 0, p[1]), K = c(1, p[2], 1)),
                  demography_control = demographic_step_control(list(equilibrium_extinct_birth_rate = 1e-8))) |>
    community_add(trait_matrix(0.5, "x"), birth_rate = 1)
}

test_that("scouting is the same, in answer and cost, under each plan", {
  by <- gk98_by_d_kappa
  des <- boundary_design(bounds(d = c(1, 2.6), kappa = c(0.5, 2)), 3)
  sequential <- boundary_scout(by, des)
  expect_setequal(sequential$status, "stable")
  expect_true(all(sequential$n %in% 1:3))
  for (type in c("multicore", "multisession")) {
    parallel <- with_two_workers(type, boundary_scout(by, des))
    expect_identical(as.data.frame(parallel), as.data.frame(sequential))
    expect_identical(attr(parallel, "evaluations"), attr(sequential, "evaluations"))
  }
})
