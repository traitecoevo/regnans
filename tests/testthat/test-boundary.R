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

test_that("a label needs bounds finite on the trait scale, and an unbounded trait has no edge", {
  # the search, the edge rule and the label's keys all measure in fractions
  # of the box's width
  open <- community_start(bounds(x = c(-Inf, Inf)), trait_scale = "linear",
                          harness = harness_dd99(sigma_C = 1.5)) |>
    community_add(trait_matrix(0.5, "x"), birth_rate = 100)
  expect_error(community_label(open), "finite on the trait scale; these are x: -Inf to Inf")
  plane <- function(lower, scale) {
    community_start(bounds(x1 = c(lower, 2), x2 = c(0.5, 2)), trait_scale = scale,
                    harness = harness_dd99_nd(x0 = c(1, 1), sigma_K = c(1, 1), sigma_C = c(1.5, 1.5))) |>
      community_add(trait_matrix(c(1, 1), c("x1", "x2")), birth_rate = 100)
  }
  expect_error(community_label(plane(-Inf, "linear")), "finite on the trait scale; these are x1: -Inf to 2")
  expect_error(community_label(plane(0, "log")), "positive lower bound for a log-scale trait")
  # an unbounded trait has no edge to refuse a birth at
  ctrl <- assembler_control(list(edge_tol = 1e-3))
  expect_null(births_refused(plane(-Inf, "linear"), trait_matrix(c(-5, 1), c("x1", "x2")), ctrl))
  expect_equal(births_refused(plane(-Inf, "linear"), trait_matrix(c(-5, 1.9999), c("x1", "x2")), ctrl)$status,
               "edge")
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

test_that("a design's columns are matched to the map's parameters by name", {
  # the map applies p by position, so a design in another order is reordered
  # rather than built with its parameters swapped
  by <- community_parameter_map(dd99_start(), c("sigma_C", "sigma_K"))
  des <- boundary_design(bounds(sigma_C = c(0.5, 1.5), sigma_K = c(0.8, 1.2)), 4)
  swapped <- des[, c("sigma_K", "sigma_C")]
  expect_identical(as.data.frame(boundary_scout(by, swapped)), as.data.frame(boundary_scout(by, des)))
  colnames(swapped) <- c("sigma_K", "x0")
  expect_error(boundary_scout(by, swapped), "must be the map's parameters, sigma_C, sigma_K")
  # a ray's ends, each by its own names
  bp <- boundary_bisect(by, c(sigma_K = 1, sigma_C = 1.5), c(sigma_K = 1.2, sigma_C = 0.5),
                        boundary_control(list(bisect_tol = 1e-4)))
  expect_equal(names(bp)[2:3], c("sigma_C", "sigma_K"))
  expect_equal(bp$sigma_C, 1.5 - 0.5 / 1.2, tolerance = 1e-3)
  expect_equal(bp$sigma_C, bp$sigma_K, tolerance = 1e-3)
})

test_that("scouting validates its inputs and names a point that fails", {
  by <- community_parameter_map(dd99_start(), "sigma_C")
  expect_error(boundary_scout("sigma_C", 1), "must be a function")
  expect_error(boundary_scout(by, c(1, NA)), "finite parameter values")
  expect_error(boundary_scout(function(p) p, 1), "must return a community")
  broken <- function(p) if (p > 1) stop("no model here") else by(p)
  expect_error(boundary_scout(broken, c(0.7, 1.5)), "failed at p = 1.5: no model here")
})

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

# ---- boundary_bisect ------------------------------------------------------------

test_that("bisection along rays recovers DD99's boundary at sigma_C = sigma_K", {
  # wide competition (one species, stable) towards narrow (the box fills).
  # The label changes where the fitness at the box's edge, 1 - exp(-2 (1 /
  # sigma_C^2 - 1 / sigma_K^2)) ~ 4 (1 - sigma_C / sigma_K), passes fitness_tol
  # (1e-8): 2.5e-9 inside the diagonal, below the bracket of 1e-7 of a ray
  # of length ~1, so the crossing is the diagonal to 1e-6.
  by <- community_parameter_map(dd99_start(), c("sigma_C", "sigma_K"))
  to <- rbind(c(0.5, 1.2), c(1.5, 2), c(2, 0.8))
  bp <- boundary_bisect(by, c(1.5, 1), to, boundary_control(list(bisect_tol = 1e-7)))
  expect_s3_class(bp, "boundary_points")
  expect_equal(bp$ray, 1:2)
  expect_equal(bp$sigma_C, bp$sigma_K, tolerance = 1e-6)
  expect_equal(bp$sigma_C, c(1.5 - 0.5 / 1.2, 1.5), tolerance = 1e-6)
  expect_true(all(bp$t_hi - bp$t_lo <= 1e-7))
  expect_equal(bp$status_lo, c("stable", "stable"))
  expect_equal(bp$status_hi, c("edge", "edge"))
  expect_true(all(bp$resolved))
  # the residents either side, for a boundary solve to start from
  expect_equal(as.numeric(bp$traits_lo[[1]]), 0, tolerance = 1e-8)
  rays <- attr(bp, "rays")
  expect_equal(rays$changes, c(1L, 1L, 0L))
  expect_equal(rays$status_end, c("edge", "edge", "stable"))
  expect_equal(attr(bp, "evaluations"), sum(rays$evaluations))
  out <- paste(utils::capture.output(print(bp)), collapse = "\n")
  expect_match(out, "2 changes on 3 rays")
  expect_match(out, "ray 2 at sigma_C = 1.5, sigma_K = 1.5: 1 stable -> 1 edge")
})

test_that("bisection recovers GK98's nested boundaries, a pitchfork and a transversal crossing", {
  # From one species at d = 1 towards three. At the pitchfork the emerging
  # peak's fitness grows as ~2.67 (d - c1)^2 (helper-boundary.R's peak), so
  # the label changes sqrt(fitness_tol / 2.67) beyond c1: 6e-7 at
  # fitness_tol = 1e-12. At the transversal crossing s(0) grows at 0.27 per
  # unit d, so the change is within fitness_tol / 0.27 of c2, and the pair's
  # gradients solved to 1e-10 leave the peak's height within the same; the
  # bracket is 1e-7 of a ray of length 1.6 to 2.
  control <- boundary_control(list(bisect_tol = 1e-7, evolve_tol = 1e-10, fitness_tol = 1e-12))
  by <- gk98_by_d_kappa
  bp <- boundary_bisect(by, c(1, 1), rbind(c(2.6, 1), c(2.6, 2)), control)
  expect_equal(nrow(bp), 4L)
  expect_true(all(bp$resolved))
  expect_equal(bp$n_lo, c(1L, 2L, 1L, 2L))
  expect_equal(bp$n_hi, c(2L, 3L, 2L, 3L))
  expect_equal(unique(c(bp$status_lo, bp$status_hi)), "stable")
  oracle <- t(vapply(bp[["p[2]"]], gk98_kappa_boundaries, numeric(2)))
  d <- bp[["p[1]"]]
  pitchfork <- d[c(1, 3)] - oracle[c(1, 3), "c1"]
  expect_true(all(pitchfork > 0 & pitchfork < 2e-6))
  expect_equal(d[c(2, 4)], oracle[c(2, 4), "c2"], tolerance = 1e-7)
  # each side's community: the pair at the closed form, two species below
  # the transversal crossing and three beyond it
  a <- gk98_kappa_pair_root(d[2], 1)
  expect_equal(sort(as.numeric(bp$traits_lo[[2]])), c(-a, a), tolerance = 1e-6)
  expect_equal(length(bp$traits_hi[[2]]), 3L)
  # with the default fitness_tol (1e-8) the pitchfork moves out by
  # sqrt(1e-8 / 2.67) = 6.1e-5, as predicted
  coarse <- boundary_bisect(by, c(1, 1), c(1.6, 1), boundary_control(list(bisect_tol = 1e-7)))
  expect_equal(coarse[["p[1]"]] - sqrt(1.5), sqrt(1e-8 / 2.67), tolerance = 0.05)
})

test_that("an unresolved label stops a bracket rather than being bisected across", {
  # the Garden of Eden model has no attractor to assemble to anywhere, so
  # both ends of the ray are unresolved and nothing is refined
  goe <- function(p) garden_of_eden_start(c = p)
  bp <- boundary_bisect(goe, 1, 2)
  expect_equal(nrow(bp), 1L)
  expect_false(bp$resolved)
  expect_equal(bp$reason, "an end of the bracket is unresolved")
  expect_equal(c(bp$t_lo, bp$t_hi), c(0, 1))
})

test_that("bisection validates its rays", {
  by <- community_parameter_map(dd99_start(), c("sigma_C", "sigma_K"))
  expect_error(boundary_bisect("x", c(1, 1), c(2, 2)), "must be a function")
  expect_error(boundary_bisect(by, c(1, 1), c(2, 2, 2)), "one column per parameter")
  expect_error(boundary_bisect(by, c(1, 1), c(1, 1)), "end away from where it starts")
  broken <- function(p) if (p[1] > 1.2) stop("no model here") else by(p)
  expect_error(boundary_bisect(broken, c(1, 1), c(2, 1)), "failed on ray 1: .*no model here")
})

test_that("bisection is the same, in answer and cost, under each plan", {
  by <- community_parameter_map(dd99_start(), c("sigma_C", "sigma_K"))
  to <- rbind(c(0.5, 1.2), c(1.5, 2))
  control <- boundary_control(list(bisect_tol = 1e-4))
  sequential <- boundary_bisect(by, c(1.5, 1), to, control)
  for (type in c("multicore", "multisession")) {
    parallel <- with_two_workers(type, boundary_bisect(by, c(1.5, 1), to, control))
    expect_identical(as.data.frame(parallel), as.data.frame(sequential))
    expect_identical(attr(parallel, "rays"), attr(sequential, "rays"))
  }
})

# ---- boundary_sweep ---------------------------------------------------------------

test_that("a sweep flags the alternative attractors of a bistable model, and only there", {
  # helper-boundary.R: both monocultures resist invasion for |q| < 1/3.
  # Started cold at x = 0.5 the assembly reaches the +1 peak unless the -1
  # monoculture can exclude it (q < -1/3); carried forward from q = -0.5 it
  # stays on -1 until that peak can be invaded (q > 1/3), and backward from
  # 0.5 on +1 until q < -1/3.
  by <- community_parameter_map(bistable_start(), "q")
  q <- c(-0.5, -0.3, -0.1, 0.1, 0.3, 0.5)
  sw <- boundary_sweep(by, q)
  expect_s3_class(sw, "boundary_sweep")
  expect_equal(sw$q, q)
  expect_equal(sw$hysteresis, abs(q) < 1 / 3)
  expect_true(all(c(sw$cold_status, sw$forward_status, sw$backward_status) == "stable"))
  expect_true(all(c(sw$cold_n, sw$forward_n, sw$backward_n) == 1L))
  peak <- function(x) vapply(x, function(t) sign(as.numeric(t)), 0)
  expect_equal(peak(sw$cold_traits), c(-1, 1, 1, 1, 1, 1))
  expect_equal(peak(sw$forward_traits), c(-1, -1, -1, -1, -1, 1))
  expect_equal(peak(sw$backward_traits), c(-1, 1, 1, 1, 1, 1))
  expect_equal(attr(sw, "evaluations"), sum(sw$evaluations))
  expect_match(paste(utils::capture.output(print(sw)), collapse = "\n"), "6 points, 4 with hysteresis")
})

test_that("a sweep finds no hysteresis where the stable community is unique", {
  sw <- boundary_sweep(gk98_by_d_kappa, cbind(c(1, 1.5, 2.2, 2.5), 1))
  expect_false(any(sw$hysteresis))
  expect_equal(sw$cold_n, c(1L, 2L, 3L, 3L))
  expect_error(boundary_sweep(gk98_by_d_kappa, cbind(1, 1)), "at least two points")
})

test_that("an unresolved route leaves hysteresis unknown rather than flagging it", {
  base <- dd99_start()
  tf <- community_trait_transform(base)
  key <- function(row) {
    comm <- base
    comm$traits <- trait_matrix(row$traits, "x")
    label_key(comm, tf)
  }
  route <- function(n, status, x) list(n = n, status = status, traits = matrix(x, ncol = 1))
  pair <- route(2L, "stable", c(-0.5, 0.5))
  lost <- route(1L, "unresolved", 0)
  expect_false(boundary_hysteresis(list(pair, pair, pair), key, 1e-3))
  expect_true(is.na(boundary_hysteresis(list(lost, pair, pair), key, 1e-3)))
  expect_true(is.na(boundary_hysteresis(list(lost, lost, pair), key, 1e-3)))
  # two resolved routes that disagree, in label or community, are hysteresis
  expect_true(boundary_hysteresis(list(lost, pair, route(1L, "stable", 0)), key, 1e-3))
  expect_true(boundary_hysteresis(list(pair, pair, route(2L, "stable", c(-0.5, 0.9))), key, 1e-3))
})

test_that("every route of a sweep searches the box the cold label searches", {
  # a harness that narrows an empty community's box to (-2, 0.5), as plant's
  # fundamental niche does: the carried routes, which start with residents,
  # must search that box too, not the map's (-2, 2)
  seen <- new.env()
  seen$max <- -Inf
  narrowed <- function(p) {
    h <- harness_dd99(sigma_C = p, sigma_K = 1)
    h$fns$viable_bounds <- function(community) {
      community$bounds[1, ] <- c(-2, 0.5)
      community
    }
    update <- h$fns$update_fitness_function
    h$fns$update_fitness_function <- function(community) {
      out <- update(community)
      f <- out$fitness_function
      out$fitness_function <- function(x) {
        seen$max <- max(seen$max, x)
        f(x)
      }
      out
    }
    community_start(bounds(x = c(-2, 2)), trait_scale = "linear", harness = h)
  }
  sw <- boundary_sweep(narrowed, c(1.4, 1.5, 1.6))
  expect_lte(seen$max, 0.5)
  expect_false(any(sw$hysteresis))
  expect_identical(sw$forward_traits, sw$cold_traits)
})

# ---- community_solve_boundary ------------------------------------------------------

# the community on the side of a bisected crossing with fewer species, and
# the newcomer on the other side: the trait farthest from those residents
boundary_seed <- function(parameter, row) {
  p <- unlist(row[setdiff(names(row), c("ray", "t_lo", "t_hi", "n_lo", "status_lo", "n_hi", "status_hi",
                                        "resolved", "reason", "traits_lo", "birth_rate_lo", "traits_hi",
                                        "birth_rate_hi"))])
  lo <- as.numeric(row$traits_lo[[1]])
  hi <- as.numeric(row$traits_hi[[1]])
  comm <- parameter(p) |> community_clear_residents() |>
    community_add(trait_matrix(lo, "x"), birth_rate = row$birth_rate_lo[[1]])
  list(community = comm, p = p, y = hi[which.max(vapply(hi, function(h) min(abs(h - lo)), 0))])
}

test_that("Newton on the augmented system lands on GK98's 2 -> 3 boundary from a bisected seed", {
  # the pair +-a* with a peak at y* = 0 exactly marginal: d = c2 and a* in
  # closed form. With every residual solved to 1e-10 and the boundary's
  # slopes order 0.1-1, d*, a* and y* are good to ~1e-9.
  by <- gk98_by_d_kappa
  bp <- boundary_bisect(by, c(1, 1), c(2.6, 1), boundary_control(list(bisect_tol = 1e-4)))
  seed <- boundary_seed(by, bp[bp$n_lo == 2L, ])
  at <- community_solve_boundary(seed$community, by, y = seed$y, p = seed$p, tol = 1e-10)
  b <- attr(at, "boundary")
  c2 <- gk98_kappa_boundaries(1)[["c2"]]
  expect_true(b$converged)
  expect_equal(b$type, "invasion")
  expect_equal(b$free, "p[1]")
  expect_equal(b$p, c(`p[1]` = c2, `p[2]` = 1), tolerance = 1e-9)
  expect_equal(sort(as.numeric(at$traits)), c(-1, 1) * gk98_kappa_pair_root(c2), tolerance = 1e-9)
  expect_lt(abs(b$y), 1e-9)
  expect_lt(max(abs(b$residual)), 1e-10)
  expect_equal(names(b$residual), c("g:x[1]", "g:x[2]", "ds/dx(y)", "s(y)"))
  # quadratic from the bracket's midpoint: one Jacobian (2mk + 2 solves)
  # and a residual per step
  expect_identical(b$iterations, 2L)
  expect_identical(b$evaluations, 9L)
  # the boundary is what it says: 1e-3 either side the labels are 2 and 3
  label <- function(d) attr(community_label(by(c(d, 1))), "label")[c("n", "status")]
  expect_equal(label(c2 - 1e-3), list(n = 2L, status = "stable"))
  expect_equal(label(c2 + 1e-3), list(n = 3L, status = "stable"))
})

test_that("every column of the boundary Jacobian matches a finer difference of the residual", {
  # The layer's columns are central differences at the derivative steps
  # (1e-3 relative), good to ~1e-6; the reference is a central difference of
  # the whole residual at 1e-5, with GK98's equilibrium iterated to 1e-12,
  # good to ~1e-7. Away from the boundary, so no column is degenerate.
  by <- gk98_by_d_kappa
  a <- gk98_kappa_pair_root(1.8, 1.3)
  pair <- by(c(1.8, 1.3)) |> community_clear_residents() |>
    community_add(trait_matrix(c(-a, a + 0.05), "x"), birth_rate = c(1.6, 1.6))
  p <- c(`p[1]` = 1.8, `p[2]` = 1.3)
  for (free in 1:2) {
    residual <- boundary_residual_fn(pair, by, p, free, "invasion", 1L)
    u <- c(-a, a + 0.05, 0.1, p[[free]])
    J <- community_boundary_jacobian(residual, u)
    reference <- oracle_jacobian(residual, u, h = 1e-5)
    expect_equal(dim(J), c(4L, 4L))
    expect_equal(colnames(J), c("x[1]", "x[2]", "y:x", names(p)[free]))
    expect_lt(max(abs(unname(J) - reference)), 2e-6)
    # 2mk solves for the residents and two for the parameter
    expect_identical(attr(J, "evaluations"), 6L)
  }
  # the branching residual's: the determinant of a resident's Hessian
  dd <- dd99_start(sigma_C = 0.9) |> community_clear_residents() |>
    community_add(trait_matrix(c(-0.3, 0.4), "x"), birth_rate = c(200, 200))
  by_dd <- community_parameter_map(dd, c("sigma_C", "sigma_K"))
  residual <- boundary_residual_fn(dd, by_dd, community_parameter_values(dd, by_dd), 1L, "branching", 2L)
  u <- c(-0.3, 0.4, 0.9)
  J <- community_boundary_jacobian(residual, u)
  expect_equal(rownames(J), c("g:x[1]", "g:x[2]", "det(H[2])"))
  expect_lt(max(abs(unname(J) - oracle_jacobian(residual, u, h = 1e-5))), 2e-6 * max(abs(J)))
})

test_that("the branching type lands on DD99's and GK98's pitchforks", {
  # DD99: x* = x0 = 0 and H = 1/sigma_C^2 - 1/sigma_K^2, zero at sigma_C =
  # sigma_K exactly. Seeded from bisection's stable side.
  by <- community_parameter_map(dd99_start(), c("sigma_C", "sigma_K"))
  bp <- boundary_bisect(by, c(1.5, 1), c(0.5, 1), boundary_control(list(bisect_tol = 1e-3)))
  seed <- by(c(bp$sigma_C, 1)) |> community_clear_residents() |>
    community_add(trait_matrix(bp$traits_lo[[1]] + 0.1, "x"), birth_rate = 400)
  at <- community_solve_boundary(seed, by, type = "branching", free = "sigma_C", tol = 1e-10)
  b <- attr(at, "boundary")
  expect_true(b$converged)
  expect_null(b$y)
  expect_equal(b$p, c(sigma_C = 1, sigma_K = 1), tolerance = 1e-9)
  expect_lt(abs(as.numeric(at$traits)), 1e-9)
  expect_equal(attr(community_label(by(c(1 + 1e-3, 1))), "label")$status, "stable")
  expect_equal(attr(community_label(by(c(1 - 1e-3, 1))), "label")$status, "edge")

  # GK98 at kappa = 2: S''(0) = 2 d^2 / (2 + kappa) - 1, zero at d = sqrt(2)
  gk <- gk98_by_d_kappa
  one <- gk(c(1.3, 2)) |> community_clear_residents() |> community_add(trait_matrix(0.05, "x"), birth_rate = 4)
  at <- community_solve_boundary(one, gk, type = "branching", p = c(1.3, 2), tol = 1e-10)
  expect_equal(attr(at, "boundary")$p[[1]], sqrt(2), tolerance = 1e-9)
  label <- function(d) attr(community_label(gk(c(d, 2))), "label")$n
  expect_equal(c(label(sqrt(2) - 1e-3), label(sqrt(2) + 1e-3)), c(1L, 2L))
})

test_that("a marginal mutant at a resident is refused, before and after the solve", {
  by <- gk98_by_d_kappa
  a <- gk98_kappa_pair_root(1.95)
  pair <- by(c(1.95, 1)) |> community_clear_residents() |>
    community_add(trait_matrix(c(-a, a), "x"), birth_rate = c(1.5, 1.5))
  # the triple handed over as the residents, its middle resident as y
  triple <- by(c(1.95, 1)) |> community_clear_residents() |>
    community_add(trait_matrix(c(-a, 0, a), "x"), birth_rate = c(1.4, 0.2, 1.4))
  expect_error(community_solve_boundary(triple, by, y = 0, p = c(1.95, 1)),
               "y starts at resident 2: .*side with fewer species")
  # from between a resident and the middle peak Newton runs onto the
  # resident, where both conditions hold for any spacing
  expect_warning(out <- community_solve_boundary(pair, by, y = 0.62, p = c(1.95, 1), tol = 1e-10),
                 "converged onto a resident")
  expect_false(attr(out, "boundary")$converged)
})

test_that("residents merged in a boundary solve are outside its domain", {
  # as for the singularity solver: no solve, a non-finite residual for Newton
  # to back off from, and a start there refused
  by <- gk98_by_d_kappa
  a <- gk98_kappa_pair_root(1.95)
  pair <- by(c(1.95, 1)) |> community_clear_residents() |>
    community_add(trait_matrix(c(-a, a), "x"), birth_rate = c(1.5, 1.5))
  residual <- boundary_residual_fn(pair, by, c(`p[1]` = 1.95, `p[2]` = 1), 1L, "invasion", 1L)
  r <- residual(c(0.3, 0.3, 1.5, 1.95))
  expect_length(r, 4L)
  expect_true(all(is.na(r)))
  expect_identical(attr(residual, "evaluations")(), 0L)
  expect_error(community_boundary_jacobian(residual, c(0.3, 0.3, 1.5, 1.95)), "not finite")
  merged <- by(c(1.95, 1)) |> community_clear_residents() |>
    community_add(trait_matrix(c(0.3, 0.3), "x"), birth_rate = c(1.5, 1.5))
  expect_error(community_solve_boundary(merged, by, y = 1.5, p = c(1.95, 1)),
               "residents 1 and 2 start within the finite-difference step")
})

test_that("a boundary point whose marginal mutant is not a peak warns", {
  # helper-boundary.R's dip: the minimum near the dip crosses zero as the dip
  # fills; Newton finds that point, which is no new peak
  start <- dip_start(0.95)
  by <- community_parameter_map(start, "h")
  expect_warning(at <- community_solve_boundary(start, by, y = 1.5, tol = 1e-10),
                 "not a strict maximum of invasion fitness")
  b <- attr(at, "boundary")
  expect_true(b$converged)
  expect_gt(as.numeric(community_fitness_hessian(at, matrix(b$y, 1))), 0)
})

test_that("community_solve_boundary validates its arguments and reports a solve that stops short", {
  by <- gk98_by_d_kappa
  a <- gk98_kappa_pair_root(1.95)
  pair <- by(c(1.95, 1)) |> community_clear_residents() |>
    community_add(trait_matrix(c(-a, a), "x"), birth_rate = c(1.5, 1.5))
  expect_error(community_solve_boundary(pair, by, p = c(1.95, 1)), "needs y")
  expect_error(community_solve_boundary(pair, by, y = c(0, 1), p = c(1.95, 1)), "one value per trait")
  expect_error(community_solve_boundary(pair, by, y = 0, type = "branching", p = c(1.95, 1)), "takes no y")
  expect_error(community_solve_boundary(pair, by, type = "branching", which = 3, p = c(1.95, 1)),
               "one of the 2 residents")
  expect_error(community_solve_boundary(pair, by, y = 0, free = "kappa", p = c(1.95, 1)),
               "free must be one of the parameters, p\\[1\\], p\\[2\\]")
  expect_error(community_solve_boundary(community_clear_residents(pair), by, y = 0, p = c(1.95, 1)),
               "at least one")
  expect_warning(out <- community_solve_boundary(pair, by, y = 0.05, p = c(1.95, 1), maxit = 1),
                 "did not converge")
  expect_false(attr(out, "boundary")$converged)
})

test_that("a boundary solve is the same, in answer and cost, under each plan", {
  by <- gk98_by_d_kappa
  a <- gk98_kappa_pair_root(1.95)
  pair <- by(c(1.95, 1)) |> community_clear_residents() |>
    community_add(trait_matrix(c(-a, a), "x"), birth_rate = c(1.5, 1.5))
  sequential <- community_solve_boundary(pair, by, y = 0.05, p = c(1.95, 1), tol = 1e-10)
  for (type in c("multicore", "multisession")) {
    parallel <- with_two_workers(type, community_solve_boundary(pair, by, y = 0.05, p = c(1.95, 1), tol = 1e-10))
    expect_identical(attr(parallel, "boundary"), attr(sequential, "boundary"))
    expect_identical(parallel$traits, sequential$traits)
  }
})
