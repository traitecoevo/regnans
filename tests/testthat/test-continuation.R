# Parameter sensitivity dx*/dp of singular strategies and coalitions, and their
# continuation along one parameter (R/continuation.R). Oracles:
#
# JJ12: x* = x_opt - a sigma^2, so dx*/d(sigma^2) = -a, dx*/da = -sigma^2,
#   dx*/dx_opt = 1, dx*/dsigma = -2 a sigma.
# DD99, one resident: x* = x0 whatever the kernels, so dx*/dx0 = 1 and
#   dx*/dsigma_C = dx*/dsigma_K = 0; the multi-trait model the same per axis.
# DD99 pair: x0 -+ a*(sigma_C, sigma_K) in closed form, differentiated by hand
#   in helper-coalition.R (dd99_pair_root_gradient).
# GM99 pair: no closed form, so a central difference of the root of its
#   independently coded gradients (gm99_pair_root).
#
# Tolerances: the sensitivity is -J^-1 G with both Jacobians central
# differences on steps of 1e-3 (relative), good to ~1e-6 relative on these
# smooth models; 1e-5 leaves room.

dd99_comm <- function(sigma_C = 0.4, sigma_K = 1, x0 = 0) {
  community_start(bounds(x = c(-2, 2)), trait_scale = "linear",
                  harness = harness_dd99(x0 = x0, sigma_K = sigma_K, sigma_C = sigma_C))
}
pair_x0 <- function(...) trait_matrix(c(...), "x")

# a caller-written p -> community, the parameter being sigma^2 rather than a
# parameter of the harness
jj12_by_sigma2 <- function(a = 0.1, x_opt = 0.5) {
  function(p) {
    community_start(bounds(x = c(-3, 3)), trait_scale = "linear",
                    harness = harness_jj12(a = a, x_opt = x_opt, sigma = sqrt(p)))
  }
}

# ---- community_parameter_map -------------------------------------------------

test_that("community_parameter_map rebuilds the model at new parameter values", {
  comm <- dd99_comm() |>
    community_add(trait_matrix(0.2, "x"), birth_rate = 100) |>
    community_demography()
  by <- community_parameter_map(comm, c("x0", "sigma_K"))
  expect_equal(attr(by, "pars"), c("x0", "sigma_K"))
  expect_equal(community_parameter_values(comm, by), c(x0 = 0, sigma_K = 1))

  moved <- by(c(0.5, 1.2))
  expect_equal(community_parameter_values(moved, by), c(x0 = 0.5, sigma_K = 1.2))
  expect_equal(moved$harness$pars$x0, 0.5)
  expect_equal(moved$harness$pars$sigma_K, 1.2)
  expect_true(isTRUE(moved$harness$fd))
  # residents and bounds kept, nothing solved under the old model carried over
  expect_equal(moved$traits, comm$traits)
  expect_equal(moved$bounds, comm$bounds)
  expect_null(moved$fitness_function)
  expect_null(moved$selection_gradient)
  expect_null(attr(moved, "converged"))
  # the model really is the new one: its singular strategy is the new x0
  expect_equal(as.numeric(community_solve_singularity(moved)$traits), 0.5, tolerance = 1e-6)
  # and the community it was built from is untouched
  expect_equal(comm$harness$pars$x0, 0)
})

test_that("community_parameter_map gives a vector parameter one entry per element", {
  comm <- community_start(bounds(x1 = c(-2, 2), x2 = c(-2, 2)), trait_scale = "linear",
                          harness = harness_dd99_nd(x0 = c(0.2, -0.4)))
  by <- community_parameter_map(comm, c("x0", "r"))
  expect_equal(community_parameter_values(comm, by), c(`x0[1]` = 0.2, `x0[2]` = -0.4, r = 1))
  expect_equal(by(c(0.1, 0.3, 2))$harness$pars$x0, c(0.1, 0.3))
})

test_that("community_parameter_map refuses what it cannot vary", {
  comm <- dd99_comm()
  expect_error(community_parameter_map(comm, "sigma"), "no parameter sigma")
  expect_error(community_parameter_map(comm, c("x0", "x0")), "distinct parameters")
  expect_error(community_parameter_map(comm, "x0")(c(1, 2)), "must have 1 value")
  plant <- community_start(bounds(lma = c(0.01, 2)), harness = harness_plant())
  expect_error(community_parameter_map(plant, "x0"), "harness_explicit")
  expect_error(community_parameter_values(comm, function(p) comm), "give p")
  jj <- community_start(bounds(x = c(-3, 3)), trait_scale = "linear", harness = harness_jj12())
  expect_error(community_parameter_values(jj, community_parameter_map(comm, "x0")),
               "not built by this parameter map")
})

test_that("a community the map built is not mistaken for one solved under the old parameters", {
  # solved at x0 = 0.3, then rebuilt at x0 = 0.6: the old gradient (zero) and
  # convergence must not stand for the new model's
  sol <- community_solve_singularity(dd99_comm(x0 = 0.3), tol = 1e-10)
  by <- community_parameter_map(sol, "x0")
  moved <- community_demography(by(0.6))
  expect_equal(as.numeric(community_selection_gradient(moved)$selection_gradient), 0.3, tolerance = 1e-8)
  expect_warning(cl <- community_classify_singularity(by(0.6)), "not at a singular point")
  expect_equal(as.numeric(cl$selection_gradient), 0.3, tolerance = 1e-8)
  # the same holds for a resident added after a solve
  added <- community_add(sol, trait_matrix(0.9, "x"), birth_rate = 10)
  expect_null(added$selection_gradient)
  expect_null(attr(added, "converged"))
})

# ---- community_selection_gradient_parameter_jacobian ------------------------

test_that("the parameter Jacobian of the selection gradient matches DD99 in closed form", {
  # one resident at x with N = K(x): g = -r (x - x0) / sigma_K^2, so
  # dg/dx0 = r / sigma_K^2, dg/dsigma_K = 2 r (x - x0) / sigma_K^3, dg/dr = g / r
  comm <- dd99_comm(x0 = 0.3, sigma_K = 1.2) |>
    community_add(trait_matrix(-0.1, "x"), birth_rate = 100) |>
    community_demography()
  G <- community_selection_gradient_parameter_jacobian(
    comm, community_parameter_map(comm, c("x0", "sigma_K", "sigma_C")))
  expect_equal(dimnames(G), list("x", c("x0", "sigma_K", "sigma_C")))
  expect_equal(as.numeric(G), c(1 / 1.2^2, 2 * (-0.4) / 1.2^3, 0), tolerance = 1e-5)
  # two equilibrium solves per parameter
  expect_equal(attr(G, "evaluations"), 6L)
})

test_that("the parameter Jacobian validates its parameter", {
  comm <- dd99_comm() |> community_add(trait_matrix(0, "x"), birth_rate = 100)
  expect_error(community_selection_gradient_parameter_jacobian(comm, function(p) comm),
               "Give p")
  expect_error(community_selection_gradient_parameter_jacobian(comm, "x0", p = 1),
               "must be a function")
  wrong <- function(p) community_start(bounds(y = c(-1, 1)), trait_scale = "linear",
                                       harness = harness_dd99(trait_name = "y"))
  expect_error(community_selection_gradient_parameter_jacobian(comm, wrong, p = 1),
               "same traits|with the traits x")
  expect_error(community_selection_gradient_parameter_jacobian(dd99_comm(), function(p) comm, p = 1),
               "at least one resident")
})

# ---- community_parameter_sensitivity ----------------------------------------

test_that("JJ12: dx*/dsigma^2 = -a, through a caller-written parameter", {
  a <- 0.1
  by <- jj12_by_sigma2(a = a)
  sol <- community_solve_singularity(by(1.96), tol = 1e-10)
  S <- community_parameter_sensitivity(sol, by, p = 1.96)
  expect_equal(as.numeric(S), -a, tolerance = 1e-5)
  expect_equal(dimnames(S), list("x", "p"))
  # 2k for the resident Jacobian (its centre taken from the community), two
  # for the parameter
  expect_equal(attr(S, "evaluations"), 4L)
  expect_equal(dim(attr(S, "jacobian")), c(1L, 1L))
  expect_equal(dim(attr(S, "parameter_jacobian")), c(1L, 1L))
})

test_that("JJ12: the sensitivity to each harness parameter at once", {
  a <- 0.1; sigma <- 1.4
  comm <- community_start(bounds(x = c(-3, 3)), trait_scale = "linear",
                          harness = harness_jj12(a = a, x_opt = 0.5, sigma = sigma))
  sol <- community_solve_singularity(comm, tol = 1e-10)
  S <- community_parameter_sensitivity(sol, community_parameter_map(comm, c("a", "x_opt", "sigma")))
  expect_equal(colnames(S), c("a", "x_opt", "sigma"))
  expect_equal(as.numeric(S), c(-sigma^2, 1, -2 * a * sigma), tolerance = 1e-5)
})

test_that("DD99: dx*/dx0 = 1 and the kernels do not move x*", {
  comm <- dd99_comm(x0 = 0.3)
  sol <- community_solve_singularity(comm, tol = 1e-10)
  S <- community_parameter_sensitivity(sol, community_parameter_map(comm, c("x0", "sigma_C", "sigma_K")))
  expect_equal(S[1, "x0"], 1, tolerance = 1e-5)
  expect_lt(max(abs(S[1, c("sigma_C", "sigma_K")])), 1e-6)
  expect_equal(attr(S, "evaluations"), 2L + 2L * 3L)
})

test_that("two-trait DD99: dx*/dx0 is the identity", {
  comm <- community_start(bounds(x1 = c(-2, 2), x2 = c(-2, 2)), trait_scale = "linear",
                          harness = harness_dd99_nd(x0 = c(0.2, -0.4), sigma_C = c(0.4, 1.5)))
  sol <- community_solve_singularity(comm, tol = 1e-10)
  S <- community_parameter_sensitivity(sol, community_parameter_map(comm, "x0"))
  expect_equal(dimnames(S), list(c("x1", "x2"), c("x0[1]", "x0[2]")))
  expect_equal(unname(S[, ]), diag(2), tolerance = 1e-5)
})

test_that("the DD99 pair moves as its closed form a*(sigma_C, sigma_K), differentiated by hand", {
  # the hand derivative itself, against a difference of the closed form
  h <- 1e-6
  expect_equal(dd99_pair_root_gradient(0.4, 1),
               c(sigma_C = (dd99_pair_root(0.4 + h, 1) - dd99_pair_root(0.4 - h, 1)) / (2 * h),
                 sigma_K = (dd99_pair_root(0.4, 1 + h) - dd99_pair_root(0.4, 1 - h)) / (2 * h)),
               tolerance = 1e-8)

  comm <- dd99_comm(x0 = 0.3)
  sol <- community_solve_singularity(comm, x0 = pair_x0(0, 0.6), tol = 1e-10)
  S <- community_parameter_sensitivity(sol, community_parameter_map(comm, c("x0", "sigma_C", "sigma_K")))
  da <- dd99_pair_root_gradient(0.4, 1)
  expect_equal(dimnames(S), list(c("x[1]", "x[2]"), c("x0", "sigma_C", "sigma_K")))
  # the pair is x0 -+ a*: both residents follow x0, and spread apart as the
  # kernels widen
  expect_equal(unname(S[, ]), unname(rbind(c(1, -da), c(1, da))), tolerance = 1e-5)
  # 2mk for the coalition Jacobian, two per parameter
  expect_equal(attr(S, "evaluations"), 4L + 6L)
})

test_that("the asymmetric GM99 pair on a log trait scale moves as its root does", {
  alpha <- 7; h <- 1e-4
  da <- (gm99_pair_root(alpha = alpha + h, beta = 15) -
           gm99_pair_root(alpha = alpha - h, beta = 15)) / (2 * h)
  gm <- community_start(bounds(x = c(0.1, 0.95)), harness = harness_gm99(alpha = alpha, beta = 15))
  sol <- community_solve_singularity(gm, x0 = pair_x0(0.5, 0.8), tol = 1e-10)
  S <- community_parameter_sensitivity(sol, community_parameter_map(gm, "alpha"))
  # raw trait units whatever the trait scale; the residents move differently
  expect_equal(as.numeric(S), da, tolerance = 1e-4)
  expect_gt(abs(diff(as.numeric(S))), 1e-3)
})

test_that("a singular resident Jacobian is refused as a fold, not answered", {
  # s(y) = theta (y - x)^3: every resident is singular and dg/dx = 0
  cubic <- harness_explicit(
    fitness = function(x_mut, x_res, n_res, pars) pars$theta * (x_mut - x_res[1])^3,
    equilibrium = function(x_res, pars) rep(1, length(x_res)),
    fitness_gradient = function(x_mut, x_res, n_res, pars)
      matrix(3 * pars$theta * (x_mut - x_res[1])^2, ncol = 1L),
    pars = list(theta = 1), trait_names = "x", label = "cubic")
  comm <- community_start(bounds(x = c(-1, 1)), trait_scale = "linear", harness = cubic) |>
    community_add(trait_matrix(0.2, "x"), birth_rate = 1) |>
    community_demography()
  expect_error(community_parameter_sensitivity(comm, community_parameter_map(comm, "theta")),
               "resident Jacobian is singular")
})

test_that("the sensitivity is taken at the parameters the community was built at", {
  # a map built at sigma_C = 0.4 asked about a pair solved at sigma_C = 0.6
  comm <- dd99_comm()
  by <- community_parameter_map(comm, "sigma_C")
  pair <- community_solve_singularity(by(0.6), x0 = pair_x0(-0.3, 0.6), tol = 1e-10)
  S <- community_parameter_sensitivity(pair, by)
  da <- dd99_pair_root_gradient(0.6, 1)[["sigma_C"]]
  expect_equal(as.numeric(S), c(-da, da), tolerance = 1e-5)
  expect_equal(as.numeric(community_parameter_sensitivity(pair, by, p = 0.6)), as.numeric(S))
  expect_error(community_parameter_sensitivity(pair, by, p = 0.4), "built at sigma_C = 0.6")
})

test_that("the sensitivity warns away from a singularity", {
  comm <- dd99_comm() |>
    community_add(trait_matrix(0.5, "x"), birth_rate = 100) |>
    community_demography()
  expect_warning(community_parameter_sensitivity(comm, community_parameter_map(comm, "x0")),
                 "not at a singular point")
  attr(comm, "converged") <- FALSE
  expect_warning(community_parameter_sensitivity(comm, community_parameter_map(comm, "x0")),
                 "marked as not converged")
})

# ---- community_continue_singularity -----------------------------------------

test_that("JJ12 is followed along sigma^2 exactly, the predictor landing on x*", {
  a <- 0.1; x_opt <- 0.5
  by <- jj12_by_sigma2(a = a, x_opt = x_opt)
  p <- c(1, 1.5, 2, 2.5)
  path <- community_continue_singularity(community_solve_singularity(by(1)), by, p, tol = 1e-10)
  expect_s3_class(path, "singularity_path")
  expect_equal(path$p, p)
  expect_equal(as.numeric(path$traits), x_opt - a * p, tolerance = 1e-8)
  expect_equal(as.numeric(path$sensitivity), rep(-a, 4), tolerance = 1e-5)
  # x* is linear in sigma^2, so the tangent predicts it to the accuracy of
  # the sensitivity
  expect_equal(path$predicted[-1, ], path$traits[-1, ], tolerance = 1e-5)
  expect_equal(path$classification, rep("CSS", 4))
  expect_equal(nrow(path$changes), 0L)
  expect_null(path$stopped)
  expect_equal(attr(path, "evaluations"), sum(path$evaluations))
  expect_equal(as.numeric(path$community$traits), x_opt - a * 2.5, tolerance = 1e-8)
})

test_that("DD99 along sigma_C: the branching point becomes a CSS as sigma_C passes sigma_K", {
  comm <- dd99_comm(sigma_C = 0.7, x0 = 0.3)
  path <- community_continue_singularity(community_solve_singularity(comm),
                                         community_parameter_map(comm, "sigma_C"),
                                         c(0.7, 0.9, 1.1, 1.3))
  expect_equal(as.numeric(path$traits), rep(0.3, 4), tolerance = 1e-6)
  expect_equal(path$classification, c("branching point", "branching point", "CSS", "CSS"))
  expect_equal(path$changes,
               data.frame(p_before = 0.9, p_after = 1.1, from = "branching point", to = "CSS",
                          stringsAsFactors = FALSE))
  # the Hessian along the path is the closed form r (1/sigma_C^2 - 1/sigma_K^2)
  expect_equal(vapply(path$classifications, function(cl) as.numeric(cl$hessian), 0),
               1 / c(0.7, 0.9, 1.1, 1.3)^2 - 1, tolerance = 1e-5)
  out <- paste(utils::capture.output(print(path)), collapse = "\n")
  expect_match(out, "4 points, p from 0.7 to 1.3")
  expect_match(out, "branching point \\(p 0.7 to 0.9\\), CSS \\(p 1.1 to 1.3\\)")
})

test_that("the DD99 pair is followed along sigma_K as its closed form", {
  comm <- dd99_comm(x0 = 0.3)
  p <- c(1, 1.1, 1.2, 1.3)
  start <- community_solve_singularity(comm, x0 = pair_x0(0, 0.6), tol = 1e-10)
  path <- community_continue_singularity(start, community_parameter_map(comm, "sigma_K"), p,
                                         tol = 1e-10)
  a <- vapply(p, function(s) dd99_pair_root(0.4, s), 0)
  da <- vapply(p, function(s) dd99_pair_root_gradient(0.4, s)[["sigma_K"]], 0)
  expect_equal(colnames(path$traits), c("x[1]", "x[2]"))
  expect_equal(unname(path$traits), unname(cbind(0.3 - a, 0.3 + a)), tolerance = 1e-8)
  expect_equal(unname(path$sensitivity), unname(cbind(-da, da)), tolerance = 1e-5)
  expect_equal(dim(path$birth_rate), c(4L, 2L))
  expect_equal(path$classification, rep("branching point", 4))
  expect_match(paste(utils::capture.output(print(path)), collapse = "\n"), "coalition of 2")
  # each point after the first: the classification's 2mk + m solves, two for
  # the sensitivity, and a corrector started from the last point's Jacobian
  # (computing its own would cost 2mk = 4 more)
  corrector <- path$evaluations[-1] - (4L + 2L + 2L)
  expect_true(all(corrector <= 5L))
})

test_that("the asymmetric GM99 pair branches again as competition grows more asymmetric", {
  # followed along alpha on a log trait scale, the small-seeded resident's
  # curvature changes sign between alpha = 7.5 and 7.75 (helper-coalition.R)
  gm <- community_start(bounds(x = c(0.1, 0.95)), harness = harness_gm99(alpha = 7, beta = 15))
  start <- community_solve_singularity(gm, x0 = pair_x0(0.5, 0.8), tol = 1e-10)
  p <- c(7, 7.5, 7.75)
  path <- community_continue_singularity(start, community_parameter_map(gm, "alpha"), p, tol = 1e-10)
  expect_equal(unname(path$traits),
               do.call(rbind, lapply(p, function(a) gm99_pair_root(alpha = a, beta = 15))),
               tolerance = 1e-7)
  curv <- vapply(path$classifications, function(cl) vapply(cl$hessian, as.numeric, 0), numeric(2))
  expect_equal(curv[, 2:3], unname(cbind(gm99_coalition_curvature(path$traits[2, ], alpha = 7.5, beta = 15),
                                  gm99_coalition_curvature(path$traits[3, ], alpha = 7.75, beta = 15))),
               tolerance = 1e-5)
  expect_equal(path$classification, c("CSS", "CSS", "branching point"))
  expect_equal(path$classifications[[3]]$resident_evolutionarily_stable, c(FALSE, TRUE))
  expect_equal(path$changes$p_before, 7.5)
  expect_equal(path$changes$p_after, 7.75)
})

test_that("a path stops, with a warning, where the DD99 pair merges", {
  # a* -> 0 as sigma_C -> sigma_K = 1, a pitchfork: beyond it there is no
  # pair, and the corrector loses the coalition. (At sigma_C = 1 itself the
  # gradients are O(a^3), so a near-merged pair can pass as a root.)
  comm <- dd99_comm()
  start <- community_solve_singularity(comm, x0 = pair_x0(-0.3, 0.6), tol = 1e-10)
  expect_warning(
    path <- community_continue_singularity(start, community_parameter_map(comm, "sigma_C"),
                                           c(0.4, 0.6, 0.8, 0.95, 1.1), tol = 1e-10),
    "stopped at p = 1.1")
  expect_equal(path$p, c(0.4, 0.6, 0.8, 0.95))
  expect_equal(path$stopped$p, 1.1)
  expect_type(path$stopped$reason, "character")
  a <- vapply(path$p, function(s) dd99_pair_root(s, 1), 0)
  expect_equal(unname(path$traits), unname(cbind(-a, a)), tolerance = 1e-8)
  expect_match(paste(utils::capture.output(print(path)), collapse = "\n"), "stopped at p = 1.1")
})

test_that("a log-scale path predicts on the log scale and clamps to the bounds", {
  # x* = 0.555, dx*/dalpha = 0.129: the step to alpha = 8.5 overshoots the upper
  # bound even on the log scale; clamped there, the corrector finds the root
  gm <- community_start(bounds(x = c(0.05, 0.95)), harness = harness_gm99(alpha = 6, beta = 15))
  by <- community_parameter_map(gm, "alpha")
  start <- community_solve_singularity(gm, tol = 1e-10)
  path <- community_continue_singularity(start, by, c(6, 8.5), tol = 1e-10)
  expect_null(path$stopped)
  expect_equal(path$predicted[2, 1], c(x = 0.95))
  expect_equal(as.numeric(path$traits[2, ]),
               as.numeric(community_solve_singularity(by(8.5), tol = 1e-10)$traits), tolerance = 1e-8)
  # a raw-scale step of the same size would have predicted a negative seed size
  S <- path$sensitivity[1, 1]
  expect_lt(start$traits[1, 1] - 3 * S * 2.5, 0)
  path <- community_continue_singularity(start, by, c(6, 6 - 2.5 * 3), tol = 1e-10)
  expect_gt(path$predicted[2, 1], 0)
})

test_that("at the pitchfork the near-merged pair is flagged, and the path stops beyond it", {
  # at sigma_C = sigma_K the corrector lands on a pair ~4e-4 apart, every
  # derivative O(a^2), so the classifier warns; beyond it there is no pair.
  # The sensitivity there is rounding noise (the Jacobian's entries are
  # ~1e-7), so how the next step fails -- a clamped prediction that merges
  # the pair, or a corrector that loses it -- depends on the platform
  comm <- dd99_comm(sigma_C = 0.95)
  start <- community_solve_singularity(comm, x0 = pair_x0(-0.3, 0.35), tol = 1e-10)
  warnings <- character(0)
  path <- withCallingHandlers(
    community_continue_singularity(start, community_parameter_map(comm, "sigma_C"),
                                   c(0.95, 1, 1.1), tol = 1e-10),
    warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    })
  expect_true(any(grepl("residents 1 and 2 are within 0.001 of the bounds' width", warnings)))
  expect_true(any(grepl("stopped at p = 1.1", warnings)))
  expect_equal(path$p, c(0.95, 1))
  expect_equal(path$stopped$p, 1.1)
})

test_that("a clamped prediction is named in the reason when the corrector then fails", {
  # x* = x0: the step from x0 = 0 to 3 predicts x = 3, clamped to the bound
  # at 2, and the singularity at 3 lies outside the bounds
  comm <- dd99_comm()
  start <- community_solve_singularity(comm, tol = 1e-10)
  expect_warning(
    path <- community_continue_singularity(start, community_parameter_map(comm, "x0"), c(0, 3)),
    "stopped at p = 3: Bounds do not include a singularity .*sensitivity of up to 1 at p = 0, was clamped to the bounds")
  expect_equal(path$p, 0)
})

test_that("a point that cannot be differentiated ends the path but keeps the points before it", {
  # DD99 below p = 1, then a model whose every point is singular with J = 0
  cubic <- harness_explicit(
    fitness = function(x_mut, x_res, n_res, pars) (x_mut - x_res[1])^3,
    equilibrium = function(x_res, pars) rep(1, length(x_res)),
    fitness_gradient = function(x_mut, x_res, n_res, pars) matrix(3 * (x_mut - x_res[1])^2, ncol = 1L),
    pars = list(), trait_names = "x", label = "cubic")
  by <- function(p) {
    h <- if (p < 1) harness_dd99(x0 = p) else cubic
    community_start(bounds(x = c(-2, 2)), trait_scale = "linear", harness = h)
  }
  start <- community_solve_singularity(by(0.2), tol = 1e-10)
  expect_warning(path <- community_continue_singularity(start, by, c(0.2, 0.5, 1.5)),
                 "stopped at p = 1.5: .*resident Jacobian is singular")
  expect_equal(path$p, c(0.2, 0.5))
  expect_equal(as.numeric(path$traits), c(0.2, 0.5), tolerance = 1e-6)
  expect_error(community_continue_singularity(community_solve_singularity(by(1.5)), by, 1.5),
               "could not classify the starting point")
})

test_that("classifier arguments reach the classifier, apart from the corrector's", {
  comm <- dd99_comm(sigma_C = 0.7)
  path <- community_continue_singularity(community_solve_singularity(comm),
                                         community_parameter_map(comm, "sigma_C"), c(0.7, 0.8),
                                         tol = 1e-9, classify = list(tol = 1e-3, speeds = "equal"))
  expect_equal(path$classifications[[1]]$tol, 1e-3)
  expect_error(community_continue_singularity(community_solve_singularity(comm),
                                              community_parameter_map(comm, "sigma_C"), 0.7,
                                              classify = list(solver = "newton")),
               "classify must be a named list")
})

test_that("community_continue_singularity validates its start", {
  comm <- dd99_comm()
  by <- community_parameter_map(comm, "sigma_C")
  expect_error(community_continue_singularity(comm, by, 0.4), "at least one resident")
  one <- community_solve_singularity(comm)
  expect_error(community_continue_singularity(one, by, numeric(0)), "finite parameter values")
  expect_error(community_continue_singularity(one, "sigma_C", 0.4), "must be a function")
  expect_error(community_continue_singularity(one, by, 0.4, solver = "bracket"))
  expect_error(community_continue_singularity(one, community_parameter_map(comm, c("sigma_C", "sigma_K")), 0.4),
               "follows one parameter; this map varies sigma_C, sigma_K")
  # no dimorphism exists when sigma_C > sigma_K
  wide <- dd99_comm(sigma_C = 1.5)
  expect_error(suppressWarnings(
    community_continue_singularity(wide |> community_add(pair_x0(-0.3, 0.6), birth_rate = c(100, 100)),
                                   community_parameter_map(wide, "sigma_C"), 1.5)),
    "could not solve for the starting point")
})
