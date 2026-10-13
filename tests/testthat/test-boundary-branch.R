# Coexistence boundaries followed through two parameters (R/boundary_branch.R).
# Oracles in helper-boundary.R: GK98's 2 -> 3 boundary d = c2(kappa) and its
# end where the pair merges, the 1 -> 2 locus d = sqrt((2 + kappa) / 2) and
# its turn to subcritical at kappa = 4; DD99's branching locus through a map
# that makes it the unit circle.

# sigma_C = 1 / |p| with sigma_K = 1: DD99's resident at 0 branches iff
# sigma_C < sigma_K, so its branching locus is the unit circle. Package
# functions only, for the workers of a plan.
dd99_circle <- function(p) {
  community_start(bounds(x = c(-2, 2)), trait_scale = "linear",
                  harness = harness_dd99(sigma_C = 1 / sqrt(p[1]^2 + p[2]^2), sigma_K = 1)) |>
    community_add(trait_matrix(0.3, "x"), birth_rate = 400)
}

# the GK98 pair at (d, kappa) = (1.93, 1), near its 2 -> 3 boundary
gk98_pair_start <- function() {
  a <- gk98_kappa_pair_root(1.93, 1)
  gk98_by_d_kappa(c(1.93, 1)) |> community_clear_residents() |>
    community_add(trait_matrix(c(-a, a), "x"), birth_rate = c(1.5, 1.5))
}

# Near where the GK98 boundaries end every curvature of the landscape is
# O(1e-6), so the labels that verify a branch there need gradients and
# fitness resolved to 1e-12; at the defaults a near-neutral coalition passes.
tight <- boundary_control(list(evolve_tol = 1e-12, fitness_tol = 1e-12))

test_that("boundary_control validates the continuation settings", {
  ctrl <- boundary_control()
  expect_equal(c(ctrl$ds_min, ctrl$ds, ctrl$ds_max), c(1e-5, 0.05, 0.2))
  expect_error(boundary_control(list(ds = 1)), "ds_min <= ds <= ds_max")
  expect_error(boundary_control(list(tangent_tol = 1)), "between 0 and 1")
  expect_error(boundary_control(list(max_steps = 0)), "positive whole number")
  expect_error(boundary_control(list(merge_tol = 0)), "single positive number")
})

test_that("DD99's branching locus is followed round the unit circle and closes", {
  # a closed curve, traced from (0.9, 0.1) back to itself: every point on the
  # circle to the corrector's tolerance, the resident at x0 = 0, the
  # arclength 2 pi less the last chord, and nothing to flag
  br <- community_continue_boundary(dd99_circle(c(0.9, 0.1)), dd99_circle, type = "branching",
                                    p = c(0.9, 0.1))
  expect_s3_class(br, "regnans_branch")
  expect_equal(unname(br$termination), "closed")
  expect_equal(nrow(br$events), 0L)
  pts <- br$points
  expect_lt(max(abs(sqrt(pts[["p[1]"]]^2 + pts[["p[2]"]]^2) - 1)), 1e-8)
  expect_lt(max(abs(pts$x)), 1e-10)
  gap <- sqrt(sum((br$u[nrow(br$u), ] - br$u[1, ])^2))
  expect_equal(max(pts$s) + gap, 2 * pi, tolerance = 1e-3)
  expect_lt(gap, 0.2)
  # labels either side, verified cold: one species, stable inside the circle
  # (a CSS), the box's edge outside (a branching point)
  v <- br$verification
  expect_true(all(v$ok))
  expect_true(all(v$residual < 1e-8))
  expect_equal(unique(sort(c(paste(v$n_minus, v$status_minus), paste(v$n_plus, v$status_plus)))),
               c("1 edge", "1 stable"))
  expect_equal(br$evaluations, sum(br$verification$evaluations) +
                 br$evaluations - sum(br$verification$evaluations))
  expect_gt(br$evaluations, 0L)
  # cubic Hermite through the points and tangents stays on the circle
  dense <- as.data.frame(br, n = 8)
  expect_equal(nrow(dense), 8L * (nrow(pts) - 1L) + 1L)
  expect_lt(max(abs(sqrt(dense[["p[1]"]]^2 + dense[["p[2]"]]^2) - 1)), 1e-6)
  out <- paste(utils::capture.output(print(br)), collapse = "\n")
  expect_match(out, "branching boundary of 1 resident in p\\[1\\] x p\\[2\\]")
  expect_match(out, "forward: closed")
})

test_that("GK98's 2 -> 3 boundary is followed as c2(kappa) to where its pair merges, and switches to the 1 -> 2 locus", {
  box <- rbind(c(1, 3), c(0.5, 4.5))
  br <- community_continue_boundary(gk98_pair_start(), gk98_by_d_kappa, y = 0.05, p = c(1.93, 1),
                                    box = box, control = tight)
  pts <- br$points
  expect_equal(unname(br$termination), c("event: separation", "left the box"))
  # on the curve to the corrector's tolerance (1e-8 on every residual) over
  # the slopes of the residual, order 0.01-0.1 while kappa < 3.5; towards the
  # end they fall with the pair's separation (a* ~ sqrt(4 - kappa)), and the
  # accuracy with them
  away <- pts[["p[2]"]] < 3.5
  c2 <- vapply(pts[["p[2]"]], function(k) gk98_kappa_boundaries(k)[["c2"]], 0)
  expect_lt(max(abs(pts[["p[1]"]] - c2)[away]), 1e-6)
  expect_lt(max(abs(pts[["p[1]"]] - c2)), 1e-4)
  a <- mapply(gk98_kappa_pair_root, c2, pts[["p[2]"]])
  expect_lt(max(abs(pts[["x[2]"]] - a)[away]), 1e-6)
  expect_lt(max(abs(pts[["x[2]"]] - a)), 1e-3)
  expect_lt(max(abs(pts[["y:x"]])), 1e-4)
  expect_true(all(pts[["birth_rate[1]"]] > 0))
  # the separation test fires where the pair is merge_tol of the bounds'
  # width (0.08) either side of y* = 0; the pair's error there (~3e-4) over
  # da*/dkappa (~1.4) places it to ~2e-4
  expect_equal(br$events$type, "separation")
  expect_lt(abs(br$events[["p[2]"]] - gk98_merge_kappa(0.08)), 5e-4)
  expect_true(all(br$verification$ok))

  # the merged pair is one resident at 0, a branching point on the 1 -> 2
  # locus, which continues to where its pitchfork turns subcritical: there
  # the landscape beyond the resident's neighbourhood (2 merge_tol of the
  # width, 0.16) rises past fitness_tol, and the locus stops being a
  # coexistence boundary
  sw <- boundary_switch(br, gk98_by_d_kappa)
  expect_equal(sw$type, "branching")
  expect_equal(sw$which, 1L)
  expect_equal(nrow(sw$community$traits), 1L)
  expect_lt(abs(sw$community$traits[1, 1]), 1e-6)
  locus <- do.call(community_continue_boundary,
                   c(sw, list(box = rbind(c(1, 3), c(0.5, 6)), control = tight)))
  lp <- locus$points
  expect_equal(locus$type, "branching")
  expect_lt(max(abs(lp[["p[1]"]] - sqrt((2 + lp[["p[2]"]]) / 2))), 1e-8)
  expect_equal(locus$events$type, "other_peak")
  expect_equal(locus$events[["p[2]"]], gk98_subcritical_kappa(0.16, 1e-12), tolerance = 1e-6)
  expect_equal(locus$termination[["backward"]], "left the box")
  expect_true(all(locus$verification$ok))

  d <- boundary_diagram(br, locus)
  expect_s3_class(d, "boundary_diagram")
  df <- as.data.frame(d, n = 4)
  expect_equal(sort(unique(df$type)), c("branching", "invasion"))
  expect_s3_class(plot(d), "ggplot")
  expect_s3_class(plot(br), "ggplot")
})

test_that("an event beyond the box is not on the branch, which leaves the box first", {
  # the separation event is at kappa = 3.96696; a box ending just short of
  # it is left on the step that would cross it
  br <- community_continue_boundary(gk98_pair_start(), gk98_by_d_kappa, y = 0.05, p = c(1.93, 1),
                                    box = rbind(c(1, 3), c(0.5, 3.9669)), direction = "forward",
                                    control = tight)
  expect_equal(unname(br$termination), "left the box")
  expect_equal(nrow(br$events), 0L)
  expect_lte(max(br$points[["p[2]"]]), 3.9669)
})

test_that("the other-peak search keeps every neighbourhood's edges", {
  # a landscape whose highest point away from the resident is exactly on the
  # edge of its neighbourhood (where a pitchfork turning subcritical first
  # shows): z +- reach can round to a hair inside reach, and the edge must
  # still be searched
  ctrl <- boundary_control()
  tf <- community_trait_transform(list(trait_scale = "linear"))
  width <- 8
  reach <- 2 * ctrl$merge_tol * width
  z <- seq(-3, 3, length.out = 401)
  fitness <- vapply(seq_along(z), function(i) {
    peak <- z[i] + if (i %% 2) reach else -reach
    comm <- list(bounds = bounds(x = c(-4, 4)), trait_names = "x",
                 fitness_function = function(x) -(as.numeric(x) - peak)^2)
    boundary_other_peak(comm, matrix(z[i]), width, tf, ctrl)$fitness
  }, 0)
  expect_equal(fitness, rep(0, length(z)))
})

test_that("the rectangular Jacobian's columns match a finer difference of the residual", {
  # both free parameters, the GK98 invasion system away from the boundary:
  # the layer's central differences at 1e-3 relative, good to ~1e-6, against
  # 1e-5 on the whole residual
  a <- gk98_kappa_pair_root(1.8, 1.3)
  pair <- gk98_by_d_kappa(c(1.8, 1.3)) |> community_clear_residents() |>
    community_add(trait_matrix(c(-a, a + 0.05), "x"), birth_rate = c(1.6, 1.6))
  p <- c(`p[1]` = 1.8, `p[2]` = 1.3)
  residual <- boundary_residual_fn(pair, gk98_by_d_kappa, p, 1:2, "invasion", 1L)
  u <- c(-a, a + 0.05, 0.1, 1.8, 1.3)
  J <- community_boundary_jacobian(residual, u)
  expect_equal(dim(J), c(4L, 5L))
  expect_equal(colnames(J), c("x[1]", "x[2]", "y:x", "p[1]", "p[2]"))
  expect_lt(max(abs(unname(J) - oracle_jacobian(residual, u, h = 1e-5))), 2e-6)
  # 2mk for the residents, two per free parameter
  expect_identical(attr(J, "evaluations"), 8L)
  # the tangent is its null vector, of unit length
  t <- boundary_tangent(unname(J))
  expect_equal(sqrt(sum(t^2)), 1)
  expect_lt(max(abs(unname(J) %*% t)), 1e-12)
})

test_that("a turning tangent rejects the step, which is halved", {
  # a tangent tolerance near one rejects any step the circle's curvature
  # turns by more than ~0.0045 rad: steps halve to fit, and the curve is
  # still followed
  ctl <- boundary_control(list(tangent_tol = 0.99999, ds = 0.05, max_steps = 6))
  br <- community_continue_boundary(dd99_circle(c(0.9, 0.1)), dd99_circle, type = "branching",
                                    p = c(0.9, 0.1), direction = "forward", control = ctl)
  expect_gt(br$rejected, 0L)
  expect_equal(unname(br$termination), "max_steps")
  expect_lt(max(abs(sqrt(br$points[["p[1]"]]^2 + br$points[["p[2]"]]^2) - 1)), 1e-8)
  expect_true(all(diff(br$points$s) < 0.0045 + 1e-9))
})

test_that("a branch is the same, in answer and cost, under each plan", {
  ctl <- boundary_control(list(max_steps = 8))
  sequential <- community_continue_boundary(dd99_circle(c(0.9, 0.1)), dd99_circle, type = "branching",
                                            p = c(0.9, 0.1), control = ctl)
  for (type in c("multicore", "multisession")) {
    parallel <- with_two_workers(type, community_continue_boundary(
      dd99_circle(c(0.9, 0.1)), dd99_circle, type = "branching", p = c(0.9, 0.1), control = ctl))
    expect_identical(parallel$points, sequential$points)
    expect_identical(parallel$u, sequential$u)
    expect_identical(parallel$verification, sequential$verification)
    expect_identical(parallel$evaluations, sequential$evaluations)
  }
})

test_that("a verified point whose labels do not change across the curve is warned of", {
  # labels that agree on both sides: the branch is no boundary there
  local_mocked_bindings(community_label = function(community, control) {
    structure(community, label = list(n = 1L, status = "stable", reason = NULL, rounds = 1L, evaluations = 1L))
  })
  ctl <- boundary_control(list(max_steps = 2))
  expect_warning(br <- community_continue_boundary(dd99_circle(c(0.9, 0.1)), dd99_circle, type = "branching",
                                                   p = c(0.9, 0.1), direction = "forward", control = ctl),
                 "verified points failed")
  expect_false(any(br$verification$ok))
})

test_that("boundary_switch seeds the boundary each kind of event starts", {
  br <- community_continue_boundary(gk98_pair_start(), gk98_by_d_kappa, y = 0.05, p = c(1.93, 1),
                                    direction = "forward",
                                    control = boundary_control(list(max_steps = 2)))
  expect_error(boundary_switch(br, gk98_by_d_kappa), "not an event \\(max_steps\\)")
  expect_error(boundary_switch(br, gk98_by_d_kappa, end = "backward"), "not an event")
  with_event <- function(type, detail) {
    br$events <- data.frame(end = "forward", type = type, s = 0, `p[1]` = 1.93, `p[2]` = 1, check.names = FALSE)
    br$events$detail <- I(list(detail))
    br
  }
  x <- br$u[nrow(br$u), 1:2]
  # a resident going extinct: the community without it, its trait marginal
  sw <- boundary_switch(with_event("density", list(density = 2L)), gk98_by_d_kappa)
  expect_equal(sw$type, "invasion")
  expect_equal(as.numeric(sw$community$traits), x[[1]])
  expect_equal(sw$y, x[[2]])
  # a resident ceasing to be an ESS: its branching locus, the same residents
  sw <- boundary_switch(with_event("resident_curvature", list(resident_curvature = 2L)), gk98_by_d_kappa)
  expect_equal(c(sw$type, sw$which), c("branching", "2"))
  expect_equal(as.numeric(sw$community$traits), unname(x))
  # another peak becoming marginal: the same residents, that peak's trait
  sw <- boundary_switch(with_event("other_peak", list(other_peak = 0.7)), gk98_by_d_kappa)
  expect_equal(c(sw$type, sw$y), c("invasion", "0.7"))
  expect_equal(sw$free, 1:2)
  expect_equal(sw$p, br$u[nrow(br$u), 4:5], ignore_attr = TRUE)
  # what starts no boundary of these kinds
  expect_error(boundary_switch(with_event("bordered", list()), gk98_by_d_kappa), "starts no boundary")
  expect_error(boundary_switch(with_event("peak_curvature", list()), gk98_by_d_kappa), "starts no boundary")
  expect_error(boundary_switch(list(), gk98_by_d_kappa), "needs a regnans_branch")
})

test_that("community_continue_boundary validates its arguments", {
  start <- dd99_circle(c(0.9, 0.1))
  expect_error(community_continue_boundary(start, dd99_circle, type = "branching", p = c(0.9, 0.1), free = 1),
               "two distinct parameters")
  expect_error(community_continue_boundary(start, dd99_circle, type = "branching", p = c(0.9, 0.1),
                                           box = rbind(c(1, 0), c(0, 1))), "lower and upper limits")
  expect_error(community_continue_boundary(start, dd99_circle, type = "branching", p = c(0.9, 0.1),
                                           direction = "sideways"), "should be one of")
  expect_error(boundary_diagram(list()), "regnans_branch")
})
