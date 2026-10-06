# Pairwise invasibility, its zero contours, mutual invasibility and the
# trait-evolution plot.
#
# DD99 oracle for a single resident x at N = K(x):
#   s(y; x) = r (1 - K(x) C(y - x) / K(y))
# The zero set is y = x together with the curve where
#   (y - x0)^2 / sigma_K^2 - (x - x0)^2 / sigma_K^2 = (y - x)^2 / sigma_C^2,
# a quadratic in y whose second root is
#   y2(x) = (x (sigma_K^2 + sigma_C^2) - 2 x0 sigma_C^2) / (sigma_K^2 - sigma_C^2)
# (for sigma_C != sigma_K), so every contour point has a closed form.

dd99_s <- function(y, x, x0 = 0, sigma_K = 1, sigma_C = 0.4, r = 1) {
  K <- function(z) exp(-(z - x0)^2 / (2 * sigma_K^2))
  r * (1 - K(x) * exp(-(y - x)^2 / (2 * sigma_C^2)) / K(y))
}
dd99_y2 <- function(x, x0 = 0, sigma_K = 1, sigma_C = 0.4) {
  (x * (sigma_K^2 + sigma_C^2) - 2 * x0 * sigma_C^2) / (sigma_K^2 - sigma_C^2)
}

dd99 <- function(sigma_C = 0.4, x0 = 0, b = c(-2, 2)) {
  community_start(bounds(x = b), trait_scale = "linear",
                  harness = harness_dd99(x0 = x0, sigma_K = 1, sigma_C = sigma_C))
}
ctrl <- function(...) pip_control(list(n_resident = 11, n_mutant = 61, refine = 0, ...))

built <- function(p) {
  expect_s3_class(p, "ggplot")
  ggplot2::ggplot_build(p)
}

test_that("pip_control validates its settings", {
  expect_equal(pip_control()$n_resident, 41L)
  expect_equal(pip_control()$seed, "interpolate")
  expect_identical(pip_control(list(refine = 3))$refine, 3L)
  expect_equal(pip_control(list(seed = "cold"))$seed, "cold")
  expect_error(pip_control(list(seed = "warm")), "should be one of")
  expect_error(pip_control(list(n = 5)), "Unknown control parameters n")
  expect_error(pip_control(list(n_resident = 1)), "at least 2")
  expect_error(pip_control(list(refine = 1.5)), "whole number")
  expect_error(pip_control(list(tol = 0)), "positive number")
  expect_error(pip_control(list(parallel = NA)), "TRUE or FALSE")
})

test_that("community_pip reproduces the DD99 surface and solves every resident", {
  pip <- community_pip(dd99(), control = ctrl())
  expect_s3_class(pip, "pip")
  expect_equal(names(pip$surface), c("resident", "mutant", "fitness"))
  expect_equal(nrow(pip$surface), 11L * 61L)
  expect_equal(range(pip$surface$resident), c(-2, 2))
  expect_equal(pip$trait_name, "x")
  expect_equal(pip$trait_scale, "linear")
  expect_equal(pip$surface$fitness, dd99_s(pip$surface$mutant, pip$surface$resident),
               tolerance = 1e-10)
  expect_equal(pip$residents$birth_rate, 500 * exp(-pip$residents$resident^2 / 2),
               tolerance = 1e-8)
  expect_equal(as.data.frame(pip), as.data.frame(pip$surface))
  expect_equal(names(pip$residents), c("resident", "birth_rate", "n_evals", "n_crossings"))
  expect_true(all(pip$residents$n_evals >= 1L))
  expect_true(is.numeric(pip$elapsed) && pip$elapsed >= 0)
  expect_match(paste(utils::capture.output(print(pip)), collapse = ""),
               "11 residents x 61 mutants")
})

test_that("seeding strategies give the same surface, and interpolated seeds cost less on an iterated model", {
  comm <- community_start(bounds(x = c(0.08, 0.9)), trait_scale = "log",
                          harness = harness_iterate_demography(harness_gm99(alpha = 7, beta = 15)))
  ctl <- function(seed) pip_control(list(n_resident = 13, n_mutant = 41, n_coarse = 4,
                                         refine = 0, seed = seed))
  # the fixed-point map is slowly convergent near the viability edge (its
  # multiplier approaches one), so a tight eps is needed for the three runs to
  # agree to the surface's precision rather than to eps / (1 - multiplier)
  comm$demography_control$equilibrium_eps <- 1e-9
  comm$demography_control$equilibrium_nsteps <- 1000
  cold <- community_pip(comm, control = ctl("cold"))
  neighbour <- community_pip(comm, control = ctl("neighbour"))
  interpolate <- community_pip(comm, control = ctl("interpolate"))
  expect_equal(neighbour$surface$fitness, cold$surface$fitness, tolerance = 1e-6)
  expect_equal(interpolate$surface$fitness, cold$surface$fitness, tolerance = 1e-6)
  expect_equal(interpolate$residents$birth_rate, cold$residents$birth_rate, tolerance = 1e-6)
  expect_lt(sum(interpolate$residents$n_evals), sum(cold$residents$n_evals))
  expect_lt(sum(neighbour$residents$n_evals), sum(cold$residents$n_evals))
  # the four coarse residents start cold; the rest start near their answer
  coarse <- unique(round(seq(1, 13, length.out = 4)))
  expect_lt(mean(interpolate$residents$n_evals[-coarse]), mean(interpolate$residents$n_evals[coarse]))
})

test_that("the zero contours are located exactly, diagonal included", {
  pip <- community_pip(dd99(), control = ctrl())
  co <- pip$contours
  # every resident has its diagonal crossing
  diag <- co[abs(co$resident - co$mutant) < 1e-10, ]
  expect_equal(sort(diag$resident), pip$residents$resident)
  # and, where the second contour lies inside the mutant range, that one too
  off <- co[abs(co$resident - co$mutant) >= 1e-10, ]
  expect_gt(nrow(off), 0L)
  expect_equal(off$mutant, dd99_y2(off$resident), tolerance = 1e-7)
  expected_in_range <- sum(abs(dd99_y2(pip$residents$resident)) <= 2 &
                           abs(dd99_y2(pip$residents$resident) - pip$residents$resident) > 1e-8)
  expect_equal(nrow(off), expected_in_range)
  # fitness really is zero there
  s <- dd99_s(off$mutant, off$resident)
  expect_lt(max(abs(s)), 1e-7)
})

test_that("the near-diagonal contour is resolved at a resident next to the singular strategy", {
  # at x = 0.01 the second root is at y2 = x (1 + 0.16)/(1 - 0.16) ~ 0.0138, which
  # a 61-point mutant grid over [-2, 2] (cell 0.067) cannot separate from x
  pip <- community_pip(dd99(), resident = c(0.01), control = ctrl())
  off <- pip$contours[abs(pip$contours$resident - pip$contours$mutant) > 1e-10, ]
  expect_equal(nrow(off), 1L)
  expect_equal(off$mutant, dd99_y2(0.01), tolerance = 1e-7)
})

test_that("refinement adds residents only where the contours disagree", {
  coarse <- community_pip(dd99(), control = pip_control(list(n_resident = 7, n_mutant = 61, refine = 0)))
  refined <- community_pip(dd99(), control = pip_control(list(n_resident = 7, n_mutant = 61, refine = 2)))
  expect_gt(nrow(refined$residents), nrow(coarse$residents))
  expect_true(all(coarse$residents$resident %in% refined$residents$resident))
  # the second contour is a straight line, so no refinement for curvature; it
  # leaves the mutant range at |y2| = 2 (|x| = 2 * 0.84 / 1.16 ~ 1.45) and meets
  # the diagonal at the singular strategy x = 0, and the new residents sit at
  # those three places only
  added <- setdiff(refined$residents$resident, coarse$residents$resident)
  near_exit <- abs(abs(added) - 1.45) < 0.67
  near_singular <- abs(added) < 0.67
  expect_true(all(near_exit | near_singular))
  expect_true(any(near_exit) && any(near_singular))
  # refined midpoints were seeded from their neighbours and still converge to K(x)
  expect_equal(refined$residents$birth_rate, 500 * exp(-refined$residents$resident^2 / 2),
               tolerance = 1e-8)
})

test_that("community_pip takes explicit grids, ignores residents, and honours a log scale", {
  comm <- dd99() |> community_add(trait_matrix(1.5, "x"), birth_rate = 100)
  pip <- community_pip(comm, resident = c(-1, 0, 1), mutant = seq(-0.5, 0.5, by = 0.5),
                       control = ctrl())
  expect_equal(nrow(pip$surface), 9L)
  expect_equal(pip$residents$resident, c(-1, 0, 1))

  h <- harness_dd99(x0 = 0.5, sigma_K = 1, sigma_C = 0.4, trait_name = "lma")
  lp <- community_start(bounds(lma = c(0.05, 2)), trait_scale = "log", harness = h) |>
    community_pip(control = pip_control(list(n_resident = 9, n_mutant = 33, refine = 0)))
  res <- lp$residents$resident
  expect_equal(diff(log(res)), rep(diff(log(c(0.05, 2))) / 8, 8), tolerance = 1e-10)
  expect_equal(lp$trait_scale, "log")
  off <- lp$contours[abs(lp$contours$resident - lp$contours$mutant) > 1e-10, ]
  expect_equal(off$mutant, dd99_y2(off$resident, x0 = 0.5), tolerance = 1e-7)
})

test_that("community_pip refuses multi-trait communities and infinite bounds", {
  two <- community_start(bounds(x1 = c(-2, 2), x2 = c(-2, 2)), trait_scale = "linear",
                         harness = harness_dd99_nd())
  expect_error(community_pip(two), "single-trait community")
  inf <- community_start(bounds_infinite("x"), trait_scale = "linear", harness = harness_dd99())
  expect_error(community_pip(inf), "finite bounds")
})

test_that("mutual invasibility compares the surface with its transpose", {
  pip <- community_pip(dd99(), control = ctrl())
  m <- pip_mutual(pip)
  x <- pip$residents$resident
  expect_equal(nrow(m), length(x)^2)
  expect_equal(m$s12, dd99_s(m$x2, m$x1), tolerance = 1e-6)
  expect_equal(m$s21, dd99_s(m$x1, m$x2), tolerance = 1e-6)
  expect_equal(m$mutual, m$x1 != m$x2 & dd99_s(m$x2, m$x1) > 0 & dd99_s(m$x1, m$x2) > 0)
  expect_true(any(m$mutual))          # a narrow kernel leaves room for dimorphisms
})

test_that("the trait-evolution plot solves each dimorphism and reports its gradients", {
  comm <- dd99()
  pip <- community_pip(comm, control = ctrl())
  tep <- community_tep(comm, pip)
  expect_s3_class(tep, "tep")
  expect_equal(names(tep), c("x1", "x2", "n1", "n2", "g1", "g2"))
  expect_true(all(tep$x1 < tep$x2))
  expect_true(all(tep$n1 > 0 & tep$n2 > 0))
  m <- pip_mutual(pip)
  expect_equal(nrow(tep) + attr(tep, "excluded"), sum(m$mutual & m$x1 < m$x2))
  # DD99 is symmetric about x0 = 0: a mirror-image pair feels mirror-image selection
  sym <- tep[abs(tep$x1 + tep$x2) < 1e-8, ]
  expect_gt(nrow(sym), 0L)
  expect_equal(sym$g1, -sym$g2, tolerance = 1e-8)
  expect_equal(sym$n1, sym$n2, tolerance = 1e-8)

  other <- community_start(bounds(y = c(-2, 2)), trait_scale = "linear",
                           harness = harness_dd99(trait_name = "y"))
  expect_error(community_tep(other, pip), "computed for trait x")
})

test_that("the raw-cell plots build, on linear and log axes", {
  pip <- community_pip(dd99(), control = ctrl())
  b <- built(plot(pip, smooth = FALSE))
  expect_equal(length(b$data), 3L)                 # cells, contour points, diagonal
  expect_equal(nrow(b$data[[1]]), 11L * 61L)
  expect_equal(nrow(b$data[[2]]), nrow(pip$contours))
  expect_equal(length(built(plot(pip, smooth = FALSE, contours = FALSE))$data), 2L)
  expect_equal(length(built(plot(pip, smooth = FALSE, fill = "fitness"))$data), 3L)
  b <- built(plot(pip, type = "mip", smooth = FALSE))
  expect_equal(nrow(b$data[[1]]), 11L * 11L)

  tep <- community_tep(dd99(), pip)
  expect_s3_class(attr(tep, "pip"), "pip")
  b <- built(plot(tep, n_arrows = 6, n_display = 30))
  expect_equal(length(b$data), 4L)                 # coexistence region, arrows, lattice points, diagonal
  expect_lte(nrow(b$data[[2]]), 6L * 6L)
  # mirrored: arrows on both sides of the diagonal (cells straddling it keep one)
  above <- sum(b$data[[2]]$x < b$data[[2]]$y)
  below <- sum(b$data[[2]]$x > b$data[[2]]$y)
  expect_gt(above, 0L)
  expect_gt(below, 0L)
  expect_lte(abs(above - below), 6L)
  expect_equal(nrow(b$data[[2]]), nrow(b$data[[3]]))

  lp <- community_start(bounds(lma = c(0.05, 2)), trait_scale = "log",
                        harness = harness_dd99(x0 = 0.5, trait_name = "lma")) |>
    community_pip(control = pip_control(list(n_resident = 7, n_mutant = 21, refine = 0)))
  expect_equal(nrow(built(plot(lp, smooth = FALSE))$data[[1]]), 7L * 21L)
})

test_that("GM99 non-viable mutants are shaded as unable to invade", {
  pip <- community_start(bounds(x = c(0.02, 0.9)), trait_scale = "log",
                         harness = harness_gm99(alpha = 7, beta = 15)) |>
    community_pip(control = pip_control(list(n_resident = 5, n_mutant = 21, refine = 0)))
  expect_true(any(is.infinite(pip$surface$fitness)))
  expect_equal(nrow(built(plot(pip, smooth = FALSE, fill = "fitness"))$data[[1]]), 5L * 21L)
  expect_true(all(is.finite(pip$contours$mutant)))
})

test_that("residents are solved in contiguous chunks and a parallel plan gives the same surface", {
  expect_equal(regnans_chunks(0), list())
  expect_equal(regnans_chunks(5, parallel = FALSE), list(1:5))
  expect_equal(regnans_map(1:3, function(i) i * 2, parallel = FALSE), list(2, 4, 6))

  skip_if_not_installed("future")
  skip_if_not_installed("future.apply")
  skip_on_os("windows")
  skip_if_not(future::supportsMulticore())
  sequential <- community_pip(dd99(), control = ctrl())
  old <- future::plan(future::multicore, workers = 2)
  on.exit(future::plan(old), add = TRUE)
  expect_equal(regnans_workers(), 2L)
  expect_equal(sort(lengths(regnans_chunks(11))), c(5L, 6L))
  expect_equal(unlist(regnans_chunks(11)), 1:11)
  parallel <- community_pip(dd99(), control = ctrl())
  expect_equal(parallel$surface, sequential$surface, tolerance = 1e-10)
  expect_equal(parallel$contours, sequential$contours, tolerance = 1e-8)
})

# ---- other models ------------------------------------------------------------
#
# GK98, symmetric three-patch (mu = (-d, 0, d), equal K), single resident x:
#   S(y; x) = log( (1/3) sum_j f_j(y) / f_j(x) )
# With u = y - x and s = y + x, (y - mu)^2 - (x - mu)^2 = u (s - 2 mu), so the
# zero set besides the diagonal is the single curve
#   s = G(u) = (2 sigma^2 / u) log( (1 + 2 cosh(u d / sigma^2)) / 3 ),
# through the singular strategy x* = 0. G is odd and increasing with
# G'(0) = 2 d^2 / (3 sigma^2) and G -> 2d as u grows, so a resident's column
# x = const meets the curve once when G'(0) < 1 and up to three times (the
# curve folds) when G'(0) > 1 -- which is d/sigma > sqrt(3/2), the branching
# condition.
#
#
# JJ12, single resident x at n*:
#   w(y; x) = (1 - p) R(y) C(y) / (R(x) C(x)) + p
# so s(y; x) > 0 iff R(y)C(y) > R(x)C(x): a strict ordering of strategies, with
# the second contour at y = 2 x* - x and no mutual invasibility anywhere.

gk98_s <- function(y, x, d, sigma = 1) {
  f <- function(z, mu) exp(-(z - mu)^2 / (2 * sigma^2))
  log(rowMeans(sapply(c(-d, 0, d), function(mu) f(y, mu) / f(x, mu))))
}
gk98_contour_sum <- function(u, d, sigma = 1) {
  (2 * sigma^2 / u) * log((1 + 2 * cosh(u * d / sigma^2)) / 3)
}

test_that("GK98 contours lie on the closed-form curve through the singular strategy", {
  for (d in c(1.0, 1.5)) {
    comm <- community_start(bounds(x = c(-3, 3)), trait_scale = "linear",
                            harness = harness_gk98(d = d, sigma = 1))
    pip <- community_pip(comm, control = pip_control(list(n_resident = 13, n_mutant = 121, refine = 0)))
    expect_equal(pip$surface$fitness, gk98_s(pip$surface$mutant, pip$surface$resident, d),
                 tolerance = 1e-10, info = paste("d =", d))
    co <- pip$contours
    expect_equal(sort(co$resident[abs(co$resident - co$mutant) < 1e-8]), pip$residents$resident,
                 info = paste("d =", d))
    off <- co[abs(co$resident - co$mutant) >= 1e-8, ]
    expect_gt(nrow(off), 0L)
    u <- off$mutant - off$resident
    expect_equal(off$mutant + off$resident, gk98_contour_sum(u, d), tolerance = 1e-7,
                 info = paste("d =", d))
    expect_true(all(pip$residents$n_crossings >= 1L), info = paste("d =", d))
    if (d < sqrt(1.5)) {
      expect_true(all(pip$residents$n_crossings <= 2L), info = paste("d =", d))
    } else {
      expect_true(any(pip$residents$n_crossings == 3L), info = paste("d =", d))
    }
  }
})

test_that("GK98 mutual invasibility and trait evolution are mirror-symmetric, with a dimorphic coalition when branching", {
  comm <- community_start(bounds(x = c(-3, 3)), trait_scale = "linear",
                          harness = harness_gk98(d = 1.5, sigma = 1))
  pip <- community_pip(comm, control = pip_control(list(n_resident = 25, n_mutant = 121, refine = 0)))
  m <- pip_mutual(pip)
  key <- function(a, b) paste(signif(a, 10), signif(b, 10))
  mm <- stats::setNames(m$mutual, key(m$x1, m$x2))
  expect_equal(unname(mm[key(m$x2, m$x1)]), unname(mm))      # s12 & s21 symmetric
  expect_equal(unname(mm[key(-m$x2, -m$x1)]), unname(mm))    # mirror symmetric
  expect_true(any(m$mutual))

  tep <- community_tep(comm, pip)
  expect_gt(nrow(tep), 0L)
  sym <- tep[abs(tep$x1 + tep$x2) < 1e-8, ]
  expect_gt(nrow(sym), 2L)
  expect_equal(sym$g1, -sym$g2, tolerance = 1e-8)
  expect_equal(sym$n1, sym$n2, tolerance = 1e-8)
  # along the symmetric line the gradient on the outer resident changes sign:
  # the dimorphism is carried to a coalition at (-a*, a*)
  expect_true(any(sym$g2 > 0) && any(sym$g2 < 0))
})

test_that("JJ12 is a CSS: the second contour is y = 2 x* - x and nothing is mutually invasible", {
  a <- 0.1; sigma <- 1.4; x_opt <- 0.5
  x_star <- x_opt - a * sigma^2
  comm <- community_start(bounds(x = c(-3, 3)), trait_scale = "linear",
                          harness = harness_jj12(a = a, x_opt = x_opt, sigma = sigma))
  pip <- community_pip(comm, control = pip_control(list(n_resident = 13, n_mutant = 121, refine = 1)))
  off <- pip$contours[abs(pip$contours$resident - pip$contours$mutant) > 1e-8, ]
  expect_gt(nrow(off), 0L)
  expect_equal(off$mutant, 2 * x_star - off$resident, tolerance = 1e-7)
  # refinement found the singular strategy, where the two contours meet
  expect_lt(min(abs(pip$residents$resident - x_star)), 6 / 12 / 2 + 1e-8)
  expect_false(any(pip_mutual(pip)$mutual))
  tep <- community_tep(comm, pip)
  expect_equal(nrow(tep), 0L)
  expect_equal(attr(tep, "excluded"), 0L)
  # the empty views say so rather than drawing a blank
  b <- built(plot(pip, type = "mip", n_display = 20))
  expect_true(any(vapply(b$data, function(d) "label" %in% names(d) && any(grepl("no mutually", d$label)), logical(1))))
  b <- built(plot(tep))
  expect_true(any(vapply(b$data, function(d) "label" %in% names(d) && any(grepl("no coexisting", d$label)), logical(1))))
})

test_that("an iterated GK98 gives the same surface and contours as the closed form", {
  h <- harness_gk98(d = 1.5, sigma = 1)
  ctl <- pip_control(list(n_resident = 9, n_mutant = 61, refine = 1))
  closed <- community_pip(community_start(bounds(x = c(-3, 3)), trait_scale = "linear", harness = h),
                          control = ctl)
  numerical <- community_start(bounds(x = c(-3, 3)), trait_scale = "linear",
                               harness = harness_iterate_demography(h))
  numerical$demography_control$equilibrium_eps <- 1e-10
  iterated <- community_pip(numerical, control = ctl)
  expect_equal(iterated$residents$resident, closed$residents$resident)
  expect_equal(iterated$surface$fitness, closed$surface$fitness, tolerance = 1e-6)
  expect_equal(iterated$contours$mutant, closed$contours$mutant, tolerance = 1e-6)
  expect_true(all(iterated$residents$n_evals >= closed$residents$n_evals))
})

# ---- smooth boundaries from the branches -------------------------------------

test_that("crossings link into branches and the derived sign matches the DD99 oracle away from the contours", {
  pip <- community_pip(dd99(), control = pip_control(list(n_resident = 13, n_mutant = 61, refine = 2)))
  br <- pip_branches(pip)
  expect_equal(names(br), c("branch", "resident", "mutant", "z_resident", "z_mutant"))
  # every crossing belongs to a branch; the one at the singular strategy to two
  expect_gte(nrow(br), nrow(pip$contours))
  # the diagonal is one branch spanning every resident, the second contour one
  # more, both passing through the intersection at x* = 0
  on_diag <- abs(br$resident - br$mutant) < 1e-8
  diag_branch <- names(which.max(table(br$branch[on_diag])))
  expect_equal(sum(br$branch == as.integer(diag_branch)), nrow(pip$residents))
  expect_equal(length(unique(br$branch)), 2L)
  other <- br[br$branch != as.integer(diag_branch), ]
  expect_equal(other$mutant, dd99_y2(other$resident), tolerance = 1e-7)

  sgn <- pip_sign_function(pip)
  g <- expand.grid(x = seq(-1.9, 1.9, length.out = 41), y = seq(-1.9, 1.9, length.out = 41))
  truth <- dd99_s(g$y, g$x)
  away <- abs(g$y - g$x) > 0.08 & abs(g$y - dd99_y2(g$x)) > 0.08
  expect_equal(sgn(g$x[away], g$y[away]), sign(truth[away]))
})

test_that("the smooth views build as ribbons and agree with the raw cells on the residents", {
  pip <- community_pip(dd99(), control = pip_control(list(n_resident = 9, n_mutant = 41, refine = 1)))
  regions <- pip_regions(pip, 40)
  expect_equal(names(regions), c("region", "x", "ymin", "ymax", "value"))
  expect_true(all(regions$ymax >= regions$ymin - 1e-12))
  # the ribbons tile each column: the lowest starts at the mutant floor, the
  # highest ends at the ceiling, and the sign alternates up the column
  col <- regions[abs(regions$x - regions$x[which.min(abs(regions$x - 0.9))]) < 1e-12, ]
  col <- col[order(col$ymin), ]
  expect_equal(min(col$ymin), -2)
  expect_equal(max(col$ymax), 2)
  expect_true(all(col$ymin[-1] == col$ymax[-nrow(col)]))
  expect_true(all(col$value[-1] != col$value[-nrow(col)]))
  b <- built(plot(pip, n_display = 40))
  expect_equal(length(b$data), 3L)                      # ribbons, branch lines, diagonal
  expect_gt(nrow(b$data[[2]]), 0L)
  b <- built(plot(pip, type = "mip", n_display = 30))
  expect_gt(nrow(b$data[[1]]), 0L)
  m <- pip_regions(pip, 30, "mip")
  expect_true(any(m$value) && any(!m$value))
  # a region may taper to nothing where two contours cross, but never inverts
  # and never has zero width throughout
  expect_true(all(m$ymax - m$ymin >= -1e-12))
  expect_true(all(tapply(m$ymax - m$ymin, m$region, max) > 0))
  # raw cells are still available
  expect_equal(nrow(built(plot(pip, smooth = FALSE))$data[[1]]), nrow(pip$surface))

  sgn <- pip_sign_function(pip)
  at <- pip$surface[abs(pip$surface$fitness) > 1e-3, ]
  expect_equal(sgn(at$resident, at$mutant), sign(at$fitness))
})

test_that("GK98's folded contour links into branches that the smooth view can draw", {
  pip <- community_pip(community_start(bounds(x = c(-3, 3)), trait_scale = "linear",
                                       harness = harness_gk98(d = 1.5)),
                       control = pip_control(list(n_resident = 25, n_mutant = 121, refine = 2)))
  br <- pip_branches(pip)
  expect_gte(nrow(br), nrow(pip$contours))
  expect_gte(length(unique(br$branch)), 2L)
  expect_equal(length(built(plot(pip, n_display = 60))$data), 3L)
})

test_that("the plots are ggplot objects that can be extended, and autoplot is the same", {
  pip <- community_pip(dd99(), control = pip_control(list(n_resident = 7, n_mutant = 31, refine = 0)))
  p <- plot(pip, n_display = 20)
  expect_s3_class(p, "ggplot")
  q <- p + ggplot2::labs(x = "resident trait", title = "DD99") + ggplot2::theme_minimal()
  b <- ggplot2::ggplot_build(q)
  expect_equal(b$plot$labels$x, "resident trait")
  expect_s3_class(ggplot2::autoplot(pip, n_display = 20), "ggplot")
  tep <- community_tep(dd99(), pip)
  expect_s3_class(ggplot2::autoplot(tep), "ggplot")
})

test_that("a fold's arms end together and are joined at an estimated vertex", {
  pip <- community_pip(community_start(bounds(x = c(-3, 3)), trait_scale = "linear",
                                       harness = harness_gk98(d = 1.5)),
                       control = pip_control(list(n_resident = 31, n_mutant = 201, refine = 2)))
  br <- pip_branches(pip)
  # the diagonal is one branch over every resident and never shares a point
  # with another branch except at the singular strategy
  on_diag <- abs(br$resident - br$mutant) < 1e-8
  diag_id <- as.integer(names(which.max(table(br$branch[on_diag]))))
  expect_equal(sum(br$branch == diag_id), nrow(pip$residents))
  shared <- br[br$branch != diag_id & on_diag, ]
  expect_true(all(abs(shared$resident) < 1e-6))
  # the two arms of the fold on x > 0 end at the same resident and get one
  # vertex each, at the same place, beyond their last resident
  arms <- br[br$branch != diag_id & br$resident > 0 & br$mutant > 0, ]
  ends <- tapply(arms$resident, arms$branch, max)
  expect_equal(length(ends), 2L)
  expect_equal(unname(ends[1]), unname(ends[2]))
  vertex <- arms[arms$resident == max(arms$resident), ]
  expect_equal(nrow(vertex), 2L)
  expect_equal(vertex$mutant[1], vertex$mutant[2])
  expect_false(vertex$resident[1] %in% pip$residents$resident)
  # the mirror fold on x < 0 has its arms starting together, and is joined too
  mirror <- br[br$branch != diag_id & br$resident < 0 & br$mutant < 0, ]
  starts <- tapply(mirror$resident, mirror$branch, min)
  expect_equal(length(starts), 2L)
  expect_equal(unname(starts[1]), unname(starts[2]))
  mv <- mirror[mirror$resident == min(mirror$resident), ]
  expect_equal(nrow(mv), 2L)
  expect_equal(mv$mutant[1], mv$mutant[2])
  expect_equal(unname(mv$resident[1]), -unname(vertex$resident[1]), tolerance = 1e-6)
  expect_equal(unname(mv$mutant[1]), -unname(vertex$mutant[1]), tolerance = 1e-6)
})

test_that("a contour leaving through the edge of the mutant range is drawn to the edge", {
  # DD99's second contour y = 1.38 x leaves the mutant range |y| <= 2 at |x| = 1.45,
  # between residents of a coarse grid
  pip <- community_pip(dd99(), control = pip_control(list(n_resident = 9, n_mutant = 41, refine = 0)))
  br <- pip_branches(pip)
  on_diag <- abs(br$resident - br$mutant) < 1e-8
  other <- br[br$branch != as.integer(names(which.max(table(br$branch[on_diag])))), ]
  expect_equal(range(other$mutant), c(-2, 2), tolerance = 1e-8)
  expect_equal(other$mutant, dd99_y2(other$resident), tolerance = 1e-6)
})
