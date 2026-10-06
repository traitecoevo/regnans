# The canonical equation with branching, against the reference models.
#
# JJ12 has a CSS at x* = x_opt - a sigma^2 and never branches. DD99 converges to
# x0 and branches there iff sigma_C < sigma_K; by symmetry about x0 the two
# daughters evolve as mirror images. GK98 at d/sigma = 1.5 branches at 0 to a
# protected dimorphism (-a*, a*), where a* can be read independently off the
# trait-evolution plot as the point on the symmetric line where the gradient on
# the outer resident changes sign.

jj12 <- function(a = 0.1, sigma = 1.4, x_opt = 0.5) {
  community_start(bounds(x = c(-3, 3)), trait_scale = "linear",
                  harness = harness_jj12(a = a, x_opt = x_opt, sigma = sigma))
}
dd99 <- function(sigma_C = 0.4, x0 = 0) {
  community_start(bounds(x = c(-2, 2)), trait_scale = "linear",
                  harness = harness_dd99(x0 = x0, sigma_K = 1, sigma_C = sigma_C))
}
gk98 <- function(d = 1.5) {
  community_start(bounds(x = c(-3, 3)), trait_scale = "linear",
                  harness = harness_gk98(d = d, sigma = 1))
}
final <- function(ce) {
  tr <- ce$trajectory
  tr[tr$time == max(tr$time), ]
}

test_that("canonical_control validates its settings", {
  expect_equal(canonical_control()$stepper, "rk23")
  expect_equal(canonical_control(list(stepper = "rosenbrock"))$stepper, "rosenbrock")
  expect_error(canonical_control(list(stepper = "euler")), "should be one of")
  expect_error(canonical_control(list(rate = 0)), "positive number")
  expect_error(canonical_control(list(max_residents = 1.5)), "whole number")
  expect_error(canonical_control(list(vcv = 1:3)), "square matrix")
  expect_error(canonical_control(list(nope = 1)), "Unknown control")
})

test_that("JJ12 converges to its CSS and stops, with no branching", {
  ce <- community_canonical_equation(jj12(), x0 = -1.5)
  expect_s3_class(ce, "canonical_equation")
  expect_equal(ce$outcome, "stable")
  x_star <- 0.5 - 0.1 * 1.4^2
  expect_equal(final(ce)$x, x_star, tolerance = 1e-3)
  expect_equal(nrow(final(ce)), 1L)
  expect_false(any(ce$events$event == "branch"))
  # monotone approach from below
  expect_true(all(diff(ce$trajectory$x) >= -1e-10))
  expect_equal(names(ce$trajectory), c("time", "lineage", "x", "density", "gradient_x"))
  expect_lt(abs(final(ce)$gradient_x), canonical_control()$gradient_tol)
  expect_match(paste(utils::capture.output(print(ce)), collapse = ""), "stable after")
})

test_that("the rate scales evolutionary time and the Rosenbrock stepper agrees with rk23", {
  slow <- community_canonical_equation(jj12(), x0 = -1.5, control = canonical_control(list(rate = 1)))
  fast <- community_canonical_equation(jj12(), x0 = -1.5, control = canonical_control(list(rate = 2)))
  # time to pass x = 0 halves when the rate doubles
  cross <- function(ce) { tr <- ce$trajectory; tr$time[which(tr$x > 0)[1]] }
  expect_equal(cross(slow) / cross(fast), 2, tolerance = 0.1)

  ros <- community_canonical_equation(jj12(), x0 = -1.5, control = canonical_control(list(stepper = "rosenbrock")))
  expect_equal(ros$outcome, "stable")
  expect_equal(final(ros)$x, final(slow)$x, tolerance = 1e-3)
  expect_gt(ros$evaluations, 0L)
})

test_that("DD99 with a wide competition kernel converges to x0 and is stable", {
  ce <- community_canonical_equation(dd99(sigma_C = 1.5, x0 = 0.3), x0 = -1)
  expect_equal(ce$outcome, "stable")
  expect_equal(final(ce)$x, 0.3, tolerance = 1e-3)
  expect_false(any(ce$events$event == "branch"))
})

test_that("DD99 with a narrow competition kernel branches at x0 into a symmetric pair", {
  ce <- community_canonical_equation(dd99(sigma_C = 0.4, x0 = 0), x0 = 0.8,
                                     control = canonical_control(list(max_residents = 2, t_max = 500)))
  br <- ce$events[ce$events$event == "branch", ]
  expect_equal(nrow(br), 1L)
  at <- ce$trajectory[ce$trajectory$time == br$time, ]
  # the record at the branching time holds the two daughters, placed
  # branch_distance either side of the singular strategy x0 = 0
  expect_equal(nrow(at), 2L)
  expect_equal(mean(at$x), 0, tolerance = 1e-3)
  expect_equal(abs(at$x), rep(0.02 * 4, 2), tolerance = 1e-2)
  fin <- final(ce)
  expect_equal(nrow(fin), 2L)
  expect_equal(sort(fin$lineage), c(1L, 2L))
  expect_equal(fin$x[1], -fin$x[2], tolerance = 1e-3)  # mirror images about x0
  expect_gt(abs(fin$x[1] - fin$x[2]), 0.1)              # and genuinely apart
  expect_equal(fin$density[1], fin$density[2], tolerance = 1e-3)
  expect_true(ce$outcome %in% c("stable", "max_residents"))
  expect_lt(max(abs(fin$gradient_x)), canonical_control()$gradient_tol)
  # the slow tail was finished by Newton on the coalition, not integrated
  expect_true(any(ce$events$event == "polish"))
  expect_lt(ce$evaluations, 600L)
})

test_that("without polishing the pair's slow approach is integrated and costs far more", {
  ce <- community_canonical_equation(dd99(sigma_C = 0.4, x0 = 0), x0 = 0.8,
                                     control = canonical_control(list(max_residents = 2, t_max = 200, polish = FALSE)))
  expect_equal(ce$outcome, "t_max")
  expect_gt(ce$evaluations, 1000L)
})

test_that("GK98 branches to the dimorphic coalition the trait-evolution plot shows", {
  comm <- gk98(d = 1.5)
  ce <- community_canonical_equation(comm, x0 = 1.2,
                                     control = canonical_control(list(max_residents = 2, t_max = 500)))
  expect_equal(sum(ce$events$event == "branch"), 1L)
  fin <- final(ce)
  expect_equal(nrow(fin), 2L)
  a_star <- mean(abs(fin$x))
  expect_equal(fin$x[1], -fin$x[2], tolerance = 1e-3)
  # the TEP's symmetric line: the gradient on the outer resident changes sign at a*
  tep <- community_tep(comm, community_pip(comm, control = pip_control(list(n_resident = 41, n_mutant = 121, refine = 1))))
  sym <- tep[abs(tep$x1 + tep$x2) < 1e-8, ]
  sym <- sym[order(sym$x2), ]
  flip <- which(diff(sign(sym$g2)) != 0)[1]
  a_tep <- (sym$x2[flip] + sym$x2[flip + 1L]) / 2
  expect_equal(a_star, a_tep, tolerance = 0.1)
})

test_that("branching can be switched off, and the limits are reported as outcomes", {
  ce <- community_canonical_equation(dd99(sigma_C = 0.4), x0 = 0.8,
                                     control = canonical_control(list(branch = FALSE)))
  expect_equal(ce$outcome, "stable")
  expect_equal(nrow(final(ce)), 1L)
  short <- community_canonical_equation(jj12(), x0 = -1.5, control = canonical_control(list(t_max = 0.5)))
  expect_equal(short$outcome, "t_max")
  expect_equal(max(short$trajectory$time), 0.5)
  few <- community_canonical_equation(jj12(), x0 = -1.5, control = canonical_control(list(max_steps = 3)))
  expect_equal(few$outcome, "max_steps")
  expect_equal(few$steps, 3L)
})

test_that("the trajectory plot builds, with a line per lineage and a mark per branching", {
  ce <- community_canonical_equation(dd99(sigma_C = 0.4), x0 = 0.8,
                                     control = canonical_control(list(max_residents = 2, t_max = 500)))
  p <- plot(ce)
  expect_s3_class(p, "ggplot")
  b <- ggplot2::ggplot_build(p)
  expect_equal(length(unique(b$data[[1]]$group)), 2L)
  expect_equal(nrow(b$data[[2]]), 1L)
  expect_s3_class(ggplot2::autoplot(ce), "ggplot")
})
