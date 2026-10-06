# The canonical equation with branching, against the reference models.
#
# JJ12 has a CSS at x* = x_opt - a sigma^2 and never branches. DD99 converges to
# x0 and branches there iff sigma_C < sigma_K; by symmetry about x0 the two
# daughters of an immediate split evolve as mirror images. GK98 at d/sigma =
# 1.5 branches at 0 to a protected dimorphism (-a*, a*), where a* can be read
# independently off the trait-evolution plot. The waiting time to branching has
# an oracle too: for a small kernel, lambda = rate * n * H * sd^2 with H the
# fitness curvature at the minimum, since p_est = 2 s and s = H delta^2 / 2.

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
events_of <- function(ce, what) ce$events[ce$events$event == what, ]

test_that("canonical_control validates its settings", {
  ctrl <- canonical_control()
  expect_equal(ctrl$stepper, "rodas")
  expect_equal(ctrl$branch, "expected")
  expect_equal(ctrl$t_max, Inf)
  expect_equal(ctrl$establishment_factor, 2)
  expect_equal(canonical_control(list(stepper = "rkck"))$stepper, "rkck")
  expect_equal(canonical_control(list(branch = FALSE))$branch, "none")      # the old logical still works
  expect_equal(canonical_control(list(branch = TRUE))$branch, "immediate")
  expect_error(canonical_control(list(stepper = "euler")), "should be one of")
  expect_error(canonical_control(list(branch = "sometimes")), "should be one of")
  expect_error(canonical_control(list(rate = 0)), "positive number")
  expect_error(canonical_control(list(t_max = -1)), "t_max must be positive")
  expect_error(canonical_control(list(max_residents = 1.5)), "whole number")
  expect_error(canonical_control(list(vcv = 1:3)), "square matrix")
  expect_error(canonical_control(list(immigration = list(pool = NULL))), "non-negative rate")
  expect_error(canonical_control(list(immigration = list(rate = 1, pool = 1:3))), "pool must be")
  expect_equal(canonical_control(list(immigration = list(rate = 1)))$immigration$establish, "deterministic")
  expect_error(canonical_control(list(nope = 1)), "Unknown control")
})

test_that("gauss_hermite integrates exactly against exp(-x^2)", {
  gh <- regnans:::gauss_hermite(12)
  expect_equal(sum(gh$w), sqrt(pi), tolerance = 1e-10)
  expect_equal(sum(gh$w * gh$x^2), sqrt(pi) / 2, tolerance = 1e-10)
  expect_equal(sum(gh$w * gh$x^4), 3 * sqrt(pi) / 4, tolerance = 1e-10)
  expect_lt(abs(sum(gh$w * gh$x)), 1e-12)
})

# ---- the equation ---------------------------------------------------------

test_that("JJ12 converges to its CSS and stops, with no branching", {
  ce <- community_canonical_equation(jj12(), x0 = -1.5)
  expect_s3_class(ce, "canonical_equation")
  expect_equal(ce$outcome, "stable")
  x_star <- 0.5 - 0.1 * 1.4^2
  expect_equal(final(ce)$x, x_star, tolerance = 1e-3)
  expect_equal(nrow(final(ce)), 1L)
  expect_equal(nrow(events_of(ce, "branch")), 0L)
  expect_true(all(diff(ce$trajectory$x) >= -1e-10))              # monotone approach from below
  expect_equal(names(ce$trajectory), c("time", "lineage", "x", "density", "gradient_x"))
  expect_lt(abs(final(ce)$gradient_x), canonical_control()$gradient_tol)
  expect_match(paste(utils::capture.output(print(ce)), collapse = ""), "stable after")
})

test_that("the rate scales evolutionary time and Cash-Karp agrees with RODAS", {
  slow <- community_canonical_equation(jj12(), x0 = -1.5, control = canonical_control(list(rate = 1)))
  fast <- community_canonical_equation(jj12(), x0 = -1.5, control = canonical_control(list(rate = 2)))
  cross <- function(ce) { tr <- ce$trajectory; tr$time[which(tr$x > 0)[1]] }
  expect_equal(cross(slow) / cross(fast), 2, tolerance = 0.1)
  ros <- community_canonical_equation(jj12(), x0 = -1.5, control = canonical_control(list(stepper = "rkck")))
  expect_equal(ros$outcome, "stable")
  expect_equal(final(ros)$x, final(slow)$x, tolerance = 1e-3)
})

test_that("DD99 with a wide competition kernel converges to x0 and is stable", {
  ce <- community_canonical_equation(dd99(sigma_C = 1.5, x0 = 0.3), x0 = -1)
  expect_equal(ce$outcome, "stable")
  expect_equal(final(ce)$x, 0.3, tolerance = 1e-3)
  expect_equal(nrow(events_of(ce, "branch")), 0L)
})

test_that("an immediate split at x0 gives mirror-image daughters, finished by Newton", {
  ce <- community_canonical_equation(dd99(sigma_C = 0.4, x0 = 0), x0 = 0.8,
                                     control = canonical_control(list(branch = "immediate", max_residents = 2)))
  br <- events_of(ce, "branch")
  expect_equal(nrow(br), 1L)
  at <- ce$trajectory[ce$trajectory$time == br$time, ]
  expect_equal(nrow(at), 2L)
  expect_equal(mean(at$x), 0, tolerance = 1e-3)
  expect_equal(abs(at$x), rep(0.02 * 4, 2), tolerance = 1e-2)
  fin <- final(ce)
  expect_equal(nrow(fin), 2L)
  expect_equal(fin$x[1], -fin$x[2], tolerance = 1e-3)
  expect_gt(abs(fin$x[1] - fin$x[2]), 0.1)
  expect_equal(fin$density[1], fin$density[2], tolerance = 1e-3)
  expect_true(ce$outcome %in% c("stable", "max_residents"))
  expect_lt(max(abs(fin$gradient_x)), canonical_control()$gradient_tol)
  expect_true(any(ce$events$event == "polish"))
  expect_lt(ce$evaluations, 600L)
})

test_that("without polishing the explicit stepper jitters at its stability limit and never converges", {
  unpolished <- community_canonical_equation(dd99(sigma_C = 0.4, x0 = 0), x0 = 0.8,
    control = canonical_control(list(stepper = "rkck", branch = "immediate", max_residents = 2, polish = FALSE, max_steps = 400)))
  expect_equal(unpolished$outcome, "max_steps")
  polished <- community_canonical_equation(dd99(sigma_C = 0.4, x0 = 0), x0 = 0.8,
    control = canonical_control(list(stepper = "rkck", branch = "immediate", max_residents = 2)))
  expect_gt(unpolished$evaluations, 5L * polished$evaluations)
  rodas <- community_canonical_equation(dd99(sigma_C = 0.4, x0 = 0), x0 = 0.8,
    control = canonical_control(list(branch = "immediate", max_residents = 2, polish = FALSE, max_steps = 400)))
  expect_equal(rodas$outcome, "max_residents")
  expect_lt(rodas$evaluations, 300L)
})

test_that("GK98 branches to the dimorphic coalition the trait-evolution plot shows", {
  comm <- gk98(d = 1.5)
  ce <- community_canonical_equation(comm, x0 = 1.2,
                                     control = canonical_control(list(branch = "immediate", max_residents = 2)))
  expect_equal(nrow(events_of(ce, "branch")), 1L)
  fin <- final(ce)
  expect_equal(nrow(fin), 2L)
  a_star <- mean(abs(fin$x))
  expect_equal(fin$x[1], -fin$x[2], tolerance = 1e-3)
  tep <- community_tep(comm, community_pip(comm, control = pip_control(list(n_resident = 41, n_mutant = 121, refine = 1))))
  sym <- tep[abs(tep$x1 + tep$x2) < 1e-8, ]
  sym <- sym[order(sym$x2), ]
  flip <- which(diff(sign(sym$g2)) != 0)[1]
  expect_equal(a_star, (sym$x2[flip] + sym$x2[flip + 1L]) / 2, tolerance = 0.1)
})

test_that("branching can be switched off, and the limits are reported as outcomes", {
  ce <- community_canonical_equation(dd99(sigma_C = 0.4), x0 = 0.8, control = canonical_control(list(branch = "none")))
  expect_equal(ce$outcome, "stable")
  expect_equal(nrow(final(ce)), 1L)
  short <- community_canonical_equation(jj12(), x0 = -1.5, control = canonical_control(list(t_max = 0.5)))
  expect_equal(short$outcome, "t_max")
  expect_equal(max(short$trajectory$time), 0.5)
  few <- community_canonical_equation(jj12(), x0 = -1.5, control = canonical_control(list(max_steps = 3)))
  expect_equal(few$outcome, "max_steps")
  expect_equal(few$steps, 3L)
})

# ---- waiting time to branching ----------------------------------------------

test_that("the branching rate matches its small-kernel oracle and the exact integral", {
  comm <- dd99(sigma_C = 0.4, x0 = 0) |> community_add(trait_matrix(0, "x"), birth_rate = 500) |> community_demography()
  base <- regnans:::community_clear_residents(comm)
  tf <- regnans:::community_trait_transform(comm)
  H <- 1 / 0.4^2 - 1                                       # r (1/sigma_C^2 - 1/sigma_K^2)
  ctrl <- canonical_control(list(rate = 1, branch_nodes = 12))
  for (sd in c(0.02, 0.1, 0.3)) {
    br <- regnans:::canonical_branch_rate(base, comm, tf, ctrl, matrix(0, 1, 1), 1L, 1, sd, 1)
    expect_equal(br$failed, 0L)
    # the exact integrand: kernel x min(1, 2 s); every mutant of a fitness
    # minimum coexists with the resident here
    exact <- stats::integrate(function(d) stats::dnorm(d, 0, sd) * pmin(1, 2 * comm$fitness_function(d)),
                              -8 * sd, 8 * sd, rel.tol = 1e-8)$value
    # within the kernel's reach s < 1/2 the integrand is smooth and twelve nodes
    # are exact to 1e-7; at sd = 0.3 the cap min(1, 2s) puts a kink in it and
    # the quadrature is good to about 1%
    expect_equal(br$rate, exact, tolerance = if (sd <= 0.1) 1e-5 else 1e-2, info = paste("sd =", sd))
    if (sd <= 0.02) expect_equal(br$rate, H * sd^2, tolerance = 2e-3)
  }
  # the rate scales with mu and with the resident's weight
  ctrl2 <- canonical_control(list(rate = 3, branch_nodes = 12))
  expect_equal(regnans:::canonical_branch_rate(base, comm, tf, ctrl2, matrix(0, 1, 1), 1L, 1, 0.02, 0.5)$rate,
               1.5 * H * 0.02^2, tolerance = 2e-3)
})

test_that("the time to branching grows with a narrower mutational kernel and a flatter minimum", {
  run <- function(sigma_C, sd) {
    community_canonical_equation(dd99(sigma_C = sigma_C, x0 = 0), x0 = 0.8,
      control = canonical_control(list(max_residents = 2, mutation_sd = sd, branch_nodes = 8)))
  }
  wait_of <- function(ce) {
    tr <- ce$trajectory
    events_of(ce, "branch")$time[1] - min(tr$time[abs(tr$gradient_x) < canonical_control()$gradient_tol])
  }
  wide <- wait_of(run(0.4, 0.10)); narrow <- wait_of(run(0.4, 0.02)); flat <- wait_of(run(0.8, 0.10))
  expect_gt(wide, 0)
  expect_gt(narrow, wide)
  expect_gt(flat, wide)
  # lambda ~ sd^2: a five-fold narrower kernel waits about 25 times longer
  expect_equal(narrow / wide, 25, tolerance = 0.3)
})

test_that("the mutant appears at the kernel's successful distance and the parent stays put", {
  ce <- community_canonical_equation(dd99(sigma_C = 0.4, x0 = 0), x0 = 0.8,
    control = canonical_control(list(max_residents = 2, mutation_sd = 0.1)))
  b <- events_of(ce, "branch")$time[1]
  at <- ce$trajectory[ce$trajectory$time == b, ]
  expect_equal(nrow(at), 2L)
  parent <- at[at$lineage == 1, ]; mutant <- at[at$lineage == 2, ]
  expect_lt(abs(parent$x), 1e-3)
  expect_gt(abs(mutant$x), 0.03); expect_lt(abs(mutant$x), 0.4)
  fin <- final(ce)
  expect_equal(nrow(fin), 2L)
  expect_lt(max(abs(fin$gradient_x)), canonical_control()$gradient_tol)
})

test_that("stochastic branching is reproducible under a seed and leaves the caller's RNG alone", {
  ctrl <- canonical_control(list(branch = "stochastic", max_residents = 2, mutation_sd = 0.1, seed = 7))
  set.seed(1); before <- runif(1)
  set.seed(1)
  a <- community_canonical_equation(dd99(sigma_C = 0.4, x0 = 0), x0 = 0.8, control = ctrl)
  expect_equal(runif(1), before)
  b <- community_canonical_equation(dd99(sigma_C = 0.4, x0 = 0), x0 = 0.8, control = ctrl)
  expect_equal(a$events, b$events)
  expect_equal(a$trajectory, b$trajectory)
  br <- events_of(a, "branch")$time[1]
  expect_true(is.finite(br) && br > 0)
  other <- community_canonical_equation(dd99(sigma_C = 0.4, x0 = 0), x0 = 0.8,
             control = canonical_control(list(branch = "stochastic", max_residents = 2, mutation_sd = 0.1, seed = 8)))
  expect_false(isTRUE(all.equal(events_of(other, "branch")$time[1], br)))
})

# ---- immigration ------------------------------------------------------------

test_that("immigrants from a pool establish when they can invade and are recorded as new lineages", {
  ce <- community_canonical_equation(dd99(sigma_C = 0.4, x0 = 0), x0 = 0.8,
    control = canonical_control(list(branch = "none", t_max = 30, max_residents = 6, seed = 3,
                                     immigration = list(rate = 1))))
  expect_equal(ce$outcome, "t_max")
  expect_gt(ce$immigration_attempts, 0L)
  established <- events_of(ce, "immigrant")
  expect_gt(nrow(established), 0L)
  expect_true(all(established$lineage > 1L))
  expect_gt(nrow(final(ce)), 1L)
  expect_match(paste(utils::capture.output(print(ce)), collapse = ""), "immigrants established")
  for (i in seq_len(nrow(established))) {
    at <- ce$trajectory[ce$trajectory$time == established$time[i], ]
    expect_true(established$lineage[i] %in% at$lineage)
  }
})

test_that("into a CSS model, arrivals fail or displace the resident, and the community stays monomorphic", {
  ce <- community_canonical_equation(dd99(sigma_C = 1.5, x0 = 0), x0 = 0.8,
    control = canonical_control(list(branch = "none", t_max = 40, seed = 5, immigration = list(rate = 1))))
  expect_gt(ce$immigration_attempts, 0L)
  n_in <- nrow(events_of(ce, "immigrant")); n_out <- nrow(events_of(ce, "extinct"))
  expect_gt(n_in, 0L)
  expect_lte(abs(n_in - n_out), 1L)                     # at most the latest arrival still being displaced
  expect_lte(nrow(final(ce)), 2L)
  expect_lt(n_in, ce$immigration_attempts)              # and some arrivals failed
})

test_that("an arrival that could establish is refused and recorded when the community is full", {
  ce <- community_canonical_equation(dd99(sigma_C = 0.4, x0 = 0), x0 = 0.8,
    control = canonical_control(list(branch = "none", t_max = 30, max_residents = 1, seed = 3,
                                     immigration = list(rate = 1))))
  expect_gt(nrow(events_of(ce, "immigrant_refused")), 0L)
  expect_equal(nrow(events_of(ce, "immigrant")), 0L)
  expect_equal(nrow(final(ce)), 1L)
  expect_match(paste(utils::capture.output(print(ce)), collapse = ""), "refused at max_residents")
})

test_that("the pool can be a matrix or a function, and a zero rate means no arrivals", {
  pool <- matrix(c(-0.5, 0.5), ncol = 1)
  ce <- community_canonical_equation(dd99(sigma_C = 0.4, x0 = 0), x0 = 1.5,
    control = canonical_control(list(branch = "none", t_max = 20, seed = 2, max_residents = 4,
                                     immigration = list(rate = 2, pool = pool))))
  imm <- events_of(ce, "immigrant")
  expect_gt(nrow(imm), 0L)
  first <- ce$trajectory[ce$trajectory$time == imm$time[1] & ce$trajectory$lineage == imm$lineage[1], ]
  expect_true(abs(abs(first$x) - 0.5) < 1e-8)
  fn_pool <- function(n) matrix(rep(0.25, n), ncol = 1)
  ce2 <- community_canonical_equation(dd99(sigma_C = 0.4, x0 = 0), x0 = 1.5,
    control = canonical_control(list(branch = "none", t_max = 20, seed = 2, max_residents = 4,
                                     immigration = list(rate = 2, pool = fn_pool, establish = "stochastic"))))
  expect_gt(ce2$immigration_attempts, 0L)
  none <- community_canonical_equation(dd99(sigma_C = 0.4, x0 = 0), x0 = 1.5,
    control = canonical_control(list(branch = "none", immigration = list(rate = 0))))
  expect_equal(none$immigration_attempts, 0L)
  expect_equal(none$outcome, "stable")
})

test_that("mutation and immigration run together: both kinds of event occur and the run ends", {
  # branching needs the resident to reach stationarity between arrivals, so
  # the mutation rate must be high relative to the immigration rate here
  ce <- community_canonical_equation(dd99(sigma_C = 0.4, x0 = 0), x0 = 0.8,
    control = canonical_control(list(rate = 200, mutation_sd = 0.1, max_residents = 6, t_max = 80, seed = 1,
                                     immigration = list(rate = 0.3))))
  expect_true(ce$outcome %in% c("t_max", "max_steps"))
  expect_gt(nrow(events_of(ce, "branch")), 0L)
  expect_gt(nrow(events_of(ce, "immigrant")), 0L)
  expect_true(all(ce$trajectory$time <= 80 + 1e-9))
})

test_that("a rare arrival during a branching wait does not hang the run (review blocker)", {
  ce <- community_canonical_equation(dd99(sigma_C = 0.4, x0 = 0), x0 = 0.8,
    control = canonical_control(list(mutation_sd = 0.1, seed = 1, max_residents = 3, max_steps = 400, t_max = 5000,
                                     immigration = list(rate = 0.02))))
  expect_true(ce$outcome %in% c("t_max", "max_steps", "max_residents"))
  expect_true(is.finite(ce$evaluations))
  expect_true(all(diff(unique(ce$trajectory$time)) > 0))
})

# ---- the record ---------------------------------------------------------------

test_that("the community can be rebuilt at any recorded time and landscapes plotted along the way", {
  ce <- community_canonical_equation(dd99(sigma_C = 0.4, x0 = 0), x0 = 0.8,
    control = canonical_control(list(branch = "immediate", max_residents = 2)))
  b <- events_of(ce, "branch")$time[1]
  times <- unique(ce$trajectory$time)
  before <- canonical_community(ce, max(times[times < b]))
  at <- canonical_community(ce, b)
  expect_equal(nrow(before$traits), 1L)
  expect_equal(nrow(at$traits), 2L)
  expect_lt(attr(before, "time"), b)
  expect_equal(attr(at, "time"), b)
  expect_lt(abs(as.numeric(before$resident_fitness)), 1e-6)
  f <- before$fitness_function
  expect_gt(f(0.3) + f(-0.3) - 2 * f(as.numeric(before$traits)), 0)   # a fitness minimum
  p <- plot(ce, type = "landscapes")
  expect_s3_class(p, "ggplot")
  expect_equal(length(unique(ggplot2::ggplot_build(p)$layout$layout$PANEL)), 3L)
  p2 <- plot(ce, type = "landscapes", times = c(0, b))
  expect_equal(length(unique(ggplot2::ggplot_build(p2)$layout$layout$PANEL)), 2L)
})

test_that("the trajectory plot builds, with a line per lineage and a mark per branching", {
  ce <- community_canonical_equation(dd99(sigma_C = 0.4), x0 = 0.8,
                                     control = canonical_control(list(branch = "immediate", max_residents = 2)))
  p <- plot(ce)
  expect_s3_class(p, "ggplot")
  b <- ggplot2::ggplot_build(p)
  expect_equal(length(unique(b$data[[1]]$group)), 2L)
  expect_equal(nrow(b$data[[2]]), 1L)
  expect_s3_class(ggplot2::autoplot(ce), "ggplot")
})
