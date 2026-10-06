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
  expect_equal(canonical_control()$branch, "expected")
  expect_error(canonical_control(list(branch = "sometimes")), "should be one of")
  expect_error(canonical_control(list(immigration = list(pool = NULL))), "non-negative rate")
  expect_error(canonical_control(list(immigration = list(rate = 1, pool = 1:3))), "pool must be")
  expect_equal(canonical_control(list(immigration = list(rate = 1)))$immigration$establish, "deterministic")
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
                                     control = canonical_control(list(branch = "immediate", max_residents = 2, t_max = 500)))
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
                                     control = canonical_control(list(branch = "immediate", max_residents = 2, t_max = 200, polish = FALSE)))
  expect_equal(ce$outcome, "t_max")
  polished <- community_canonical_equation(dd99(sigma_C = 0.4, x0 = 0), x0 = 0.8,
                                           control = canonical_control(list(branch = "immediate", max_residents = 2, t_max = 200)))
  expect_gt(ce$evaluations, 5L * polished$evaluations)
})

test_that("GK98 branches to the dimorphic coalition the trait-evolution plot shows", {
  comm <- gk98(d = 1.5)
  ce <- community_canonical_equation(comm, x0 = 1.2,
                                     control = canonical_control(list(branch = "immediate", max_residents = 2, t_max = 500)))
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
                                     control = canonical_control(list(branch = "none")))
  expect_equal(ce$outcome, "stable")
  expect_equal(canonical_control(list(branch = FALSE))$branch, "none")      # the old logical still works
  expect_equal(canonical_control(list(branch = TRUE))$branch, "immediate")
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
                                     control = canonical_control(list(branch = "immediate", max_residents = 2, t_max = 500)))
  p <- plot(ce)
  expect_s3_class(p, "ggplot")
  b <- ggplot2::ggplot_build(p)
  expect_equal(length(unique(b$data[[1]]$group)), 2L)
  expect_equal(nrow(b$data[[2]]), 1L)
  expect_s3_class(ggplot2::autoplot(ce), "ggplot")
})

# ---- waiting time to branching ----------------------------------------------

test_that("the time to branching grows with a narrower mutational kernel and a flatter minimum", {
  run <- function(sigma_C, sd) {
    community_canonical_equation(dd99(sigma_C = sigma_C, x0 = 0), x0 = 0.8,
      control = canonical_control(list(max_residents = 2, t_max = 1e6, mutation_sd = sd, branch_nodes = 8)))
  }
  arrive <- function(ce) {
    tr <- ce$trajectory
    c(stationary = min(tr$time[abs(tr$gradient_x) < canonical_control()$gradient_tol]),
      branch = ce$events$time[ce$events$event == "branch"][1])
  }
  wide <- arrive(run(0.4, 0.10))
  narrow <- arrive(run(0.4, 0.02))
  flat <- arrive(run(0.8, 0.10))
  # the community waits at the singular strategy before branching
  expect_gt(wide[["branch"]], wide[["stationary"]])
  # smaller mutations: longer wait; a shallower fitness minimum: longer wait
  expect_gt(narrow[["branch"]] - narrow[["stationary"]], wide[["branch"]] - wide[["stationary"]])
  expect_gt(flat[["branch"]] - flat[["stationary"]], wide[["branch"]] - wide[["stationary"]])
})

test_that("the mutant appears at the kernel's successful distance and the parent stays put", {
  ce <- community_canonical_equation(dd99(sigma_C = 0.4, x0 = 0), x0 = 0.8,
    control = canonical_control(list(max_residents = 2, t_max = 1e6, mutation_sd = 0.1)))
  b <- ce$events$time[ce$events$event == "branch"][1]
  at <- ce$trajectory[ce$trajectory$time == b, ]
  expect_equal(nrow(at), 2L)
  parent <- at[at$lineage == 1, ]; mutant <- at[at$lineage == 2, ]
  expect_lt(abs(parent$x), 1e-3)                                    # parent at x* = 0
  expect_gt(abs(mutant$x), 0.03); expect_lt(abs(mutant$x), 0.4)    # mutant within a few sd
  fin <- final(ce)
  expect_equal(nrow(fin), 2L)
  expect_lt(max(abs(fin$gradient_x)), canonical_control()$gradient_tol)
})

test_that("stochastic branching is reproducible under a seed and leaves the caller's RNG alone", {
  ctrl <- canonical_control(list(branch = "stochastic", max_residents = 2, t_max = 1e6, mutation_sd = 0.1, seed = 7))
  set.seed(1); before <- runif(1)
  set.seed(1)
  a <- community_canonical_equation(dd99(sigma_C = 0.4, x0 = 0), x0 = 0.8, control = ctrl)
  after <- runif(1)
  expect_equal(after, before)                                  # our seed did not disturb the stream
  b <- community_canonical_equation(dd99(sigma_C = 0.4, x0 = 0), x0 = 0.8, control = ctrl)
  expect_equal(a$events, b$events)
  expect_equal(a$trajectory, b$trajectory)
  br <- a$events$time[a$events$event == "branch"][1]
  expect_true(is.finite(br) && br > 0)
  other <- community_canonical_equation(dd99(sigma_C = 0.4, x0 = 0), x0 = 0.8,
             control = canonical_control(list(branch = "stochastic", max_residents = 2, t_max = 1e6, mutation_sd = 0.1, seed = 8)))
  expect_false(isTRUE(all.equal(other$events$time[other$events$event == "branch"][1], br)))
})

test_that("a CSS resident is reported unbranchable only when no coexisting mutant is reachable", {
  # DD99 with sigma_C > sigma_K: the singular strategy is an ESS, so no
  # branching candidate at all -- stable, no waiting
  ce <- community_canonical_equation(dd99(sigma_C = 1.5, x0 = 0.3), x0 = -1)
  expect_equal(ce$outcome, "stable")
  expect_false(any(ce$events$event %in% c("branch", "unbranchable")))
})

# ---- immigration ------------------------------------------------------------

test_that("immigrants from a pool establish when they can invade and are recorded as new lineages", {
  ce <- community_canonical_equation(dd99(sigma_C = 0.4, x0 = 0), x0 = 0.8,
    control = canonical_control(list(branch = "none", t_max = 30, max_residents = 6, seed = 3,
                                     immigration = list(rate = 1))))
  expect_gt(ce$immigration_attempts, 0L)
  established <- ce$events[ce$events$event == "immigrant", ]
  expect_gt(nrow(established), 0L)
  expect_true(all(established$lineage > 1L))
  expect_gt(nrow(final(ce)), 1L)
  expect_match(paste(utils::capture.output(print(ce)), collapse = ""), "immigrants established")
  # a trajectory record exists at each arrival, with the new lineage present
  for (i in seq_len(nrow(established))) {
    at <- ce$trajectory[ce$trajectory$time == established$time[i], ]
    expect_true(established$lineage[i] %in% at$lineage)
  }
})

test_that("into a CSS model, arrivals either fail or replace the resident, and the community stays monomorphic", {
  ce <- community_canonical_equation(dd99(sigma_C = 1.5, x0 = 0), x0 = 0.8,
    control = canonical_control(list(branch = "none", t_max = 40, seed = 5,
                                     immigration = list(rate = 1))))
  expect_gt(ce$immigration_attempts, 0L)
  expect_equal(nrow(final(ce)), 1L)
  n_in <- sum(ce$events$event == "immigrant"); n_out <- sum(ce$events$event == "extinct")
  expect_equal(n_in, n_out)                                     # every successful arrival displaced one
  expect_lt(n_in, ce$immigration_attempts)                      # and some arrivals failed
  expect_lt(abs(final(ce)$x), 0.2)                              # the survivor sits near x* = 0
})

test_that("the pool can be a matrix or a function, and a zero rate means no arrivals", {
  pool <- matrix(c(-0.5, 0.5), ncol = 1)
  ce <- community_canonical_equation(dd99(sigma_C = 0.4, x0 = 0), x0 = 1.5,
    control = canonical_control(list(branch = "none", t_max = 20, seed = 2, max_residents = 4,
                                     immigration = list(rate = 2, pool = pool))))
  imm <- ce$events[ce$events$event == "immigrant", ]
  expect_gt(nrow(imm), 0L)
  first <- ce$trajectory[ce$trajectory$time == imm$time[1] & ce$trajectory$lineage == imm$lineage[1], ]
  expect_true(abs(abs(first$x) - 0.5) < 1e-8)
  fn_pool <- function(n) matrix(rep(0.25, n), ncol = 1)
  ce2 <- community_canonical_equation(dd99(sigma_C = 0.4, x0 = 0), x0 = 1.5,
    control = canonical_control(list(branch = "none", t_max = 20, seed = 2, max_residents = 4,
                                     immigration = list(rate = 2, pool = fn_pool, establish = "stochastic"))))
  expect_gt(ce2$immigration_attempts, 0L)
  none <- community_canonical_equation(dd99(sigma_C = 0.4, x0 = 0), x0 = 1.5,
    control = canonical_control(list(branch = "none", t_max = 20, immigration = list(rate = 0))))
  expect_equal(none$immigration_attempts, 0L)
})

# ---- the record ---------------------------------------------------------------

test_that("the community can be rebuilt at any recorded time and landscapes plotted along the way", {
  ce <- community_canonical_equation(dd99(sigma_C = 0.4, x0 = 0), x0 = 0.8,
    control = canonical_control(list(branch = "immediate", max_residents = 2, t_max = 500)))
  b <- ce$events$time[ce$events$event == "branch"][1]
  times <- unique(ce$trajectory$time)
  before <- canonical_community(ce, max(times[times < b]))
  at <- canonical_community(ce, b)
  expect_equal(nrow(before$traits), 1L)
  expect_equal(nrow(at$traits), 2L)
  expect_lt(attr(before, "time"), b)
  expect_equal(attr(at, "time"), b)
  expect_lt(abs(as.numeric(before$resident_fitness)), 1e-6)
  # the landscape at the singular strategy is a fitness minimum: curvature up
  f <- before$fitness_function
  expect_gt(f(0.3) + f(-0.3) - 2 * f(as.numeric(before$traits)), 0)

  p <- plot(ce, type = "landscapes")
  expect_s3_class(p, "ggplot")
  built <- ggplot2::ggplot_build(p)
  expect_equal(length(unique(built$layout$layout$PANEL)), 3L)   # before, at, after the branching
  p2 <- plot(ce, type = "landscapes", times = c(0, b))
  expect_equal(length(unique(ggplot2::ggplot_build(p2)$layout$layout$PANEL)), 2L)
})
