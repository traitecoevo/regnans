# util_nlsolve (issue #27): thin wrapper over nleqslv::nleqslv and BB::dfsane.
# Previously untested because it was only reached through the equilibrium_solve_*
# demography path. Tested here directly on a small system with a known root, plus
# the non-convergence error path. Pure numerics -- no model, no SCM.

# System: x^2 + y^2 = 4, x*y = 1. The root near (1.5, 1) is
# (x, y) = ((1+sqrt(3))/sqrt(2), (sqrt(3)-1)/sqrt(2)) ~= (1.93185, 0.51764).
sys <- function(z) c(z[1]^2 + z[2]^2 - 4, z[1] * z[2] - 1)
root <- c((1 + sqrt(3)) / sqrt(2), (sqrt(3) - 1) / sqrt(2))

test_that("util_nlsolve finds a known root with every solver", {
  for (s in c("nleqslv", "dfsane", "newton")) {
    r <- util_nlsolve(c(1.5, 1), sys, solver = s)
    expect_equal(as.numeric(r), root, tolerance = 1e-4, info = s)
    expect_lt(max(abs(sys(r))), 1e-5)            # residual ~0
    expect_true(attr(r, "converged"), info = s)
    expect_equal(attr(r, "solver"), s)
  }
})

test_that("util_nlsolve respects the requested tolerance", {
  loose <- util_nlsolve(c(1.5, 1), sys, tol = 1e-3)
  tight <- util_nlsolve(c(1.5, 1), sys, tol = 1e-10)
  expect_lte(max(abs(sys(tight))), max(abs(sys(loose))) + 1e-12)
})

test_that("util_nlsolve rejects an unknown solver", {
  expect_error(util_nlsolve(c(1, 1), sys, solver = "secant"), "should be one of")
})

# ---- newton ------------------------------------------------------------------

sys_jac <- function(z) matrix(c(2 * z[1], z[2], 2 * z[2], z[1]), 2, 2)

test_that("newton converges quadratically, counts its evaluations and returns its Jacobian", {
  r <- util_nlsolve(c(1.5, 1), sys, solver = "newton", tol = 1e-12)
  expect_equal(as.numeric(r), root, tolerance = 1e-10)
  expect_lte(attr(r, "iter"), 8L)
  expect_equal(dim(attr(r, "jacobian")), c(2L, 2L))
  # the Broyden-updated Jacobian is a usable approximation of the true one
  # (exact only along the directions it has stepped in)
  expect_equal(attr(r, "jacobian"), sys_jac(root), tolerance = 0.2)
  # evaluations: one at the start, two for the finite-difference Jacobian, then
  # about one per step
  expect_lte(attr(r, "feval"), 3L + 2L * attr(r, "iter"))
})

test_that("an exact Jacobian or a carried-over one saves the finite-difference pass", {
  fd <- util_nlsolve(c(1.5, 1), sys, solver = "newton", tol = 1e-10)
  exact <- util_nlsolve(c(1.5, 1), sys, solver = "newton", tol = 1e-10, jac = sys_jac)
  carried <- util_nlsolve(c(1.5, 1), sys, solver = "newton", tol = 1e-10, jac = NULL,
                          J0 = attr(fd, "jacobian"))
  expect_equal(as.numeric(exact), root, tolerance = 1e-8)
  expect_equal(as.numeric(carried), root, tolerance = 1e-8)
  expect_lt(attr(exact, "feval"), attr(fd, "feval"))
  expect_lt(attr(carried, "feval"), attr(fd, "feval"))
  # a Jacobian of the wrong shape is ignored, not used
  wrong <- util_nlsolve(c(1.5, 1), sys, solver = "newton", J0 = matrix(1, 3, 3))
  expect_equal(as.numeric(wrong), root, tolerance = 1e-4)
  # nleqslv takes the exact Jacobian too
  expect_equal(as.numeric(util_nlsolve(c(1.5, 1), sys, jac = sys_jac)), root, tolerance = 1e-4)
})

test_that("newton reports failure honestly", {
  no_root <- function(z) c(z[1]^2 + 1, z[2]^2 + 1)
  expect_error(util_nlsolve(c(0.5, 0.5), no_root, solver = "newton", maxit = 10), "Solver has likely failed")
  r <- util_nlsolve(c(0.5, 0.5), no_root, solver = "newton", maxit = 10, require_converged = FALSE)
  expect_false(attr(r, "converged"))
  expect_true(attr(r, "code") > 1L)
})

test_that("util_nlsolve errors when the system has no solution", {
  # x^2 = -1, y^2 = -1 has no real root; the solver cannot converge.
  no_root <- function(z) c(z[1]^2 + 1, z[2]^2 + 1)
  expect_error(util_nlsolve(c(0, 0), no_root, maxit = 5), "Solver has likely failed")
})
