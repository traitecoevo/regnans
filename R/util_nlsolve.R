##' Solve a nonlinear system.
##'
##' One wrapper over three root finders, so the equilibrium code chooses by
##' name: \code{nleqslv} (Broyden or Newton with its own line search),
##' \code{dfsane} (derivative-free spectral method) and \code{newton}, the
##' package's own Newton--Broyden iteration. \code{newton} evaluates a
##' finite-difference Jacobian once (or takes an exact one from \code{jac}, or
##' a starting approximation from \code{J0}), then updates it by Broyden's
##' rank-one formula after every step and refreshes it by finite differences
##' only when a line search fails or every \code{refresh} steps. Each Newton
##' step is damped by backtracking on the residual norm. Its final Jacobian is
##' returned in \code{attr(., "jacobian")}, which lets the next solve of a
##' nearby problem start from it and skip the finite-difference pass.
##'
##' @title Solve a nonlinear system
##' @param x Starting point.
##' @param fn Function to solve: returns a residual vector of the same length.
##' @param tol Tolerance on both the residual and the step (for
##' \code{nleqslv} this is \code{xtol} and \code{ftol}).
##' @param maxit Maximum number of iterations. The number of function
##' evaluations will likely exceed this.
##' @param solver \code{"nleqslv"}, \code{"dfsane"} or \code{"newton"}.
##' @param require_converged If \code{TRUE} (the default) a solver that fails to
##' converge is an error. Set \code{FALSE} to get the best point found back with
##' \code{attr(., "converged") == FALSE} and decide what to do with it.
##' @param jac Optional exact Jacobian function of \code{fn}, used by
##' \code{newton} and \code{nleqslv}.
##' @param J0 Optional starting Jacobian approximation for \code{newton}, e.g.
##' the one a previous nearby solve returned.
##' @param refresh For \code{newton}: recompute the Jacobian by finite
##' differences every this many steps (Broyden updates in between).
##' @param max_step For \code{newton}: the largest change in any coordinate a
##' single step may make, before the line search. Far from the root a Newton
##' step can be huge (in log density, a factor of \eqn{e^{1000}}); capping it
##' keeps the iterates where the residual is well behaved.
##' @return The solution with attributes \code{y} (residual), \code{iter},
##' \code{feval}, \code{code}, \code{message}, \code{converged},
##' \code{solver}, and for \code{newton} also \code{jacobian}.
##' @export
util_nlsolve <- function(x, fn, tol=1e-6, maxit=100, solver="nleqslv",
                         require_converged=TRUE, jac = NULL, J0 = NULL,
                         refresh = 10L, max_step = Inf) {
  solver <- match.arg(solver, c("nleqslv", "dfsane", "newton"))

  res <- switch(solver,
                nleqslv=util_nlsolve_nleqslv(x, fn, tol, maxit, jac = jac),
                dfsane=util_nlsolve_dfsane(x, fn, tol, maxit),
                newton=util_nlsolve_newton(x, fn, tol, maxit, jac = jac, J0 = J0,
                                           refresh = refresh, max_step = max_step),
                stop("Unknown solver ", solver))

  if (require_converged && !attr(res, "converged")) {
    stop(sprintf("Solver has likely failed: code=%d, msg: %s",
                 attr(res, "code"), attr(res, "message")),
         immediate.=TRUE)
  }

  res
}

util_nlsolve_nleqslv <- function(x, fn, tol=1e-6, maxit=100, jac = NULL) {
  control <- list(xtol=tol, ftol=tol, maxit=maxit)
  sol <- nleqslv::nleqslv(x, fn, jac = jac, global="none", control=control)
  code <- sol$termcd
  res <- sol$x
  attributes(res) <- util_nlsolve_nleqslv_attr(sol)
  res
}

util_nlsolve_nleqslv_attr <- function(sol) {
  list(y=sol$fvec, # different to dfsane
       iter=sol$iter,
       feval=sol$nfcnt, # does not include jacobian evals
       code=sol$termcd,
       message=sol$message,
       converged=!(sol$termcd > 2 || sol$termcd < 0),
       solver="nleqslv")
}

util_nlsolve_dfsane <- function(x, fn, tol=1e-6, maxit=100) {
  control <- list(tol=tol, maxit=maxit, trace=FALSE)
  ## This works around `is.vector`, which returns FALSE if x has any
  ## attribute, which confuses dfsane.
  fn_vector <- function(x) as.numeric(fn(x))
  sol <- BB::dfsane(x, fn_vector, control=control, quiet=TRUE)
  res <- sol$par
  attributes(res) <- util_nlsolve_dfsane_attr(sol)
  res
}

util_nlsolve_dfsane_attr <- function(sol) {
  list(y=sol$residual,
       iter=sol$iter,
       feval=sol$feval,
       code=sol$convergence,
       message=sol$message,
       converged=sol$convergence == 0,
       solver="dfsane")
}

failed <- function(x) {
  inherits(x, "try-error")
}

## Finite-difference Jacobian of fn at x, forward differences, one evaluation
## per coordinate; `fx` is fn(x) if already known.
util_fd_jacobian <- function(fn, x, fx = fn(x), h = 1e-6 * pmax(abs(x), 1)) {
  n <- length(x)
  J <- matrix(NA_real_, length(fx), n)
  for (j in seq_len(n)) {
    xj <- x
    xj[j] <- xj[j] + h[j]
    J[, j] <- (fn(xj) - fx) / h[j]
  }
  J
}

util_nlsolve_newton <- function(x, fn, tol = 1e-6, maxit = 100, jac = NULL,
                                J0 = NULL, refresh = 10L, max_step = Inf) {
  x <- as.numeric(x)
  n <- length(x)
  feval <- 0L
  f <- function(z) { feval <<- feval + 1L; as.numeric(fn(z)) }
  fx <- f(x)
  if (length(fx) != n) stop("util_nlsolve: fn must return one residual per unknown")

  jacobian <- function(z, fz) {
    if (is.function(jac)) {
      matrix(as.numeric(jac(z)), n, n)
    } else {
      feval <<- feval + n
      util_fd_jacobian(fn, z, fz)
    }
  }
  J <- if (!is.null(J0) && identical(dim(J0), c(n, n)) && all(is.finite(J0))) J0 else jacobian(x, fx)
  since_refresh <- 0L
  code <- 1L
  message <- "converged"
  iter <- 0L
  converged <- max(abs(fx)) <= tol

  while (!converged && iter < maxit) {
    iter <- iter + 1L
    step <- tryCatch(-solve(J, fx), error = function(e) NULL)
    if (is.null(step) || any(!is.finite(step))) {
      J <- jacobian(x, fx); since_refresh <- 0L
      step <- tryCatch(-solve(J, fx), error = function(e) NULL)
      if (is.null(step) || any(!is.finite(step))) { code <- 4L; message <- "singular Jacobian"; break }
    }
    if (max(abs(step)) > max_step) step <- step * (max_step / max(abs(step)))
    ## backtracking on the residual norm; a step that cannot reduce it with the
    ## current (Broyden) Jacobian gets one more try with a fresh one
    norm0 <- sum(fx^2)
    lambda <- 1
    accepted <- FALSE
    for (k in seq_len(12L)) {
      x_new <- x + lambda * step
      f_new <- f(x_new)
      if (all(is.finite(f_new)) && sum(f_new^2) <= (1 - 1e-4 * lambda) * norm0) { accepted <- TRUE; break }
      lambda <- lambda / 2
    }
    if (!accepted) {
      if (since_refresh > 0L) {
        J <- jacobian(x, fx); since_refresh <- 0L
        next
      }
      code <- 3L; message <- "line search failed to reduce the residual"; break
    }
    dx <- x_new - x
    df <- f_new - fx
    x <- x_new
    fx <- f_new
    since_refresh <- since_refresh + 1L
    if (max(abs(fx)) <= tol || max(abs(dx)) <= tol) { converged <- TRUE; break }
    if (since_refresh >= refresh) {
      J <- jacobian(x, fx); since_refresh <- 0L
    } else {
      J <- J + ((df - J %*% dx) %*% t(dx)) / sum(dx^2)
    }
  }
  if (!converged && code == 1L) { code <- 2L; message <- "iteration limit reached" }

  res <- x
  attributes(res) <- list(y = fx, iter = iter, feval = feval, code = code,
                          message = message, converged = converged,
                          solver = "newton", jacobian = J)
  res
}
