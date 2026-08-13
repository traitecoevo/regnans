# Internal finite-difference gradient helpers.
#
# Vendored from richfitz/grader (gradient_points / gradient_extrapolate) so the
# package does not depend on that tiny, non-CRAN package. Used by
# community_selection_gradient() in solve_attractors.R.
#
# gradient_points() builds the set of evaluation points (a sequence of
# successively halved step sizes per dimension) needed for a central-difference
# estimate with optional Richardson extrapolation. gradient_extrapolate() turns
# the fitness values at those points into the gradient.

gradient_points <- function(x, eps = 1e-04, d = 1e-04, r = 4,
                            zero_tol = sqrt(.Machine$double.eps / 7e-07)) {
  v <- 2
  n <- length(x)
  h <- abs(d * x) + eps * (abs(x) < zero_tol)
  pts <- vector("list", r * n)
  dim(pts) <- c(r, n)
  dx <- matrix(NA, r, n)
  j <- seq_len(n)
  for (k in seq_len(r)) {
    for (i in seq_len(n)) {
      dx_i <- h * (i == j)
      pts[[k, i]] <- rbind(x + dx_i, x - dx_i)
    }
    dx[k, ] <- 2 * h
    h <- h / v
  }
  ret <- do.call("rbind", pts)
  attr(ret, "dim_y") <- c(2L, r * n)
  attr(ret, "dx") <- dx
  attr(ret, "n") <- n
  attr(ret, "r") <- r
  class(ret) <- "gradient_points"
  ret
}

gradient_extrapolate <- function(y, pts) {
  dx <- attr(pts, "dx")
  r <- attr(pts, "r")
  n <- attr(pts, "n")
  a <- (y[1, ] - y[2, ]) / dx
  for (m in seq_len(r - 1L)) {
    four_m <- 4^m
    a_next <- matrix(NA, r - m, n)
    for (i in seq_len(nrow(a_next))) {
      a_next[i, ] <- (a[i + 1L, ] * four_m - a[i, ]) / (four_m - 1)
    }
    a <- a_next
  }
  drop(a)
}

## Step sizes used by the second-order helpers below. Same convention as
## gradient_points(): relative step `d` on each coordinate, falling back to the
## absolute step `eps` for coordinates that are (numerically) zero -- which is
## the usual case for the toy harnesses, whose singular strategies sit at 0.
util_fd_step <- function(x, d, eps,
                         zero_tol = sqrt(.Machine$double.eps / 7e-07)) {
  h <- abs(d * x) + eps * (abs(x) < zero_tol)
  if (any(h <= 0)) {
    stop("Non-positive finite-difference step; increase `eps`")
  }
  h
}

## Richardson extrapolation over a list of estimates computed at successively
## halved step sizes (all O(h^2) accurate), as in gradient_extrapolate().
util_richardson <- function(est) {
  while (length(est) > 1L) {
    nxt <- vector("list", length(est) - 1L)
    for (i in seq_along(nxt)) {
      nxt[[i]] <- (4 * est[[i + 1L]] - est[[i]]) / 3
    }
    est <- nxt
  }
  est[[1L]]
}

##' Hessian of a scalar function by central differences.
##'
##' Companion to \code{gradient_points()}/\code{gradient_extrapolate()}: the
##' second-order analogue used by \code{\link{community_classify_singularity}}
##' to get the curvature of invasion fitness in the mutant direction. All
##' evaluation points are assembled first and \code{f} is called \emph{once},
##' because a community's \code{fitness_function} is vectorised over mutants
##' (one SCM mutant sweep against a cached resident environment) and calling it
##' per point would be far slower.
##'
##' Diagonal entries use the standard three-point second difference and
##' off-diagonals the four-point mixed difference; both are O(h^2), so estimates
##' at \code{r} successively halved steps are combined by Richardson
##' extrapolation.
##'
##' @title Finite-difference Hessian
##' @param f Function taking a matrix of points (one per row) and returning a
##' numeric vector of values.
##' @param x Point at which to evaluate the Hessian.
##' @param d Relative step size.
##' @param eps Absolute step size, used for coordinates that are ~0.
##' @param r Number of successively halved step sizes to extrapolate over
##' (\code{r = 1} disables extrapolation).
##' @return A symmetric \code{length(x)} by \code{length(x)} matrix.
##' @author Daniel Falster
##' @noRd
util_hessian <- function(f, x, d = 1e-3, eps = 1e-3, r = 2L) {
  x <- as.numeric(x)
  n <- length(x)
  h0 <- util_fd_step(x, d, eps)
  hs <- lapply(seq_len(r), function(k) h0 / 2^(k - 1L))

  unit <- diag(1, n)

  ## Points: the centre once, then per level the +/- pairs (diagonal) and the
  ## four-point stencils (off-diagonal).
  pts <- list(x)
  for (h in hs) {
    for (i in seq_len(n)) {
      for (j in seq_len(i)) {
        hi <- h[i] * unit[i, ]
        if (i == j) {
          pts <- c(pts, list(x + hi, x - hi))
        } else {
          hj <- h[j] * unit[j, ]
          pts <- c(pts, list(x + hi + hj, x + hi - hj,
                             x - hi + hj, x - hi - hj))
        }
      }
    }
  }

  P <- do.call(rbind, pts)
  y <- as.numeric(f(P))
  if (length(y) != nrow(P)) {
    stop("util_hessian: `f` must return one value per row of its argument")
  }
  if (!all(is.finite(y))) {
    stop("util_hessian: non-finite fitness at a finite-difference point")
  }

  y0 <- y[[1L]]
  pos <- 1L
  est <- vector("list", r)
  for (k in seq_len(r)) {
    h <- hs[[k]]
    H <- matrix(NA_real_, n, n)
    for (i in seq_len(n)) {
      for (j in seq_len(i)) {
        if (i == j) {
          H[i, i] <- (y[pos + 1L] - 2 * y0 + y[pos + 2L]) / (h[i] * h[i])
          pos <- pos + 2L
        } else {
          H[i, j] <- H[j, i] <-
            (y[pos + 1L] - y[pos + 2L] - y[pos + 3L] + y[pos + 4L]) /
            (4 * h[i] * h[j])
          pos <- pos + 4L
        }
      }
    }
    est[[k]] <- H
  }

  util_richardson(est)
}

##' Jacobian of a vector-valued function by central differences.
##'
##' Used by \code{\link{community_classify_singularity}} for the derivative of
##' the selection gradient with respect to the resident trait(s) --- the
##' convergence-stability condition. Unlike \code{\link{util_hessian}} each
##' evaluation here is expensive (one demographic equilibrium solve per point),
##' so \code{g} is called point by point and \code{r = 1} (a plain central
##' difference, \code{2 * length(x)} evaluations) is the default.
##'
##' @title Finite-difference Jacobian
##' @param g Function taking a numeric vector and returning a numeric vector.
##' @param x Point at which to evaluate the Jacobian.
##' @param d Relative step size.
##' @param eps Absolute step size, used for coordinates that are ~0.
##' @param r Number of successively halved step sizes to extrapolate over.
##' @return A matrix with \code{J[i, j] = d g_i / d x_j}.
##' @author Daniel Falster
##' @noRd
util_jacobian <- function(g, x, d = 1e-3, eps = 1e-3, r = 1L) {
  x <- as.numeric(x)
  n <- length(x)
  h0 <- util_fd_step(x, d, eps)
  unit <- diag(1, n)

  est <- vector("list", r)
  for (k in seq_len(r)) {
    h <- h0 / 2^(k - 1L)
    J <- NULL
    for (j in seq_len(n)) {
      hj <- h[j] * unit[j, ]
      up <- as.numeric(g(x + hj))
      dn <- as.numeric(g(x - hj))
      col <- (up - dn) / (2 * h[j])
      if (is.null(J)) {
        J <- matrix(NA_real_, length(col), n)
      }
      J[, j] <- col
    }
    est[[k]] <- J
  }

  util_richardson(est)
}
