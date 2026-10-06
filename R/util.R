copy_attributes <- function(from, to, exclude=character(0)) {
  at_from <- attributes(from)
  for (i in setdiff(names(at_from), c("dim", "dimnames"))) {
    attr(to, i) <- at_from[[i]]
  }
  to
}


util_rbind_list <- function(x) do.call("rbind", as.list(x))

## TODO: Get some testing in here so I can drop the asserts and merge
## with closest2
##
## TODO:
closest <- function(pt, A, r=NULL) {
  ## stopifnot(ncol(A) == length(pt))
  if (!is.null(r)) {
    ## stopifnot(is.matrix(r) &&
    ##                         ncol(r) == 2L &&
    ##                         nrow(r) == ncol(A))
    dr <- r[,2] - r[,1]
    A <- t((t(A) - r[,1]) / dr)
    pt <- (pt - r[,1]) / dr
    ## TODO:?
    ## A <- rescale(A, r)
    ## pt <- rescale(pt, r)
  }
  d <- colSums((t(A) - drop(pt))^2)
  i <- which.min(d)
  attr(i, "distance") <- sqrt(d[i])
  i
}

closest_log <- function(pt, A, r=NULL) {
  closest(log(pt), log(A), if (is.null(r)) NULL else log(r))
}

## Possible rewrite?
## closest2 <- function(pt, A) {
##   d <- colSums((t(A) - drop(pt))^2)
##   i <- which.min(d)
##   attr(i, "distance") <- sqrt(d[i])
##   i
## }

## Same rescaling issue here as for determining distance to points in
## the maximisation.
rescale <- function(A, r=NULL) {
  if (is.null(r)) {
    r <- t(apply(A, 2, range))
    colnames(r) <- c("lower", "upper")
  }
  stopifnot(is.matrix(r) &&
                          ncol(r) == 2L &&
                          nrow(r) == ncol(A))
  dr <- r[,2] - r[,1]
  t((t(A) - r[,1]) / dr)
}

unrescale <- function(A, r) {
  stopifnot(is.matrix(r) &&
            ncol(r) == 2L &&
            nrow(r) == ncol(A))
  dr <- r[,2] - r[,1]
  t(t(A) * dr + r[,1])
}

maximize_logspace <- function(f, x, bounds, tol) {
  control <- list(tol=tol, maximize=TRUE)
  log_bounds <- log(bounds)
  fit <- dfoptim::nmkb(log(x), function(x) f(exp(x)),
                       lower=log_bounds[,1], upper=log_bounds[,2],
                       control=list(tol=tol, maximize=TRUE))
  fit$par <- exp(fit$par)
  fit
}

## Trait-space transform used when spacing/searching during assembly. Strictly
## positive traits (plant) are best handled on a log scale; traits that span
## zero (e.g. DD99/GK98) need a linear scale. Returns forward/inverse maps that
## are applied to trait values and bounds throughout the assembly machinery
## (fitness landscape grid, nearest-resident distance, fitness maximisation).
community_trait_transform <- function(community) {
  scale <- if (is.null(community$trait_scale)) "log" else community$trait_scale
  if (identical(scale, "linear")) {
    list(fwd = function(x) x, inv = function(x) x, scale = "linear")
  } else {
    list(fwd = log, inv = exp, scale = "log")
  }
}

## The package's local maximiser: maximise f over `bounds`, searching in the
## coordinates of the trait transform `tf`. One trait is a golden-section
## search over the whole of `bounds` to tolerance `tol` (`x` is not needed);
## more traits are a coarse bounded Nelder-Mead from `x` (whose initial simplex
## steps off a stationary start, such as a resident at a singular point, where
## a gradient method would stay) polished by L-BFGS-B. Values of f that are not
## finite (inviable traits) count as very low fitness.
## Returns list(par, value) with par on the raw trait scale.
maximize_scaled <- function(f, x, bounds, tf, tol = .Machine$double.eps^0.25) {
  lb <- tf$fwd(bounds)
  g <- function(z) {
    w <- f(tf$inv(z))
    if (is.finite(w)) w else -1e6
  }
  if (nrow(lb) == 1L) {
    out <- optimise(g, lb[1, ], maximum = TRUE, tol = tol)
    return(list(par = tf$inv(out$maximum), value = out$objective))
  }
  ## nmkb needs a strictly-interior start; a resident sitting on a bound (common
  ## early in assembly) would otherwise produce NaNs. Inset slightly.
  span <- lb[, 2] - lb[, 1]
  z0 <- pmin(pmax(tf$fwd(x), lb[, 1] + 1e-6 * span), lb[, 2] - 1e-6 * span)
  coarse <- dfoptim::nmkb(z0, g, lower = lb[, 1], upper = lb[, 2],
                          control = list(tol = 1e-2, maximize = TRUE))
  fit <- optim(coarse$par, g, method = "L-BFGS-B",
               lower = lb[, 1], upper = lb[, 2],
               control = list(fnscale = -1, factr = 1e10))
  list(par = tf$inv(fit$par), value = fit$value)
}

has_attr <- function(x, which) {
  !is.null(attr(x, which, exact=TRUE))
}

norm2 <- function(x) {
  sqrt(sum(x^2))
}

vnapply <- function(X, FUN, ...) {
  vapply(X, FUN, numeric(1), ...)
}

last <- function(x) {
  x[[length(x)]]
}


