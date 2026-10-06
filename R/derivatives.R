# Derivatives of invasion fitness and of the selection gradient.
#
# Every method that needs a derivative -- the selection gradient, the mutant
# Hessian behind the ESS test, the resident Jacobian behind convergence
# stability -- obtains it here and never from a stencil of its own. This is the
# one place that decides where a derivative comes from, so that a model able to
# supply exact derivatives (automatic differentiation in plant, closed forms in
# the reference models) can do so without any consumer changing. Today the only
# source is the finite-difference machinery in util_gradient.R.
#
# Conventions: derivatives are with respect to raw trait values, and consumers
# apply any trait-scale chain rule themselves. `y` is a set of mutant traits
# (one row per point) evaluated against the community's current residents; the
# resident Jacobian differentiates with respect to the resident trait, which
# means re-solving the demographic equilibrium at each stencil point.

##' Control how derivatives are computed.
##'
##' Returns the list of settings read by the derivative functions
##' (\code{\link{community_fitness_gradient}},
##' \code{\link{community_fitness_hessian}},
##' \code{\link{community_selection_gradient_jacobian}}). Pass a list of values
##' to override the defaults; the result lives at
##' \code{community$derivative_control} (see \code{\link{community_start}}).
##'
##' Finite-difference steps are relative (\code{d_*}) with an absolute
##' fallback (\code{eps_*}) for trait values that are numerically zero, as the
##' singular strategies of the reference models are. First derivatives of
##' fitness use the \code{*_gradient} settings; the Hessian and the resident
##' Jacobian use \code{*_second}. The \code{r_*} settings are the number of
##' successively halved steps combined by Richardson extrapolation; each level of
##' the Jacobian costs \code{2k} demographic equilibrium solves, so it stays at
##' one by default.
##'
##' @title Derivative settings
##' @param control A list of values to modify from the defaults.
##' @return A list with elements \code{d_gradient}, \code{eps_gradient},
##' \code{r_gradient}, \code{d_second}, \code{eps_second}, \code{r_hessian},
##' \code{r_jacobian}.
##' @author Daniel Falster
##' @export
derivative_control <- function(control = NULL) {
  defaults <- list(
    d_gradient   = 1e-4,
    eps_gradient = 1e-4,
    r_gradient   = 1L,
    d_second     = 1e-3,
    eps_second   = 1e-3,
    r_hessian    = 2L,
    r_jacobian   = 1L
  )

  control <- as.list(control)
  extra <- setdiff(names(control), names(defaults))
  if (length(extra) > 0L) {
    stop("Unknown control parameters ", paste(extra, collapse = ", "))
  }
  ret <- modifyList(defaults, control)

  for (nm in c("d_gradient", "eps_gradient", "d_second", "eps_second")) {
    v <- ret[[nm]]
    if (!is.numeric(v) || length(v) != 1L || !is.finite(v) || v <= 0) {
      stop(nm, " must be a single positive number")
    }
  }
  for (nm in c("r_gradient", "r_hessian", "r_jacobian")) {
    v <- ret[[nm]]
    if (!is.numeric(v) || length(v) != 1L || v < 1 || v != round(v)) {
      stop(nm, " must be a positive whole number")
    }
    ret[[nm]] <- as.integer(v)
  }
  ret
}

## Communities built before derivative_control existed, or by hand in tests,
## carry no settings; they get the defaults.
community_derivative_control <- function(community) {
  ctrl <- community$derivative_control
  if (is.null(ctrl)) derivative_control() else ctrl
}

## Mutant points as a matrix with one row per point, following the convention
## of fitness_function(): a bare vector is one point per element for a
## single-trait community and a single point otherwise.
derivative_points <- function(y, trait_names) {
  k <- length(trait_names)
  if (!is.matrix(y)) {
    y <- if (k == 1L) matrix(y, ncol = 1L) else matrix(y, nrow = 1L)
  }
  if (ncol(y) != k) {
    stop("Expected points with ", k, " trait column(s), got ", ncol(y))
  }
  if (nrow(y) == 0L) {
    stop("No trait points to differentiate at")
  }
  colnames(y) <- trait_names
  y
}

community_fitness_function <- function(community) {
  f <- community$fitness_function
  if (!is.function(f)) {
    stop("The community has no fitness function; run community_demography() first")
  }
  f
}

##' Gradient of invasion fitness in the mutant direction.
##'
##' The derivative of invasion fitness \eqn{s(y)} with respect to the mutant
##' trait \eqn{y}, with the residents and their environment held fixed.
##' Evaluated at the residents' own traits this is the selection gradient
##' (\code{\link{community_selection_gradient}}); evaluated elsewhere it is the
##' slope of the fitness landscape. All stencil points are evaluated in a
##' single call to the community's vectorised fitness function.
##'
##' @title Fitness gradient in the mutant direction
##' @param community A \code{community} solved to demographic equilibrium.
##' @param y Trait points to differentiate at, one per row (a bare vector is
##' one point per element for a single trait, or a single point otherwise).
##' Defaults to the residents' traits.
##' @return A \code{nrow(y)} by \code{k} matrix of gradients with the fitness
##' at each point in \code{attr(., "value")}.
##' @author Daniel Falster
##' @export
community_fitness_gradient <- function(community, y = community$traits) {
  y <- derivative_points(y, community$trait_names)
  ctrl <- community_derivative_control(community)
  g <- fd_fitness_gradient(community_fitness_function(community), y, ctrl)
  dimnames(g) <- list(NULL, community$trait_names)
  g
}

fd_fitness_gradient <- function(f, y, ctrl) {
  n <- nrow(y)
  k <- ncol(y)
  stencils <- lapply(seq_len(n), function(i) {
    gradient_points(y[i, ], eps = ctrl$eps_gradient, d = ctrl$d_gradient,
                    r = ctrl$r_gradient)
  })
  pts <- rbind(y, do.call(rbind, stencils))
  ff <- as.numeric(f(pts))
  if (length(ff) != nrow(pts)) {
    stop("fitness_function must return one value per row of its argument")
  }
  if (!all(is.finite(ff))) {
    stop("Non-finite fitness at a finite-difference point")
  }

  g <- matrix(NA_real_, n, k)
  pos <- n
  for (i in seq_len(n)) {
    m <- nrow(stencils[[i]])
    yi <- ff[pos + seq_len(m)]
    dim(yi) <- attr(stencils[[i]], "dim_y")
    g[i, ] <- gradient_extrapolate(yi, stencils[[i]])
    pos <- pos + m
  }
  attr(g, "value") <- ff[seq_len(n)]
  g
}

##' Hessian of invasion fitness in the mutant direction.
##'
##' The curvature of invasion fitness \eqn{s(y)} with respect to the mutant
##' trait, residents held fixed. At a singular strategy a negative-definite
##' Hessian is the condition for evolutionary stability
##' (\code{\link{community_classify_singularity}}). All stencil points are
##' evaluated in one call to the fitness function.
##'
##' @title Fitness Hessian in the mutant direction
##' @param community A \code{community} solved to demographic equilibrium.
##' @param y A single trait point; defaults to the sole resident's traits.
##' @return A symmetric \code{k} by \code{k} matrix.
##' @author Daniel Falster
##' @export
community_fitness_hessian <- function(community, y = community$traits) {
  y <- derivative_points(y, community$trait_names)
  if (nrow(y) != 1L) {
    stop("community_fitness_hessian takes a single trait point; got ", nrow(y))
  }
  ctrl <- community_derivative_control(community)
  H <- util_hessian(community_fitness_function(community), y[1, ],
                    d = ctrl$d_second, eps = ctrl$eps_second, r = ctrl$r_hessian)
  dimnames(H) <- list(community$trait_names, community$trait_names)
  H
}

##' Jacobian of the selection gradient with respect to the resident trait.
##'
##' How the selection gradient changes as the (single) resident moves: the
##' matrix \eqn{J_{ij} = d g_i / d x_j}. Eigenvalues with negative real parts
##' mean selection carries a nearby resident towards the point --- convergence
##' stability. Each stencil point introduces the perturbed trait as the sole
##' resident and solves the community to demographic equilibrium, so this is
##' the expensive derivative: \code{2k} equilibrium solves per Richardson level.
##'
##' @title Resident Jacobian of the selection gradient
##' @param community A \code{community} with exactly one resident.
##' @param birth_rate Birth rate to start each equilibrium solve from; defaults
##' to the resident's own, with each solve then warm-starting from the last.
##' @return A \code{k} by \code{k} matrix with the selection gradient at the
##' resident in \code{attr(., "selection_gradient")}.
##' @author Daniel Falster
##' @export
community_selection_gradient_jacobian <- function(community, birth_rate = NULL) {
  trait_names <- community$trait_names
  if (nrow(community$traits) != 1L) {
    stop("community_selection_gradient_jacobian needs exactly one resident; ",
         "this community has ", nrow(community$traits))
  }
  x <- as.numeric(community$traits[1, ])
  ctrl <- community_derivative_control(community)

  gradient <- singularity_gradient_fn(community, birth_rate = birth_rate)
  g0 <- gradient(x)
  J <- util_jacobian(gradient, x, d = ctrl$d_second, eps = ctrl$eps_second,
                     r = ctrl$r_jacobian)
  dimnames(J) <- list(trait_names, trait_names)
  names(g0) <- trait_names
  attr(J, "selection_gradient") <- g0
  J
}
