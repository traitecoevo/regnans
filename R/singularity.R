# Singular strategies in N dimensions: solving and classifying.
#
# community_solve_singularity_1D() (R/solve_attractors.R) brackets the scalar
# selection gradient with uniroot(). That is robust but strictly one-trait. The
# functions here are dimension-agnostic:
#
#   community_solve_singularity()     multivariate root-find on the selection
#                                     gradient (which is already dimension
#                                     agnostic, see community_selection_gradient)
#   community_classify_singularity()  second-order conditions at a singular
#                                     point: is it a CSS, a branching point, a
#                                     repeller, or a Garden of Eden?
#
# Both work through the harness connectors only, so they run on the fast toy
# harnesses (milliseconds, analytic answers) exactly as they run on the plant
# SCM.

## Strip residents from a community, keeping bounds/controls/harness. Used so
## that the singularity machinery always evaluates a *monomorphic* resident at
## the candidate trait, however the incoming community was built.
community_clear_residents <- function(community) {
  community$traits <- trait_matrix(numeric(0), community$trait_names)
  community$birth_rate <- numeric(0)
  community_reset(community)
}

## Every candidate trait is solved to equilibrium from some starting birth rate.
## Warm-starting matters: `birth_rate_initial` is a fixed 1e-3, but equilibrium
## birth rates are model- and parameterisation-dependent and can be five orders
## of magnitude away from it, in which case every solve crawls up from scratch
## and repeatedly trips the `equilibrium_large_birth_rate_change` schedule reset.
## Successive candidates sit close together in trait space, so the previous
## solve's equilibrium is a far better guess than a constant.
singularity_seed_birth_rate <- function(community, birth_rate = NULL) {
  if (!is.null(birth_rate)) {
    return(as.numeric(birth_rate))
  }
  br <- community$birth_rate
  if (nrow(community$traits) == 1L && length(br) == 1L &&
      is.finite(br) && br > 0) as.numeric(br) else NULL
}

## Closure: resident trait vector -> selection gradient at that resident.
## Each call introduces the trait as the sole resident, solves the community to
## demographic equilibrium, and differentiates invasion fitness in the mutant
## direction. The community behind the most recent call is kept so callers can
## return it rather than re-solving.
singularity_gradient_fn <- function(community, dx = 1e-4, birth_rate = NULL) {
  base <- community_clear_residents(community)
  trait_names <- community$trait_names
  last_community <- NULL
  seed_birth_rate <- singularity_seed_birth_rate(community, birth_rate)

  fn <- function(x) {
    out <- base |>
      community_add(trait_matrix(x, trait_names),
                    birth_rate = seed_birth_rate) |>
      community_demography() |>
      community_selection_gradient(dx = dx)
    last_community <<- out

    ## carry this equilibrium forward as the next candidate's starting point
    br <- as.numeric(out$birth_rate)
    if (length(br) == 1L && is.finite(br) && br > 0) {
      seed_birth_rate <<- br
    }

    as.numeric(out$selection_gradient)
  }
  attr(fn, "last") <- function() last_community
  fn
}

## Normalise a bounds argument to a k x 2 matrix. Accepts what
## community_solve_singularity_1D() accepts (a bare length-2 vector) when there
## is a single trait.
singularity_bounds <- function(bounds, trait_names) {
  k <- length(trait_names)
  if (!is.matrix(bounds)) {
    if (k != 1L || length(bounds) != 2L) {
      stop("bounds must be a ", k, " x 2 matrix (one row per trait)")
    }
    bounds <- matrix(bounds, nrow = 1L)
  }
  if (nrow(bounds) != k || ncol(bounds) != 2L) {
    stop("bounds must be a ", k, " x 2 matrix (one row per trait)")
  }
  rownames(bounds) <- trait_names
  colnames(bounds) <- c("lower", "upper")
  bounds
}

##' Find a singular strategy in any number of trait dimensions.
##'
##' A singular strategy is a resident trait combination at which the selection
##' gradient vanishes. \code{\link{community_solve_singularity_1D}} finds one in
##' a single trait by bracketing the scalar gradient with \code{uniroot}; this
##' function generalises it to \code{k} traits by handing the vector-valued
##' selection gradient to a multivariate root finder (\code{nleqslv} or
##' \code{dfsane}, via \code{\link{util_nlsolve}}).
##'
##' Each residual evaluation introduces the candidate trait combination as the
##' sole resident, solves the community to demographic equilibrium, and
##' differentiates invasion fitness in the mutant direction --- so it costs one
##' full demographic solve plus a gradient stencil. Any residents on the
##' incoming community are discarded (the singular point is monomorphic); if
##' there is exactly one, its traits are used as the starting point.
##'
##' The search runs on the community's trait scale (see
##' \code{community_start(trait_scale = )}): for \code{"log"} traits the root is
##' sought in \code{log(x)} and the residual is the gradient with respect to
##' \code{log(x)}, which is far better conditioned for strictly positive
##' biological traits. The root is the same either way. Candidate points are
##' clamped to \code{bounds}, and hitting a bound produces a warning --- as with
##' the 1-D solver, that means the search region did not contain a singularity.
##'
##' @title Solve for a singular strategy (N-dimensional)
##' @param community A \code{community} object to search within.
##' @param x0 Starting trait combination. Defaults to the traits of a single
##' resident if the community has one, otherwise the midpoint of \code{bounds}
##' on the community's trait scale.
##' @param bounds A \code{k} by 2 matrix of lower/upper bounds (a length-2
##' vector is accepted when there is a single trait). Defaults to the
##' community's bounds.
##' @param solver Root finder: \code{"nleqslv"} (default) or \code{"dfsane"}.
##' @param tol Convergence tolerance passed to the solver.
##' @param maxit Maximum solver iterations.
##' @param dx Step size for the selection-gradient finite differences.
##' @param birth_rate Birth rate to start each candidate's equilibrium solve
##' from. Defaults to the resident's own birth rate if the community has one,
##' otherwise \code{birth_rate_initial}; thereafter each solve warm-starts from
##' the previous one. Worth setting when the model's equilibrium birth rates are
##' far from \code{birth_rate_initial}.
##' @param edge_ok Is it (not) an error if the solution lands on the edge of
##' \code{bounds}?
##' @return The \code{community} at the singular strategy: one resident, solved
##' to demographic equilibrium, with \code{selection_gradient} set. Attributes
##' \code{converged}, \code{solver} and \code{singularity} (the root, named by
##' trait) record the solve.
##' @seealso \code{\link{community_classify_singularity}} to determine whether
##' the point found is a CSS, a branching point, a repeller or a Garden of Eden.
##' @author Daniel Falster
##' @export
community_solve_singularity <- function(community, x0 = NULL, bounds = NULL,
                                        solver = c("nleqslv", "dfsane"),
                                        tol = 1e-6, maxit = 100, dx = 1e-4,
                                        birth_rate = NULL, edge_ok = TRUE) {
  solver <- match.arg(solver)
  trait_names <- community$trait_names
  k <- length(trait_names)

  if (is.null(bounds)) {
    bounds <- community$bounds
  }
  bounds <- singularity_bounds(bounds, trait_names)

  tf <- community_trait_transform(community)
  z_lo <- tf$fwd(bounds[, 1])
  z_hi <- tf$fwd(bounds[, 2])
  if (any(!is.finite(z_lo) | !is.finite(z_hi))) {
    stop("community_solve_singularity needs finite bounds on the trait scale")
  }

  if (is.null(x0)) {
    if (nrow(community$traits) == 1L) {
      x0 <- as.numeric(community$traits[1, ])
    } else {
      x0 <- tf$inv((z_lo + z_hi) / 2)
    }
  }
  x0 <- as.numeric(x0)
  if (length(x0) != k) {
    stop("x0 must have one value per trait (", k, ")")
  }

  plant_log_assembler(sprintf(
    "Solving %dD singularity for [%s] from [%s] using %s",
    k, paste(trait_names, collapse = ", "),
    paste(signif(x0, 5), collapse = ", "), solver))

  gradient <- singularity_gradient_fn(community, dx = dx,
                                      birth_rate = birth_rate)

  ## Residual in search coordinates z: dS/dz = dS/dx * dx/dz. For a log trait
  ## scale dx/dz = x, so the residual is the gradient with respect to log(x) --
  ## the natural scale on which to ask whether selection has stopped.
  residual <- function(z) {
    z <- pmin(pmax(z, z_lo), z_hi)
    x <- tf$inv(z)
    dxdz <- if (identical(tf$scale, "log")) x else rep(1, k)
    gradient(x) * dxdz
  }

  z0 <- tf$fwd(x0)
  sol <- util_nlsolve(z0, residual, tol = tol, maxit = maxit, solver = solver,
                      require_converged = FALSE)
  converged <- isTRUE(attr(sol, "converged"))

  z_root <- pmin(pmax(as.numeric(sol), z_lo), z_hi)
  on_edge <- (z_root <= z_lo) | (z_root >= z_hi)
  if (any(on_edge)) {
    ## Candidates are clamped to the bounds, so a search region that does not
    ## contain a singularity ends up pinned against an edge -- which also makes
    ## the residual locally flat and the solver report failure. Report the
    ## informative cause, not the symptom.
    msg <- paste("Bounds do not include a singularity for trait(s)",
                 paste(trait_names[on_edge], collapse = ", "))
    if (edge_ok) warning(msg) else stop(msg)
  } else if (!converged) {
    warning(sprintf("community_solve_singularity did not converge (%s: %s)",
                    solver, attr(sol, "message")))
  }

  ## Re-evaluate at the root so the returned community is exactly the one at
  ## the singular point (the solver's last probe need not be).
  residual(z_root)
  out <- attr(gradient, "last")()

  root <- tf$inv(z_root)
  names(root) <- trait_names

  attr(out, "converged") <- converged
  attr(out, "solver") <- solver
  attr(out, "singularity") <- root

  plant_log_assembler(sprintf(
    "Solved! %dD singularity for [%s] is [%s]",
    k, paste(trait_names, collapse = ", "),
    paste(signif(root, 6), collapse = ", ")))

  out
}

##' Classify a singular strategy.
##'
##' Finding a singular strategy is only half the question: the interesting part
##' is what happens near it. Two independent second-order conditions decide
##' that (Geritz et al. 1998; Leimar 2009):
##'
##' \describe{
##'   \item{Evolutionary stability}{the curvature of invasion fitness in the
##'     \emph{mutant} direction, with the resident held at the singular point.
##'     In one trait this is the scalar
##'     \eqn{\partial^2 s / \partial y^2}; in \code{k} traits it is the
##'     \code{k x k} Hessian \eqn{H}. The point is an ESS (uninvadable) when
##'     \eqn{H} is negative definite, i.e. every eigenvalue is negative.}
##'   \item{Convergence stability}{whether selection carries a nearby resident
##'     towards the point. This is the derivative of the selection gradient with
##'     respect to the \emph{resident} trait: the scalar \eqn{dG/dx} in one
##'     trait, the Jacobian \eqn{J} in \code{k}. Eigenvalues of \eqn{J} with
##'     negative real parts give convergence stability under the canonical
##'     equation with an isotropic mutational covariance; a negative-definite
##'     symmetric part \eqn{(J + J^T)/2} gives \emph{strong} convergence
##'     stability, which holds for any mutational covariance matrix.}
##' }
##'
##' Crossing the two gives the standard four-way classification:
##'
##' \tabular{lll}{
##'   \tab \strong{convergence stable} \tab \strong{not convergence stable} \cr
##'   \strong{ESS} \tab CSS \tab Garden of Eden \cr
##'   \strong{not ESS} \tab branching point \tab repeller
##' }
##'
##' The full eigen-decompositions are returned, not just the verdict: when the
##' point is not an ESS the leading eigenvector of \eqn{H} is the direction in
##' trait space along which the population disruptively splits, which is itself
##' the scientific result in a multi-trait problem.
##'
##' Cost: the Hessian needs \code{1 + 4k^2} invasion-fitness evaluations, all
##' made in a single vectorised call against the cached resident environment, so
##' it is cheap. The Jacobian needs \code{2k} \emph{resident} evaluations, each a
##' full demographic equilibrium solve, so it dominates.
##'
##' @title Classify a singular strategy (1-D and N-D)
##' @param community A \code{community} with a single resident at (or very near)
##' a singular strategy --- typically the output of
##' \code{\link{community_solve_singularity}} or
##' \code{\link{community_solve_singularity_1D}}.
##' @param dx Step size for the selection-gradient finite differences.
##' @param d Relative finite-difference step for the Hessian and Jacobian.
##' @param eps Absolute finite-difference step, used for traits that are ~0
##' (the usual case for the toy harnesses, whose singular strategies sit at 0).
##' @param r_hessian Number of halved step sizes to Richardson-extrapolate the
##' Hessian over. Cheap, so 2 by default.
##' @param r_jacobian As \code{r_hessian} for the Jacobian. Each level costs
##' \code{2k} demographic solves, so 1 by default.
##' @param birth_rate Birth rate to start each resident equilibrium solve from;
##' see \code{\link{community_solve_singularity}}. Defaults to the singular
##' resident's own equilibrium birth rate, which is normally what you want.
##' @param tol Magnitude below which an eigenvalue counts as zero, making the
##' classification degenerate rather than forcing a verdict.
##' @return An object of class \code{singularity_classification}: a list with
##' the traits and selection gradient at the point, \code{hessian} and
##' \code{jacobian} with their eigen-decompositions
##' (\code{hessian_eigen}, \code{jacobian_eigen}, \code{jacobian_symmetric_eigen}),
##' the logical verdicts \code{evolutionarily_stable},
##' \code{convergence_stable} and \code{strongly_convergence_stable}, the
##' \code{branching_direction} where the point is invadable, and the four-way
##' \code{classification}.
##' @author Daniel Falster
##' @export
community_classify_singularity <- function(community, dx = 1e-4,
                                           d = 1e-3, eps = 1e-3,
                                           r_hessian = 2L, r_jacobian = 1L,
                                           birth_rate = NULL, tol = 1e-8) {

  trait_names <- community$trait_names
  k <- length(trait_names)

  if (nrow(community$traits) != 1L) {
    stop("community_classify_singularity needs exactly one resident (the ",
         "singular strategy); this community has ", nrow(community$traits))
  }
  x <- as.numeric(community$traits[1, ])

  if (is.null(community$fitness_function)) {
    community <- community_demography(community)
  }

  plant_log_assembler(sprintf(
    "Classifying %dD singularity at [%s] = [%s]",
    k, paste(trait_names, collapse = ", "),
    paste(signif(x, 6), collapse = ", ")))

  ## --- evolutionary stability: curvature in the mutant direction ------------
  ## community$fitness_function holds the resident fixed at x, so this is the
  ## Hessian of s(y; x) with respect to the mutant trait y, evaluated at y = x.
  H <- util_hessian(community$fitness_function, x, d = d, eps = eps,
                    r = r_hessian)
  dimnames(H) <- list(trait_names, trait_names)
  H_sym <- (H + t(H)) / 2
  H_eigen <- eigen(H_sym, symmetric = TRUE)

  ## --- convergence stability: how the gradient responds to the resident -----
  gradient <- singularity_gradient_fn(community, dx = dx,
                                      birth_rate = birth_rate)
  g0 <- gradient(x)
  J <- util_jacobian(gradient, x, d = d, eps = eps, r = r_jacobian)
  dimnames(J) <- list(trait_names, trait_names)
  J_eigen <- eigen(J)
  J_sym <- (J + t(J)) / 2
  J_sym_eigen <- eigen(J_sym, symmetric = TRUE)

  ## --- verdicts -------------------------------------------------------------
  ev_H <- H_eigen$values
  ev_J <- Re(J_eigen$values)
  ev_Js <- J_sym_eigen$values

  ess <- all(ev_H < -tol)
  cs <- all(ev_J < -tol)
  scs <- all(ev_Js < -tol)

  degenerate <- any(abs(ev_H) <= tol) || any(abs(ev_J) <= tol)

  classification <-
    if (degenerate) {
      "degenerate"
    } else if (ess && cs) {
      "CSS"
    } else if (!ess && cs) {
      "branching point"
    } else if (ess && !cs) {
      "Garden of Eden"
    } else {
      "repeller"
    }

  ## Where the point is invadable, the leading eigenvector of the Hessian is
  ## the trait-space direction along which fitness curves upward -- the
  ## direction a branching population splits along.
  branching_direction <- NULL
  if (!ess && max(ev_H) > tol) {
    branching_direction <- H_eigen$vectors[, which.max(ev_H)]
    ## eigen() fixes the sign arbitrarily; make it reproducible by pointing the
    ## largest component positive (the direction is an axis, not an arrow).
    branching_direction <-
      branching_direction * sign(branching_direction[which.max(abs(branching_direction))])
    names(branching_direction) <- trait_names
  }

  names(x) <- trait_names
  names(g0) <- trait_names

  ret <- list(
    trait_names = trait_names,
    traits = x,
    selection_gradient = g0,
    hessian = H,
    hessian_eigen = H_eigen,
    jacobian = J,
    jacobian_eigen = J_eigen,
    jacobian_symmetric_eigen = J_sym_eigen,
    evolutionarily_stable = ess,
    convergence_stable = cs,
    strongly_convergence_stable = scs,
    branching_direction = branching_direction,
    degenerate = degenerate,
    classification = classification,
    tol = tol
  )
  class(ret) <- "singularity_classification"

  plant_log_assembler(sprintf("Classified! Singularity at [%s] is a %s",
                              paste(signif(x, 6), collapse = ", "),
                              classification))

  ret
}

##' @param x A \code{singularity_classification} object.
##' @param ... Ignored.
##' @rdname community_classify_singularity
##' @export
print.singularity_classification <- function(x, ...) {
  cat(sprintf("<singularity_classification: %s>\n", x$classification))
  cat(sprintf("  traits: %s\n",
              paste(sprintf("%s = %s", x$trait_names, signif(x$traits, 6)),
                    collapse = ", ")))
  cat(sprintf("  selection gradient: %s\n",
              paste(signif(x$selection_gradient, 3), collapse = ", ")))
  cat(sprintf("  evolutionarily stable: %s (Hessian eigenvalues %s)\n",
              x$evolutionarily_stable,
              paste(signif(x$hessian_eigen$values, 4), collapse = ", ")))
  cat(sprintf("  convergence stable:    %s (Jacobian eigenvalues %s)\n",
              x$convergence_stable,
              paste(signif(Re(x$jacobian_eigen$values), 4), collapse = ", ")))
  cat(sprintf("  strongly conv. stable: %s\n", x$strongly_convergence_stable))
  if (!is.null(x$branching_direction)) {
    cat(sprintf("  branching direction:   %s\n",
                paste(sprintf("%s = %s", x$trait_names,
                              signif(x$branching_direction, 4)),
                      collapse = ", ")))
  }
  invisible(x)
}
