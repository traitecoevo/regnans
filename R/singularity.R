# Singular strategies in any number of dimensions: solving and classifying.
#
#   community_solve_singularity()     root-find on the selection gradient
#                                     (which is already dimension agnostic, see
#                                     community_selection_gradient)
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
singularity_seed_birth_rate <- function(community, m = 1L, birth_rate = NULL) {
  if (!is.null(birth_rate)) {
    return(as.numeric(birth_rate))
  }
  br <- community$birth_rate
  if (nrow(community$traits) == m && length(br) == m &&
      all(is.finite(br) & br > 0)) as.numeric(br) else NULL
}

## Residents stacked into one vector as an m x k trait matrix is stored, trait
## by trait; these label its entries ("x" for one resident, "x[1]", "x[2]", ...
## for several).
resident_labels <- function(trait_names, m) {
  if (m == 1L) trait_names else paste0(rep(trait_names, each = m), "[", seq_len(m), "]")
}

## Closure: the stacked traits of m residents -> the selection gradient of each
## at its own traits, stacked the same way. Each call introduces the traits as
## the residents, solves the community to demographic equilibrium (warm-started
## from the previous call's densities and solver state), and differentiates
## invasion fitness in the mutant direction; a call at the traits of the
## previous one returns its answer without solving again. Residents that have
## merged (residents_merged()) are outside the domain of a coalition: the
## answer is non-finite and nothing is solved, so a solver backtracks. The
## community behind the most recent solve is kept so callers can return it
## rather than re-solving, and the solves are counted.
##
## attr(fn, "points") evaluates a list of points, as a stencil asks for them:
## each solved from the densities and solver state of the last solve here
## (the stencil's centre), not from the previous point, so the points are
## independent and go through the parallel map, and the answer is the same
## with or without a plan. Each is solved to the stencil's tolerance
## (stencil_community()) and carries its equilibrium birth rates. The centre
## remains the starting point of the next call.
singularity_gradient_fn <- function(community, m = 1L, birth_rate = NULL) {
  base <- community_clear_residents(community)
  trait_names <- community$trait_names
  k <- length(trait_names)
  last_community <- NULL
  last_x <- NULL
  last_g <- NULL
  state <- NULL
  evaluations <- 0L
  refused <- 0L
  seed_birth_rate <- singularity_seed_birth_rate(community, m, birth_rate)
  ctrl <- community_derivative_control(community)

  fn <- function(x) {
    x <- as.numeric(x)
    if (identical(x, last_x)) {
      return(last_g)
    }
    if (nrow(residents_merged(x, m, k, ctrl)) > 0L) {
      refused <<- refused + 1L
      return(rep(NA_real_, m * k))
    }
    out <- singularity_solve_point(x, base, m, seed_birth_rate, state)
    evaluations <<- evaluations + 1L
    remember(x, out)
  }
  points <- function(points, step) {
    points <- lapply(points, as.numeric)
    merged <- vapply(points, function(x) nrow(residents_merged(x, m, k, ctrl)) > 0L, logical(1))
    out <- rep(list(structure(rep(NA_real_, m * k), birth_rate = rep(NA_real_, m))), length(points))
    stencil <- stencil_community(base, step)
    out[!merged] <- regnans_map(points[!merged], singularity_gradient_point, base = stencil, m = m,
                                birth_rate = seed_birth_rate, state = state)
    failed <- vapply(out, inherits, logical(1), "error")
    refused <<- refused + sum(merged)
    evaluations <<- evaluations + sum(!merged & !failed)
    if (any(failed)) stop(out[[which(failed)[1L]]])
    stencil_check_converged(out[!merged], stencil)
    out
  }
  ## keep a solved community as the answer at x and carry its equilibrium
  ## forward as the next candidate's starting point
  remember <- function(x, out) {
    last_community <<- out
    state <<- out$demography_state
    br <- as.numeric(out$birth_rate)
    if (length(br) == m && all(is.finite(br) & br > 0)) {
      seed_birth_rate <<- br
    }
    last_x <<- x
    last_g <<- as.numeric(out$selection_gradient)
    last_g
  }
  attr(fn, "points") <- points
  attr(fn, "last") <- function() last_community
  attr(fn, "evaluations") <- function() evaluations
  attr(fn, "refused") <- function() refused
  ## Take a community already at demographic equilibrium as the answer at its
  ## own traits, so a caller holding one does not pay to solve it again.
  attr(fn, "prime") <- function(community) {
    if (is.null(community$selection_gradient)) {
      community <- community_selection_gradient(community)
    }
    remember(as.numeric(community$traits), community)
  }
  fn
}

## The solve behind one call of the gradient closure: m residents at the
## stacked traits x added to the community without residents `base`, solved
## to demographic equilibrium from seed densities and the equilibrium solver's
## state, with the selection gradient at each.
singularity_solve_point <- function(x, base, m, birth_rate, state) {
  trait_names <- base$trait_names
  comm <- base |>
    community_add(trait_matrix(matrix(x, m, length(trait_names)), trait_names),
                  birth_rate = birth_rate)
  comm$demography_state <- state
  comm |>
    community_demography() |>
    community_selection_gradient()
}

## The same for a stencil point, perhaps on a worker, returning only numbers
## -- the gradient, with the equilibrium birth rates and whether the solve
## converged (a solved plant community cannot cross back to this process) --
## or the error, so that the solves which did finish are still counted.
singularity_gradient_point <- function(x, base, m, birth_rate, state) {
  tryCatch({
    comm <- singularity_solve_point(x, base, m, birth_rate, state)
    structure(as.numeric(comm$selection_gradient), birth_rate = as.numeric(comm$birth_rate),
              converged = isTRUE(attr(comm, "converged")))
  }, error = function(e) e)
}

## Has this community been solved to a demographic equilibrium that can stand
## for a fresh solve at its residents: its fitness function in place, the
## solve converged and every resident at a positive density.
community_at_equilibrium <- function(community) {
  br <- as.numeric(community$birth_rate)
  !is.null(community$fitness_function) && isTRUE(attr(community, "converged")) &&
    length(br) == nrow(community$traits) && all(is.finite(br) & br > 0)
}

## Normalise a bounds argument to a k x 2 matrix, accepting a bare length-2
## vector when there is a single trait.
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

## Root of a scalar residual bracketed between z_lo and z_hi. Without a sign
## change there is no root to bracket; the bound the residual points towards is
## returned unconverged, which the caller reports as a missed singularity.
singularity_bracket <- function(residual, z_lo, z_hi, tol, maxit) {
  f_lo <- residual(z_lo)
  f_hi <- residual(z_hi)
  if (!isTRUE(f_lo * f_hi <= 0)) {
    z <- if (f_lo < 0) z_lo else z_hi
    return(structure(z, converged = FALSE,
                     message = "no sign change across the bounds"))
  }
  out <- uniroot(residual, lower = z_lo, upper = z_hi, f.lower = f_lo,
                 f.upper = f_hi, tol = tol, maxiter = maxit)
  structure(out$root, converged = out$iter < maxit,
            message = sprintf("uniroot, %d iterations", out$iter))
}

##' Find a singular strategy or a singular coalition.
##'
##' A singular strategy is a resident trait combination at which the selection
##' gradient vanishes; a singular coalition is a set of coexisting residents at
##' each of which its own selection gradient vanishes, in the environment they
##' make together (the residents a population reaches after evolutionary
##' branching). The solver finds either: the stacked selection gradients of
##' all residents are handed to a root finder, the community being re-solved
##' to demographic equilibrium at every evaluation. One resident is the
##' monomorphic singular strategy.
##'
##' Root finders: \code{nleqslv} (Broyden with Newton restarts), \code{newton}
##' (the package's Newton--Broyden iteration with a backtracking line search,
##' see \code{\link{util_nlsolve}}) and \code{dfsane} (derivative-free), all
##' local searches from \code{x0}; or, for one resident with a single trait,
##' \code{"bracket"}, which brackets the scalar gradient between the bounds with
##' \code{uniroot} and so needs a sign change across them but then cannot miss
##' the root inside. \code{nleqslv} and \code{newton} take the Jacobian of the
##' gradients with respect to the residents from
##' \code{\link{community_selection_gradient_jacobian}} (finite differences
##' across the residents, \code{2mk} equilibrium solves, until the model
##' supplies equilibrium sensitivities), refreshing it only when their
##' rank-one updates fail; under a \code{future::plan()} its solves run in
##' parallel (\code{\link{derivative_control}}). A Jacobian already in hand for a nearby point
##' (\code{jacobian}, as in a continuation) answers their first request instead,
##' saving those solves.
##'
##' Each residual evaluation costs one demographic equilibrium solve, warm
##' started from the last, plus a gradient of invasion fitness at each
##' resident. The residents searched are those of \code{x0}, by default the
##' community's own; a community without residents is searched from the
##' midpoint of \code{bounds} for a monomorphic singular strategy.
##'
##' The search runs on the community's trait scale (see
##' \code{community_start(trait_scale = )}): for \code{"log"} traits the root is
##' sought in \code{log(x)} and the residual is the gradient with respect to
##' \code{log(x)}, which is far better conditioned for strictly positive
##' biological traits. The root is the same either way. Candidate points are
##' clamped to \code{bounds}, and hitting a bound produces a warning: the search
##' region did not contain a singularity. With \code{"bracket"} the same happens
##' when the gradient has one sign across the bounds, and the returned point is
##' the bound selection pushes towards. A coalition that loses a resident on
##' the way --- one whose equilibrium density is no longer positive, or two that
##' merge --- is not a coalition of that size: the solve is reported as not
##' converged, with a warning naming the residents.
##'
##' @title Solve for a singular strategy or coalition (N-dimensional)
##' @param community A \code{community} object to search within.
##' @param x0 Starting traits, one row per resident (a bare vector is one
##' resident with one value per trait, whatever the number of traits; several
##' residents need a matrix, e.g. \code{trait_matrix(c(-0.3, 0.6), "x")}).
##' Defaults to the community's residents, or, if it has none, the
##' midpoint of \code{bounds} on the community's trait scale.
##' @param bounds A \code{k} by 2 matrix of lower/upper bounds (a length-2
##' vector is accepted when there is a single trait). Defaults to the
##' community's bounds.
##' @param solver Root finder: \code{"nleqslv"} (default), \code{"newton"},
##' \code{"dfsane"}, or \code{"bracket"} (one resident, one trait).
##' @param tol Convergence tolerance passed to the solver: on the residual for
##' \code{"nleqslv"}, \code{"newton"} and \code{"dfsane"}, on the root (on the
##' trait scale) for \code{"bracket"}.
##' @param maxit Maximum solver iterations.
##' @param birth_rate Birth rates to start the first equilibrium solve from,
##' one per resident. Defaults to the residents' own if the community's
##' residents are being searched, otherwise \code{birth_rate_initial};
##' thereafter each solve warm-starts from the previous one. Worth setting when
##' the model's equilibrium birth rates are far from \code{birth_rate_initial}.
##' @param edge_ok Is it (not) an error if the solution lands on the edge of
##' \code{bounds}?
##' @param jacobian A resident Jacobian of the selection gradients, \code{mk}
##' by \code{mk} in raw trait units as
##' \code{\link{community_selection_gradient_jacobian}} returns it, to start
##' \code{"nleqslv"} or \code{"newton"} from in place of computing one at
##' \code{x0}: typically that of a nearby singular point. Later Jacobians are
##' computed as usual.
##' @return The \code{community} at the singular strategy or coalition, solved
##' to demographic equilibrium, with \code{selection_gradient} set. Attributes
##' record the solve: \code{converged}, \code{solver}, \code{evaluations} (the
##' equilibrium solves it cost) and \code{singularity} (the root: named by
##' trait for one resident, a residents-by-traits matrix for several).
##' @seealso \code{\link{community_classify_singularity}} to determine whether
##' the point found is a CSS, a branching point, a repeller or a Garden of Eden.
##' @examples
##' comm <- community_start(bounds(x = c(-2, 2)), trait_scale = "linear",
##'                         harness = harness_dd99(sigma_C = 0.4))
##' # the monomorphic singular strategy, then the dimorphic coalition it
##' # branches into
##' one <- community_solve_singularity(comm, x0 = 0.5)
##' pair <- community_solve_singularity(comm, x0 = trait_matrix(c(-0.3, 0.6), "x"))
##' pair$traits
##' community_classify_singularity(pair)
##' @author Daniel Falster
##' @export
community_solve_singularity <- function(community, x0 = NULL, bounds = NULL,
                                        solver = c("nleqslv", "newton", "dfsane", "bracket"),
                                        tol = 1e-6, maxit = 100,
                                        birth_rate = NULL, edge_ok = TRUE,
                                        jacobian = NULL) {
  solver <- match.arg(solver)
  trait_names <- community$trait_names
  k <- length(trait_names)

  if (solver == "bracket" && k != 1L) {
    stop("solver = \"bracket\" needs a single trait; this community has ", k)
  }

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
    x0 <- if (nrow(community$traits) >= 1L) community$traits else tf$inv((z_lo + z_hi) / 2)
  }
  ## a vector is one resident whatever the number of traits, so the same call
  ## cannot mean a coalition in one model and a single resident in another
  if (!is.matrix(x0)) {
    if (length(x0) != k) {
      stop("x0 as a vector is one resident, with one value per trait (", k, "); ",
           "give several residents as a matrix with one row per resident, e.g. ",
           "trait_matrix(c(...), \"", trait_names[1], "\")")
    }
    x0 <- matrix(x0, nrow = 1L)
  }
  if (ncol(x0) != k || nrow(x0) < 1L) {
    stop("x0 must have one column per trait (", k, ") and a row per resident")
  }
  m <- nrow(x0)
  if (solver == "bracket" && m != 1L) {
    stop("solver = \"bracket\" needs a single resident; x0 has ", m)
  }
  labels <- resident_labels(trait_names, m)
  if (!is.null(jacobian) && !solver %in% c("nleqslv", "newton")) {
    stop("jacobian is a starting point for the \"nleqslv\" and \"newton\" solvers; ",
         "\"", solver, "\" takes none")
  }
  if (!is.null(jacobian) &&
      (!is.numeric(jacobian) || !identical(dim(jacobian), c(m * k, m * k)) ||
       !all(is.finite(jacobian)))) {
    stop("jacobian must be a finite ", m * k, " x ", m * k,
         " matrix, the resident Jacobian of x0's residents")
  }
  ## the bounds, stacked as the residents' traits are
  z_lo <- rep(z_lo, each = m)
  z_hi <- rep(z_hi, each = m)

  plant_log_assembler(sprintf(
    "Solving %dD singularity%s for [%s] from [%s] using %s",
    k, if (m > 1L) sprintf(" (%d residents)", m) else "",
    paste(trait_names, collapse = ", "),
    paste(signif(x0, 5), collapse = ", "), solver))

  ctrl <- community_derivative_control(community)
  ## the search starts inside the bounds, as every candidate is clamped there;
  ## residents that a start outside them puts on the same edge have merged
  z0 <- pmin(pmax(as.numeric(tf$fwd(x0)), z_lo), z_hi)
  merged <- residents_merged(tf$inv(z0), m, k, ctrl)
  if (nrow(merged) > 0L) {
    stop("x0 has residents ", merged[1, 1], " and ", merged[1, 2],
         " at the same traits", if (!isTRUE(all.equal(tf$inv(z0), as.numeric(x0)))) " once clamped to the bounds",
         "; a coalition needs distinct residents")
  }
  gradient <- singularity_gradient_fn(community, m, birth_rate = birth_rate)
  log_scale <- identical(tf$scale, "log")

  ## Residual in search coordinates z: dS/dz = dS/dx * dx/dz. For a log trait
  ## scale dx/dz = x, so the residual is the gradient with respect to log(x) --
  ## the natural scale on which to ask whether selection has stopped.
  ## `refused_before` is the closure's count of refused (merged) probes as of
  ## the last finite residual, so an error can be traced to a refusal since.
  refused_before <- 0L
  residual <- function(z) {
    z <- pmin(pmax(z, z_lo), z_hi)
    x <- tf$inv(z)
    dxdz <- if (log_scale) x else rep(1, m * k)
    r <- gradient(x) * dxdz
    if (all(is.finite(r))) refused_before <<- attr(gradient, "refused")()
    r
  }
  ## Its Jacobian, the resident Jacobian carried onto the trait scale: on a
  ## log scale d(g_i x_i)/dz_j = x_i J_ij x_j + [i = j] g_i x_i. The solvers
  ## ask for it where they have just evaluated the residual, so the gradient
  ## there is already in hand. A Jacobian handed in is carried onto the trait
  ## scale at the start, as the solver's starting approximation.
  residual_jacobian <- function(z, J = NULL) {
    x <- tf$inv(pmin(pmax(z, z_lo), z_hi))
    g <- gradient(x)
    if (is.null(J)) J <- resident_jacobian(gradient, x, ctrl)
    if (log_scale) J <- J * outer(x, x) + diag(g * x, length(x))
    J
  }

  sol <- if (solver == "bracket") {
    singularity_bracket(residual, z_lo, z_hi, tol = tol, maxit = maxit)
  } else {
    ## a solver that cannot step back from a refused (merged) probe errors
    ## out of it; that is the coalition being lost, not a failure of the solve.
    ## Any other error is re-raised, carrying the equilibrium solves spent.
    tryCatch(
      util_nlsolve(z0, residual, tol = tol, maxit = maxit,
                   solver = solver, require_converged = FALSE,
                   jac = if (solver == "dfsane") NULL else residual_jacobian,
                   J0 = if (is.null(jacobian)) NULL else residual_jacobian(z0, unname(jacobian))),
      error = function(e) {
        if (attr(gradient, "refused")() == refused_before) {
          e$evaluations <- attr(gradient, "evaluations")()
          stop(e)
        }
        structure(as.numeric(tf$fwd(attr(gradient, "last")()$traits)), converged = FALSE,
                  merged = TRUE, message = conditionMessage(e))
      })
  }
  converged <- isTRUE(attr(sol, "converged"))

  z_root <- pmin(pmax(as.numeric(sol), z_lo), z_hi)
  on_edge <- (z_root <= z_lo) | (z_root >= z_hi)

  ## Re-evaluate at the root so the returned community is exactly the one at
  ## the singular point (the solver's last probe need not be). A coalition
  ## whose residents merged has no community there; the last one solved is
  ## returned, marked as not converged.
  residual(z_root)
  out <- attr(gradient, "last")()
  lost <- if (m > 1L) coalition_lost(out, tf$inv(z_root), m, k, ctrl) else NULL
  if (is.null(lost) && isTRUE(attr(sol, "merged"))) {
    lost <- "residents merged during the search"
  }
  if (!is.null(lost)) {
    converged <- FALSE
    z_root <- as.numeric(tf$fwd(out$traits))
  }

  if (any(on_edge)) {
    ## Candidates are clamped to the bounds, so a search region that does not
    ## contain a singularity ends up pinned against an edge -- which also makes
    ## the residual locally flat and the solver report failure. Report the
    ## informative cause, not the symptom.
    msg <- paste("Bounds do not include a singularity for trait(s)",
                 paste(labels[on_edge], collapse = ", "))
    if (edge_ok) warning(msg) else stop(msg)
  } else if (!is.null(lost)) {
    warning("community_solve_singularity lost the coalition: ", lost)
  } else if (!converged) {
    warning(sprintf("community_solve_singularity did not converge (%s: %s)",
                    solver, attr(sol, "message")))
  }

  root <- tf$inv(z_root)
  if (m == 1L) {
    names(root) <- trait_names
  } else {
    root <- trait_matrix(root, trait_names)
  }

  attr(out, "converged") <- converged
  attr(out, "solver") <- solver
  attr(out, "evaluations") <- attr(gradient, "evaluations")()
  attr(out, "singularity") <- root

  plant_log_assembler(sprintf(
    "Solved! %dD singularity for [%s] is [%s]",
    k, paste(trait_names, collapse = ", "),
    paste(signif(root, 6), collapse = ", ")))

  out
}

## Pairs of residents (rows i < j) that have merged: closer, in every trait,
## than the resident Jacobian's finite-difference step, so that no derivative
## across the residents can tell them apart. At coincidence the community is
## degenerate -- two copies of one strategy, its split between them undefined
## -- and a model's equilibrium may not even be finite there. `reach` widens
## the test to that many steps (2: a stencil point of one resident can land
## within a step of the other).
residents_merged <- function(x, m, k, ctrl, reach = 1) {
  none <- matrix(integer(0), 0L, 2L)
  if (m < 2L) return(none)
  xm <- matrix(x, m, k)
  h <- reach * matrix(util_fd_step(x, ctrl$d_second, ctrl$eps_second), m, k)
  pairs <- which(upper.tri(diag(m)), arr.ind = TRUE)
  close <- apply(pairs, 1, function(p) all(abs(xm[p[1], ] - xm[p[2], ]) < pmax(h[p[1], ], h[p[2], ])))
  pairs[close, , drop = FALSE]
}

## Why a solved set of m residents is not a coalition of m, or NULL if it is:
## two residents that have merged, or one without a positive equilibrium
## density.
coalition_lost <- function(community, x, m, k, ctrl) {
  merged <- residents_merged(x, m, k, ctrl)
  if (nrow(merged) > 0L) {
    return(paste(sprintf("residents %d and %d have merged", merged[, 1], merged[, 2]), collapse = "; "))
  }
  n <- as.numeric(community$birth_rate)
  dead <- which(!is.finite(n) | n <= 0)
  if (length(dead) > 0L) {
    return(sprintf("resident%s %s ha%s no positive equilibrium density",
                   if (length(dead) > 1L) "s" else "", paste(dead, collapse = ", "),
                   if (length(dead) > 1L) "ve" else "s"))
  }
  NULL
}

##' Classify a singular strategy or coalition.
##'
##' Finding a singular strategy is only half the question: the interesting part
##' is what happens near it. Two independent second-order conditions decide
##' that (Geritz et al. 1998; Leimar 2009):
##'
##' \describe{
##'   \item{Evolutionary stability}{the curvature of invasion fitness in the
##'     \emph{mutant} direction, with the residents held where they are.
##'     In one trait this is the scalar
##'     \eqn{\partial^2 s / \partial y^2}; in \code{k} traits it is the
##'     \code{k x k} Hessian \eqn{H}. A resident is uninvadable by nearby
##'     mutants when \eqn{H} at its traits is negative definite, i.e. every
##'     eigenvalue is negative; a coalition is an ESS when every resident is.}
##'   \item{Convergence stability}{whether selection carries nearby residents
##'     towards the point. This is the derivative of the selection gradients
##'     with respect to the \emph{resident} traits: the scalar \eqn{dG/dx} in
##'     one trait, the Jacobian \eqn{J} in \code{k}, and for \code{m} residents
##'     the \code{mk x mk} Jacobian of every resident's gradient with respect to
##'     every resident's traits, through the environment they share. Under
##'     the canonical equation with an isotropic mutational covariance each
##'     resident moves at its own speed \eqn{w_i} along its gradient, so the
##'     point is convergence stable when the eigenvalues of \eqn{D J}, with
##'     \eqn{D} the speeds on the diagonal, have negative real parts. For one
##'     resident the speed only rescales \eqn{J}; for a coalition it does not,
##'     and a coalition stable at equal speeds can be unstable at the speeds
##'     its densities give (the canonical equation's speed is proportional to
##'     the resident's equilibrium density). A negative-definite symmetric
##'     part \eqn{(J + J^T)/2} gives \emph{strong} convergence stability,
##'     which holds for any mutational covariance matrix and any (positive)
##'     speeds of the residents.}
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
##' For a coalition, "CSS" is an evolutionarily stable coalition that selection
##' reaches, and "branching point" one that selection reaches but where at
##' least one resident sits at a fitness minimum and will branch again
##' (\code{resident_evolutionarily_stable} says which).
##'
##' A coalition is first checked for protected coexistence: each resident,
##' made rare, must be able to invade the community the others settle at
##' without it. Positive densities alone do not make a coalition, and one a
##' resident lost by chance cannot return to is not an evolutionary outcome:
##' it is classified \code{"unprotected"} whatever its second-order
##' conditions. For two Lotka--Volterra residents with positive densities,
##' protection is exactly the stability of their equilibrium (priority
##' effects give positive densities at an unstable one). For three or more it
##' is a stronger condition than stability: a stable coalition can have a
##' resident that, once lost, stays lost, the others settling where it cannot
##' invade. The test is made on invasion fitness, so it means the same
##' whatever the model's time step.
##'
##' The community the others settle at is their equilibrium with any of them
##' that has no positive density there removed and the rest re-solved (a
##' model's own equilibrium can be infeasible, with negative densities, once a
##' resident is left out). Where that is not determined --- an equilibrium
##' solve fails, no resident persists, or a removed resident could invade the
##' rest back --- that resident's \code{invasion_fitness} is \code{NA};
##' \code{protected_coexistence} is then \code{NA} unless another resident
##' already fails, with a warning, and the classification is decided by the
##' second-order conditions.
##'
##' The full eigen-decompositions are returned, not just the verdict: when a
##' resident is not an ESS the leading eigenvector of its \eqn{H} is the
##' direction in trait space along which it disruptively splits, which is itself
##' the scientific result in a multi-trait problem.
##'
##' Cost: each Hessian (\code{\link{community_fitness_hessian}}) needs
##' \code{1 + 4k^2} invasion-fitness evaluations, all made in a single
##' vectorised call against the cached resident environment, so it is cheap. The
##' Jacobian (\code{\link{community_selection_gradient_jacobian}}) needs
##' \code{2mk} \emph{resident} evaluations, each a full demographic equilibrium
##' solve, so it dominates. Under a \code{future::plan()} with several workers
##' these solves, and those of the coexistence test, run in parallel, with
##' the same answer. Protected coexistence costs \code{m} more, one
##' per resident left out (and one more for each re-solve after removing a
##' resident without a positive density). Derivative settings come from
##' \code{\link{derivative_control}}.
##'
##' A community marked as not converged (\code{attr(., "converged")}
##' \code{FALSE}, set by \code{community_solve_singularity} or by the
##' equilibrium solve) is classified with a warning: the verdicts then describe
##' the point the solve stopped at. So is one whose residents are visibly away
##' from a singularity, a Newton step on their selection gradients moving some
##' trait by more than a thousandth of the bounds' width on the trait scale
##' (\code{newton_reach}): the second-order conditions describe a singularity,
##' and away from one they describe nothing in particular. A coalition with two
##' residents within a thousandth of the bounds' width of each other also
##' warns: at a pitchfork the coalition shrinks to a point, and a near-merged
##' pair can pass as a root whatever its verdict.
##'
##' @title Classify a singular strategy or coalition (1-D and N-D)
##' @param community A \code{community} whose residents are at (or very near) a
##' singular strategy or coalition --- typically the output of
##' \code{\link{community_solve_singularity}}.
##' @param birth_rate Birth rates to start each resident equilibrium solve from;
##' see \code{\link{community_solve_singularity}}. Defaults to the residents'
##' own equilibrium birth rates, which is normally what you want.
##' @param tol Magnitude below which an eigenvalue counts as zero, making the
##' classification degenerate rather than forcing a verdict.
##' @param speeds The residents' relative speeds of evolution, which decide
##' convergence stability for a coalition: \code{"density"} (default) their
##' equilibrium densities, as in \code{\link{community_canonical_equation}}
##' with \code{canonical_control(density = TRUE)}; \code{"equal"}; or a
##' positive vector with one value per resident. Ignored for one resident.
##' @param invasion_tol The invasion fitness when rare each resident of a
##' coalition must exceed for protected coexistence. It is in units of
##' invasion fitness, not of an eigenvalue, so worth setting apart from
##' \code{tol} when the model's fitness carries numerical noise (on plant,
##' its SCM's error).
##' @return An object of class \code{singularity_classification}: a list with
##' the traits and selection gradient at the point, \code{hessian} and
##' \code{jacobian} with their eigen-decompositions
##' (\code{hessian_eigen}, \code{jacobian_eigen}, \code{jacobian_symmetric_eigen}),
##' \code{speeds} (normalised to mean one) with \code{jacobian_weighted_eigen},
##' the eigen-decomposition of \eqn{D J} on which convergence stability is
##' decided (the same as \code{jacobian_eigen} for one resident or equal
##' speeds), \code{invasion_fitness} (each resident's when rare against the
##' others; \code{NULL} for one resident) with the verdict
##' \code{protected_coexistence} (\code{NA} where undetermined),
##' \code{newton_reach} (how far the point is from a singularity, as a
##' fraction of the bounds' width), the logical verdicts
##' \code{evolutionarily_stable} (with \code{resident_evolutionarily_stable},
##' one per resident), \code{convergence_stable} and
##' \code{strongly_convergence_stable}, the \code{branching_direction} where
##' the point is invadable, the \code{classification} (the four-way verdict,
##' \code{"degenerate"}, or \code{"unprotected"} for a coalition without
##' protected coexistence), and \code{evaluations} (the equilibrium solves of
##' the Jacobian and the coexistence test). For one resident \code{traits}, \code{selection_gradient} and
##' \code{branching_direction} are vectors named by trait and \code{hessian} a
##' matrix; for a coalition they are residents-by-traits matrices (a row of
##' \code{NA} in \code{branching_direction} for a resident that is not at a
##' fitness minimum) and \code{hessian} and \code{hessian_eigen} are lists with
##' one entry per resident.
##' @author Daniel Falster
##' @export
community_classify_singularity <- function(community, birth_rate = NULL,
                                           tol = 1e-8, speeds = "density",
                                           invasion_tol = tol) {

  trait_names <- community$trait_names
  k <- length(trait_names)
  m <- nrow(community$traits)

  if (m < 1L) {
    stop("community_classify_singularity needs at least one resident (the ",
         "singular strategy or coalition)")
  }
  x <- community$traits
  unconverged <- isFALSE(attr(community, "converged"))
  if (unconverged) {
    warning("community_classify_singularity: the community is marked as not converged ",
            "(its singularity or its equilibrium solve), so this classifies the point ",
            "it stopped at, which need not be a singularity")
  }

  if (is.null(community$fitness_function)) {
    community <- community_demography(community)
  }

  plant_log_assembler(sprintf(
    "Classifying %dD singularity at [%s] = [%s]",
    k, paste(trait_names, collapse = ", "),
    paste(signif(x, 6), collapse = ", ")))

  ## --- evolutionary stability: each resident's curvature in the mutant
  ## direction ----------------------------------------------------------------
  H <- lapply(seq_len(m), function(i) community_fitness_hessian(community, x[i, , drop = FALSE]))
  H_eigen <- lapply(H, function(h) eigen((h + t(h)) / 2, symmetric = TRUE))

  ## --- convergence stability: how the gradients respond to the residents ----
  J <- community_selection_gradient_jacobian(community, birth_rate = birth_rate)
  g0 <- as.numeric(attr(J, "selection_gradient"))
  evaluations <- attr(J, "evaluations")
  attr(J, "selection_gradient") <- NULL
  attr(J, "evaluations") <- NULL
  J_eigen <- eigen(J)
  J_sym <- (J + t(J)) / 2
  J_sym_eigen <- eigen(J_sym, symmetric = TRUE)
  ## each resident's rows scaled by its speed (stacked trait by trait,
  ## residents within each)
  w <- classification_speeds(speeds, community, m)
  J_w_eigen <- if (all(w == 1)) J_eigen else eigen(rep(w, times = k) * J)

  ## a community not marked as unconverged can still be away from a
  ## singularity (built by hand, or solved to a loose tolerance)
  reach <- singularity_newton_reach(community, x, g0, J)
  if (!unconverged && isTRUE(reach > singularity_reach_tol)) {
    warning(sprintf(paste0(
      "community_classify_singularity: the residents are not at a singular point; ",
      "a Newton step on the selection gradients moves a trait by %s of the bounds' width"),
      signif(reach, 2)))
  }
  ## at a pitchfork a coalition shrinks to a point, and every derivative with
  ## it, so a near-merged pair can pass as a root with any verdict
  close <- coalition_near_merged(community, x)
  if (nrow(close) > 0L) {
    warning(sprintf(paste0(
      "community_classify_singularity: %s within %s of the bounds' width of each other; ",
      "a coalition this close to merging is at or near a bifurcation, where the verdict ",
      "is not reliable"),
      paste(sprintf("residents %d and %d are", close[, 1], close[, 2]), collapse = "; "),
      singularity_reach_tol))
  }

  ## --- protected coexistence: can each resident invade the others? ---------
  invasion <- NULL
  protected <- TRUE
  if (m > 1L) {
    invasion <- coalition_invasion_fitness(community, birth_rate = birth_rate)
    evaluations <- evaluations + attr(invasion, "evaluations")
    reasons <- attr(invasion, "reasons")
    attr(invasion, "evaluations") <- NULL
    attr(invasion, "reasons") <- NULL
    known <- !is.na(invasion)
    protected <- if (any(known & invasion <= invasion_tol)) FALSE else if (all(known)) TRUE else NA
    if (is.na(protected)) {
      i <- which(!known)
      warning(sprintf(paste0(
        "community_classify_singularity: protected coexistence is undetermined; ",
        "the invasion fitness of %s could not be found (%s)"),
        paste("resident", i, collapse = ", "),
        paste(unique(unlist(reasons[i])), collapse = "; ")))
    }
  }

  ## --- verdicts -------------------------------------------------------------
  ev_H <- lapply(H_eigen, `[[`, "values")
  ev_J <- Re(J_w_eigen$values)
  ev_Js <- J_sym_eigen$values

  ess_i <- vapply(ev_H, function(v) all(v < -tol), logical(1))
  ess <- all(ess_i)
  cs <- all(ev_J < -tol)
  scs <- all(ev_Js < -tol)

  degenerate <- any(abs(unlist(ev_H)) <= tol) || any(abs(ev_J) <= tol)

  classification <-
    if (isFALSE(protected)) {
      "unprotected"
    } else if (degenerate) {
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

  ## Where a resident is invadable, the leading eigenvector of its Hessian is
  ## the trait-space direction along which fitness curves upward -- the
  ## direction it splits along. eigen() fixes the sign arbitrarily; make it
  ## reproducible by pointing the largest component positive (the direction is
  ## an axis, not an arrow).
  direction <- lapply(seq_len(m), function(i) {
    e <- H_eigen[[i]]
    if (ess_i[i] || max(e$values) <= tol) return(rep(NA_real_, k))
    v <- e$vectors[, which.max(e$values)]
    v * sign(v[which.max(abs(v))])
  })
  branching_direction <- trait_matrix(do.call(rbind, direction), trait_names)

  if (m == 1L) {
    traits <- stats::setNames(as.numeric(x[1, ]), trait_names)
    selection_gradient <- stats::setNames(g0, trait_names)
    branching_direction <- if (anyNA(branching_direction)) NULL else
      stats::setNames(as.numeric(branching_direction), trait_names)
    H <- H[[1L]]
    H_eigen <- H_eigen[[1L]]
  } else {
    traits <- trait_matrix(as.numeric(x), trait_names)
    selection_gradient <- trait_matrix(g0, trait_names)
    if (all(is.na(branching_direction))) branching_direction <- NULL
  }

  ret <- list(
    trait_names = trait_names,
    traits = traits,
    selection_gradient = selection_gradient,
    hessian = H,
    hessian_eigen = H_eigen,
    jacobian = J,
    jacobian_eigen = J_eigen,
    jacobian_symmetric_eigen = J_sym_eigen,
    speeds = w,
    jacobian_weighted_eigen = J_w_eigen,
    invasion_fitness = invasion,
    protected_coexistence = protected,
    newton_reach = reach,
    evolutionarily_stable = ess,
    resident_evolutionarily_stable = ess_i,
    convergence_stable = cs,
    strongly_convergence_stable = scs,
    branching_direction = branching_direction,
    degenerate = degenerate,
    classification = classification,
    evaluations = evaluations,
    tol = tol
  )
  class(ret) <- "singularity_classification"

  plant_log_assembler(sprintf("Classified! Singularity at [%s] is a %s",
                              paste(signif(x, 6), collapse = ", "),
                              classification))

  ret
}

## The invasion fitness of `invader` (a 1 x k trait matrix), rare, in the
## community the `residents` settle at: their equilibrium, solved from
## `birth_rate`, with any resident that has no positive density there removed
## and the rest re-solved. A model's own equilibrium (the "model" solver) can
## return negative densities, which describe no community at all, and an
## invader's fitness against them can have either sign. The community left must
## be one the removed residents cannot invade back, or it is not where the
## residents settle; that, an equilibrium solve that fails, or no resident
## persisting, leaves the answer undetermined: NA, with the reason. The
## equilibrium solves are counted. `base` is a community without residents.
community_invasion_when_rare <- function(base, residents, birth_rate, invader) {
  trait_names <- base$trait_names
  evaluations <- 0L
  removed <- residents[0, , drop = FALSE]
  undetermined <- function(why) structure(NA_real_, evaluations = evaluations, reason = why)
  repeat {
    comm <- tryCatch(
      base |>
        community_add(trait_matrix(residents, trait_names), birth_rate = birth_rate) |>
        community_demography(),
      error = function(e) NULL)
    evaluations <- evaluations + 1L
    if (is.null(comm) || !isTRUE(attr(comm, "converged"))) {
      return(undetermined("an equilibrium solve failed"))
    }
    n <- as.numeric(comm$birth_rate)
    dead <- !is.finite(n) | n <= 0
    if (!any(dead)) break
    if (all(dead)) {
      return(undetermined("no resident persists"))
    }
    removed <- rbind(removed, residents[dead, , drop = FALSE])
    residents <- residents[!dead, , drop = FALSE]
    birth_rate <- n[!dead]
  }
  f <- community_fitness_function(comm)
  if (nrow(removed) > 0L && any(as.numeric(f(removed)) > 0)) {
    return(undetermined("a resident removed for want of a positive density could invade back"))
  }
  structure(as.numeric(f(invader)), evaluations = evaluations, reason = NULL)
}

## One invasion test, a list of residents, birth_rate and invader, against
## the community without residents `base`; the tests of a coalition or a
## branching are independent, so they go through the parallel map.
invasion_when_rare_test <- function(test, base) {
  community_invasion_when_rare(base, test$residents, test$birth_rate, test$invader)
}

## Each resident's invasion fitness when rare in the community the others
## settle at without it (community_invasion_when_rare), with the equilibrium
## solves it cost and, for an undetermined one (NA), the reason. Every one
## positive is protected coexistence: no resident lost by chance stays lost.
coalition_invasion_fitness <- function(community, birth_rate = NULL) {
  m <- nrow(community$traits)
  base <- community_clear_residents(community)
  n <- if (is.null(birth_rate)) as.numeric(community$birth_rate) else as.numeric(birth_rate)
  x <- community$traits
  tests <- lapply(seq_len(m), function(i) {
    list(residents = x[-i, , drop = FALSE], birth_rate = n[-i], invader = x[i, , drop = FALSE])
  })
  out <- regnans_map(tests, invasion_when_rare_test, base = base)
  s <- vapply(out, as.numeric, numeric(1))
  reasons <- lapply(out, function(o) {
    if (is.null(attr(o, "reason")) && is.na(o)) "the model's invasion fitness is not a number" else attr(o, "reason")
  })
  attr(s, "evaluations") <- sum(vapply(out, attr, integer(1), "evaluations"))
  attr(s, "reasons") <- reasons
  s
}

## A Newton step from a singularity found to the solver's tolerance is many
## orders of magnitude below this fraction of the bounds' width; one above it
## means the point classified is not a singularity.
singularity_reach_tol <- 1e-3

## The largest move, as a fraction of the bounds' width on the trait scale, of
## any resident's trait under one Newton step on the selection gradients from
## here: how far the point is from a singularity, by the Jacobian's own measure
## (NA where the Jacobian is singular or the bounds infinite).
singularity_newton_reach <- function(community, x, g, J) {
  dx <- tryCatch(solve(J, g), error = function(e) NULL)
  if (is.null(dx)) return(NA_real_)
  tf <- community_trait_transform(community)
  bounds <- singularity_bounds(community$bounds, community$trait_names)
  width <- rep(tf$fwd(bounds[, 2]) - tf$fwd(bounds[, 1]), each = nrow(x))
  xs <- as.numeric(x)
  dz <- if (identical(tf$scale, "log")) dx / xs else dx
  r <- abs(dz) / width
  if (!all(is.finite(r))) return(NA_real_)
  max(r)
}

## Pairs of residents (rows i < j) closer than singularity_reach_tol of the
## bounds' width, on the trait scale, in every trait.
coalition_near_merged <- function(community, x) {
  m <- nrow(x)
  none <- matrix(integer(0), 0L, 2L)
  if (m < 2L) return(none)
  tf <- community_trait_transform(community)
  bounds <- singularity_bounds(community$bounds, community$trait_names)
  width <- tf$fwd(bounds[, 2]) - tf$fwd(bounds[, 1])
  if (!all(is.finite(width))) return(none)
  z <- matrix(tf$fwd(as.numeric(x)), m)
  pairs <- which(upper.tri(diag(m)), arr.ind = TRUE)
  close <- apply(pairs, 1, function(p) all(abs(z[p[1], ] - z[p[2], ]) < singularity_reach_tol * width))
  pairs[close, , drop = FALSE]
}

## The residents' relative speeds for the convergence-stability verdict,
## normalised to mean one (so one resident, or equal speeds, is all ones).
classification_speeds <- function(speeds, community, m) {
  if (is.character(speeds)) {
    speeds <- match.arg(speeds, c("density", "equal"))
    if (m == 1L || speeds == "equal") return(rep(1, m))
    speeds <- as.numeric(community$birth_rate)
    if (length(speeds) != m || any(!is.finite(speeds) | speeds <= 0)) {
      stop("speeds = \"density\" needs a positive equilibrium density for every resident")
    }
  }
  if (!is.numeric(speeds) || length(speeds) != m || any(!is.finite(speeds) | speeds <= 0)) {
    stop("speeds must be \"density\", \"equal\" or ", m, " positive numbers")
  }
  if (m == 1L) return(1)
  speeds / mean(speeds)
}

##' @param x A \code{singularity_classification} object.
##' @param ... Ignored.
##' @rdname community_classify_singularity
##' @export
print.singularity_classification <- function(x, ...) {
  m <- NROW(x$traits)
  fmt_res <- function(v) paste(sprintf("%s = %s", x$trait_names, signif(v, 6)), collapse = ", ")
  cat(sprintf("<singularity_classification: %s%s>\n", x$classification,
              if (m > 1L) sprintf(" (coalition of %d)", m) else ""))
  if (m == 1L) {
    cat(sprintf("  traits: %s\n", fmt_res(x$traits)))
    cat(sprintf("  selection gradient: %s\n",
                paste(signif(x$selection_gradient, 3), collapse = ", ")))
    cat(sprintf("  evolutionarily stable: %s (Hessian eigenvalues %s)\n",
                x$evolutionarily_stable,
                paste(signif(x$hessian_eigen$values, 4), collapse = ", ")))
  } else {
    for (i in seq_len(m)) {
      cat(sprintf("  resident %d: %s; evolutionarily stable: %s (Hessian eigenvalues %s)\n",
                  i, fmt_res(x$traits[i, ]), x$resident_evolutionarily_stable[i],
                  paste(signif(x$hessian_eigen[[i]]$values, 4), collapse = ", ")))
    }
    cat(sprintf("  largest selection gradient: %s\n",
                signif(max(abs(x$selection_gradient)), 3)))
    cat(sprintf("  protected coexistence: %s (invasion fitness when rare %s)\n",
                x$protected_coexistence, paste(signif(x$invasion_fitness, 4), collapse = ", ")))
  }
  weighted <- max(abs(x$speeds - 1)) > 1e-6
  cat(sprintf("  convergence stable:    %s (%s eigenvalues %s)\n",
              x$convergence_stable,
              if (weighted) "speed-weighted Jacobian" else "Jacobian",
              paste(signif(Re(x$jacobian_weighted_eigen$values), 4), collapse = ", ")))
  if (weighted) {
    cat(sprintf("  resident speeds:       %s\n", paste(signif(x$speeds, 4), collapse = ", ")))
  }
  cat(sprintf("  strongly conv. stable: %s\n", x$strongly_convergence_stable))
  if (!is.null(x$branching_direction)) {
    if (m == 1L) {
      cat(sprintf("  branching direction:   %s\n",
                  paste(sprintf("%s = %s", x$trait_names,
                                signif(x$branching_direction, 4)),
                        collapse = ", ")))
    } else {
      for (i in which(!is.na(x$branching_direction[, 1]))) {
        cat(sprintf("  branching direction of resident %d: %s\n", i,
                    paste(sprintf("%s = %s", x$trait_names,
                                  signif(x$branching_direction[i, ], 4)),
                          collapse = ", ")))
      }
    }
  }
  invisible(x)
}
