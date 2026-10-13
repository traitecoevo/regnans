# Coexistence boundaries followed through a plane of two parameters (#60,
# stages 4 and 5).
#
#   community_continue_boundary()  pseudo-arclength continuation of the
#                                  augmented system of community_solve_boundary()
#                                  with two parameters free, test functions
#                                  monitored along it
#   boundary_switch()              the start of the branch an event seeds
#   boundary_diagram()             branches collected for plotting
#
# With two parameters free the augmented system has one more unknown than
# equations, and its solutions form a curve. Each step predicts along the
# curve's tangent (the null vector of the system's Jacobian) and corrects back
# onto it by Newton on the system bordered by the tangent, so the step length
# along the curve rather than any one coordinate is fixed and the curve can turn
# back in either parameter. A step whose tangent turns sharply is rejected: the
# corrector has jumped to another branch, which is how continuation fails
# silently. Test functions are monitored at every accepted point, each positive
# while the boundary means what it says; where one changes sign the step is
# searched by false position to locate the event, the branch ends there, and boundary_switch()
# says which branch starts.

##' Follow a coexistence boundary through a plane of two parameters.
##'
##' Starting from a point on a boundary (one from
##' \code{\link{community_solve_boundary}}, or a bracket from
##' \code{\link{boundary_bisect}}), follows the curve that the augmented system
##' traces when two parameters are free: \eqn{nk + k + 1} equations in
##' \eqn{nk + k + 2} unknowns for an invasion boundary, \eqn{nk + 1} in
##' \eqn{nk + 2} for a branching one. Each step predicts along the tangent and
##' corrects by Newton on the system bordered by it (pseudo-arclength
##' continuation), with the Jacobian from the derivative layer
##' (\code{2nk + 4} equilibrium solves, through the parallel map). Steps adapt
##' to the corrector's iterations between \code{ds_min} and \code{ds_max}; a
##' step that does not converge, or whose tangent turns so that its dot product
##' with the last falls below \code{tangent_tol}, is halved and retried.
##'
##' At every accepted point test functions are evaluated, each positive while
##' the boundary describes the coexistence it was found for:
##' \describe{
##'   \item{\code{density}}{the rarest resident's share of the density, less
##'     \code{density_tol}: a resident going extinct, where the boundary meets
##'     the one of a community without it.}
##'   \item{\code{resident_curvature}}{minus the largest eigenvalue of the
##'     residents' Hessians (all residents of an invasion boundary, the others
##'     of a branching one): a resident ceasing to be an ESS, where the boundary
##'     meets that resident's branching locus.}
##'   \item{\code{peak_curvature}}{(invasion) minus the largest eigenvalue of
##'     the Hessian at \eqn{y^*}: the marginal peak annihilating with a valley.}
##'   \item{\code{separation}}{the smallest distance between residents (and
##'     \eqn{y^*}), as a fraction of the bounds' width on the trait scale, less
##'     \code{merge_tol}: residents merging, as at a pitchfork.}
##'   \item{\code{other_peak}}{(one trait) \code{fitness_tol} less the
##'     highest fitness on the landscape farther than twice \code{merge_tol}
##'     from the residents and \eqn{y^*} (searched on \code{n_landscape}
##'     points, refined): another peak becoming marginal first, so that the
##'     boundary is now that peak's; on a branching locus, also the
##'     landscape beside the resident turning upward, where the pitchfork
##'     turns subcritical and the locus stops being a coexistence boundary.}
##'   \item{\code{bordered}}{the determinant of the Jacobian bordered by the
##'     tangent, which changes sign where another branch crosses this one.}
##' }
##' Where one changes sign the step is searched by false position, the trial
##' points corrected without the Jacobian that only \code{bordered} needs, until
##' the event is bracketed to \code{event_tol} in arclength, the located point
##' ends the branch, and the event is recorded (unless it lies beyond
##' \code{box} or the trait bounds: the branch has left them first, and ends
##' there); \code{\link{boundary_switch}} gives the start of
##' the branch it seeds. A branch also ends where it closes on itself, leaves
##' \code{box}, takes a resident or \eqn{y^*} out of the trait bounds, can no
##' longer step (\code{ds} below \code{ds_min}) or reaches \code{max_steps}.
##'
##' The two directions from the start are traced as two jobs of the parallel
##' map, with the same answer under any \code{future::plan()}. Afterwards every
##' \code{verify_every}-th point and both ends (the point before an event,
##' where two boundaries meet and labels cannot tell them apart) are verified
##' cold, also through
##' the map: the residual re-solved from the model's initial densities, with no
##' warm start, and the labels (\code{\link{community_label}}) either side of
##' the curve: the boundary's \code{n} species, stable, on one side and
##' \code{n + 1} on the other, or a count censored at \code{n} or more (as
##' where DD99's box fills beyond its branching locus). A point that
##' fails is reported in \code{verification} and warned of: systematic drift,
##' or a boundary that has stopped being one, shows here.
##'
##' @title Continue a coexistence boundary in two parameters
##' @inheritParams community_solve_boundary
##' @param free The two parameters that vary, by name or position in \code{p};
##' the others are held. The start is first corrected onto the curve with the
##' second held.
##' @param box Optional: a two-row matrix of lower and upper limits (columns
##' \code{lower}, \code{upper}) for the free parameters, in the order of
##' \code{free}; the branch ends where it leaves it.
##' @param direction \code{"both"} (from the start in both directions),
##' \code{"forward"} (the second free parameter increasing at the start) or
##' \code{"backward"}.
##' @param control See \code{\link{boundary_control}}.
##' @return An object of class \code{regnans_branch}: a list with
##' \code{points} (a data frame: the arclength \code{s}, the free parameters,
##' the residents' traits, \eqn{y^*}, the densities, the test functions and the
##' corrector's iterations, one row per point in order along the curve), the
##' unknowns and unit tangents in \code{u} and \code{tangent} (matrices, one row
##' per point), \code{start} (the row of the start), \code{events} (a data
##' frame: \code{end}, \code{type}, \code{s}, the free parameters, \code{detail}),
##' \code{termination} (why each direction ended), \code{verification},
##' \code{rejected} steps, \code{evaluations} (equilibrium solves) and what
##' \code{\link{boundary_switch}} needs to start a branch from an event.
##' @examples
##' \donttest{
##' # GK98's three patches, -d, 0 and d, the middle one kappa times the others;
##' # the map seeds one resident, where every label starts
##' gk <- function(p) {
##'   community_start(bounds(x = c(-4, 4)), trait_scale = "linear",
##'                   harness = harness_gk98(mu = c(-p[1], 0, p[1]), K = c(1, p[2], 1))) |>
##'     community_add(trait_matrix(0.5, "x"), birth_rate = 1)
##' }
##' pair <- community_start(bounds(x = c(-4, 4)), trait_scale = "linear",
##'                         harness = harness_gk98(mu = c(-1.93, 0, 1.93))) |>
##'   community_add(trait_matrix(c(-1.27, 1.27), "x"), birth_rate = c(1.5, 1.5))
##' # where a third species can invade the pair, as the middle patch grows
##' branch <- community_continue_boundary(pair, gk, y = 0, p = c(1.93, 1),
##'                                       box = rbind(c(1, 3), c(0.5, 3)),
##'                                       direction = "forward")
##' branch
##' plot(branch)
##' }
##' @author Daniel Falster
##' @export
community_continue_boundary <- function(community, parameter, y = NULL,
                                        type = c("invasion", "branching"), free = 1:2,
                                        which = 1L, p = NULL, box = NULL,
                                        direction = c("both", "forward", "backward"),
                                        control = boundary_control()) {
  type <- match.arg(type)
  direction <- match.arg(direction)
  control <- boundary_control(control)
  p <- check_parameter_values(parameter, p, community)
  free <- if (is.character(free)) match(free, names(p)) else as.integer(free)
  if (length(free) != 2L || anyNA(free) || anyDuplicated(free) || any(free < 1L | free > length(p))) {
    stop("free must name two distinct parameters of ", paste(names(p), collapse = ", "))
  }
  box <- boundary_check_box(box, names(p)[free])
  trait_names <- community$trait_names

  ## onto the curve, the second free parameter held
  seed <- community_solve_boundary(community, parameter, y = y, type = type, free = free[1L],
                                   which = which, p = p, tol = control$corrector_tol)
  b <- attr(seed, "boundary")
  if (!isTRUE(b$converged)) {
    stop("community_continue_boundary could not put the start on the boundary: ", b$message)
  }
  start <- community_clear_residents(parameter_community(parameter, b$p, trait_names)) |>
    community_add(trait_matrix(seed$traits, trait_names), birth_rate = as.numeric(seed$birth_rate))
  residual <- boundary_residual_fn(start, parameter, b$p, free, type, which)
  u0 <- c(as.numeric(seed$traits), b$y, b$p[free])
  J0 <- community_boundary_jacobian(residual, u0)
  evaluations <- b$evaluations + attr(residual, "evaluations")() + attr(J0, "evaluations")
  t0 <- boundary_tangent(unname(J0))
  ## forward is the second free parameter increasing, or failing that the first
  lead <- if (abs(t0[length(t0)]) > 1e-12) t0[length(t0)] else t0[length(t0) - 1L]
  if (lead < 0) t0 <- -t0

  signs <- switch(direction, both = c(1, -1), forward = 1, backward = -1)
  traces <- regnans_map(as.list(signs), boundary_trace, start = start, parameter = parameter,
                        p = b$p, free = free, type = type, which = which, u0 = u0, t0 = t0,
                        J0 = unname(J0), box = box, control = control)
  failed <- vapply(traces, inherits, logical(1), "error")
  if (any(failed)) {
    stop("community_continue_boundary failed tracing the branch: ", conditionMessage(traces[[which(failed)[1L]]]))
  }
  names(traces) <- c("forward", "backward")[match(signs, c(1, -1))]
  ## a closed curve is the whole branch in one direction
  if (direction == "both" && identical(traces$forward$termination, "closed")) {
    traces$backward <- NULL
  }
  branch <- boundary_branch(traces, info = attr(residual, "info"), seed = seed, start = start,
                            box = box, control = control)
  branch$evaluations <- branch$evaluations + evaluations
  branch$verification <- boundary_verify(branch, start, parameter, control)
  branch$evaluations <- branch$evaluations + sum(branch$verification$evaluations)
  bad <- !branch$verification$ok
  if (any(bad)) {
    warning(sprintf(paste0("community_continue_boundary: %d of %d verified points failed (at s = %s): ",
                           "the branch has drifted, or stopped being the boundary it was found as"),
                    sum(bad), length(bad), paste(signif(branch$verification$s[bad], 4), collapse = ", ")))
  }
  branch
}

## A box for the free parameters as a 2 x 2 matrix, rows in the order of
## free; NULL for none.
boundary_check_box <- function(box, names) {
  if (is.null(box)) return(NULL)
  if (!is.matrix(box) || !identical(dim(box), c(2L, 2L)) || !is.numeric(box) || anyNA(box) ||
      any(box[, 1] >= box[, 2])) {
    stop("box must be a 2 x 2 matrix of lower and upper limits, one row per free parameter")
  }
  if (!is.null(rownames(box)) && !identical(rownames(box), names)) {
    if (!setequal(rownames(box), names)) {
      stop("box's rows must be the free parameters, ", paste(names, collapse = ", "))
    }
    box <- box[names, , drop = FALSE]
  }
  dimnames(box) <- list(names, c("lower", "upper"))
  box
}

## The unit tangent of a curve whose Jacobian is J (one more column than
## rows): its null vector, oriented along `along` where given.
boundary_tangent <- function(J, along = NULL) {
  n <- ncol(J)
  q <- qr.Q(qr(t(J)), complete = TRUE)[, n]
  if (!is.null(along) && sum(q * along) < 0) q <- -q
  q
}

## One direction of a branch, on a worker: returns numbers only (the points,
## the event and why it ended), or the error.
boundary_trace <- function(sign, start, parameter, p, free, type, which, u0, t0, J0, box, control) {
  tryCatch({
    residual <- boundary_residual_fn(start, parameter, p, free, type, which)
    info <- attr(residual, "info")
    jacobian_solves <- 0L
    jacobian <- function(u) {
      J <- community_boundary_jacobian(residual, u)
      jacobian_solves <<- jacobian_solves + attr(J, "evaluations")
      unname(J)
    }
    u <- u0
    t <- sign * t0
    J <- J0
    tests <- boundary_test_functions(residual, u, J, t, info, control)
    s <- 0
    ds <- control$ds
    points <- list(boundary_trace_point(residual, u, t, s, tests, 0L))
    event <- NULL
    termination <- NULL
    rejected <- 0L
    step <- 0L
    while (is.null(termination)) {
      if (step >= control$max_steps) {
        termination <- "max_steps"
        break
      }
      out <- boundary_corrector(residual, jacobian, u, t, J, ds, control)
      if (!is.null(out$u) && sum(out$t * t) < control$tangent_tol) {
        out <- list(u = NULL, reason = sprintf("the tangent turned (dot product %.3g)", sum(out$t * t)))
      }
      if (is.null(out$u)) {
        rejected <- rejected + 1L
        ds <- ds / 2
        if (ds < control$ds_min) {
          termination <- paste("the step fell below ds_min:", out$reason)
        }
        next
      }
      step <- step + 1L
      new_tests <- boundary_test_functions(residual, out$u, out$J, out$t, info, control)
      out$tests <- new_tests
      crossed <- boundary_crossed(tests, new_tests)
      if (length(crossed) > 0L) {
        out <- boundary_locate_event(residual, jacobian, u, t, J, tests, crossed, ds, out, info, control)
        ## an event beyond the box or the trait bounds is not on this branch:
        ## the branch leaves them first
        if (is.null(boundary_outside(out$u, info, box, start))) {
          crossed <- boundary_crossed(tests, out$tests)
          s <- s + sqrt(sum((out$u - u)^2))
          points[[length(points) + 1L]] <- boundary_trace_point(residual, out$u, out$t, s, out$tests,
                                                                out$iterations)
          event <- list(type = crossed[1L], s = s, detail = attr(out$tests, "detail")[crossed])
          termination <- paste("event:", paste(crossed, collapse = ", "))
          break
        }
      }
      termination <- boundary_outside(out$u, info, box, start)
      if (!is.null(termination)) {
        break
      }
      s <- s + sqrt(sum((out$u - u)^2))
      u <- out$u
      t <- out$t
      J <- out$J
      tests <- new_tests
      points[[length(points) + 1L]] <- boundary_trace_point(residual, u, t, s, tests, out$iterations)
      if (step >= 3L && sqrt(sum((u - u0)^2)) < ds && sum(t * sign * t0) > 0) {
        termination <- "closed"
        break
      }
      ds <- if (out$iterations <= 2L) min(ds * 1.5, control$ds_max)
            else if (out$iterations >= 5L) max(ds / 2, control$ds_min)
            else ds
    }
    list(points = points, event = event, termination = termination, rejected = rejected,
         evaluations = attr(residual, "evaluations")() + jacobian_solves)
  }, error = function(e) e)
}

## A point as numbers: the unknowns, tangent, arclength, densities and test
## functions at the community the residual last solved.
boundary_trace_point <- function(residual, u, t, s, tests, iterations) {
  residual(u)
  list(u = u, t = t, s = s, birth_rate = as.numeric(attr(residual, "last")()$birth_rate),
       tests = as.numeric(tests), test_names = names(tests), detail = attr(tests, "detail"),
       iterations = as.integer(iterations))
}

## One pseudo-arclength step of length ds from u along t: Newton on the
## system bordered by the tangent, from the prediction, starting from the
## last point's Jacobian; then, with `tangent`, the Jacobian and tangent at
## the new point (2nk + 4 solves, which a trial point of an event's location
## does without).
boundary_corrector <- function(residual, jacobian, u, t, J, ds, control, tangent = TRUE) {
  u_pred <- u + ds * t
  bordered <- function(v) c(residual(v), sum(t * (v - u_pred)))
  sol <- tryCatch(util_nlsolve(u_pred, bordered, tol = control$corrector_tol, maxit = control$corrector_maxit,
                               solver = "newton", require_converged = FALSE, J0 = rbind(J, t),
                               jac = function(v) rbind(jacobian(v), t)),
                  error = function(e) e)
  if (inherits(sol, "error")) {
    return(list(u = NULL, reason = conditionMessage(sol)))
  }
  if (!isTRUE(attr(sol, "converged"))) {
    return(list(u = NULL, reason = attr(sol, "message")))
  }
  v <- as.numeric(sol)
  if (!tangent) {
    return(list(u = v, iterations = attr(sol, "iter")))
  }
  Jv <- jacobian(v)
  list(u = v, J = Jv, t = boundary_tangent(Jv, t), iterations = attr(sol, "iter"))
}

## The test functions at u, each positive while the boundary describes the
## coexistence it was found for, with what an event would need in
## attr(., "detail"): which resident, where another peak is. Without J
## (a trial point of an event's location) there is no `bordered`.
boundary_test_functions <- function(residual, u, J, t, info, control) {
  residual(u)
  comm <- attr(residual, "last")()
  v <- boundary_unknowns(u, info)
  m <- info$m
  k <- info$k
  n <- as.numeric(comm$birth_rate)
  share <- n / sum(n)
  values <- c(density = min(share) - control$density_tol)
  detail <- list(density = which.min(share))

  ess <- if (identical(info$type, "invasion")) seq_len(m) else setdiff(seq_len(m), info$which)
  if (length(ess) > 0L) {
    top <- vapply(ess, function(i) {
      H <- community_fitness_hessian(comm, comm$traits[i, , drop = FALSE])
      max(eigen(H, symmetric = TRUE, only.values = TRUE)$values)
    }, numeric(1))
    values["resident_curvature"] <- -max(top)
    detail$resident_curvature <- ess[which.max(top)]
  }
  if (identical(info$type, "invasion")) {
    H <- community_fitness_hessian(comm, matrix(v$y, nrow = 1L))
    values["peak_curvature"] <- -max(eigen(H, symmetric = TRUE, only.values = TRUE)$values)
  }

  tf <- community_trait_transform(comm)
  lb <- tf$fwd(comm$bounds)
  width <- lb[, 2] - lb[, 1]
  z <- rbind(matrix(tf$fwd(as.numeric(comm$traits)), m, k),
             if (identical(info$type, "invasion")) matrix(tf$fwd(v$y), 1L, k))
  if (nrow(z) >= 2L) {
    pairs <- which(upper.tri(diag(nrow(z))), arr.ind = TRUE)
    gap <- apply(pairs, 1, function(pr) max(abs(z[pr[1], ] - z[pr[2], ]) / width))
    values["separation"] <- min(gap) - control$merge_tol
    ## the collapse is one event: when the nearest two reach merge_tol, a pair
    ## straddling y* is twice that apart, so the clusters reach further
    detail$separation <- boundary_clusters(z, width, 2.5 * control$merge_tol)
  }

  if (k == 1L) {
    other <- boundary_other_peak(comm, z, width, tf, control)
    values["other_peak"] <- control$fitness_tol - other$fitness
    detail$other_peak <- other$at
  }
  if (!is.null(J)) {
    values["bordered"] <- det(rbind(J, t))
  }
  attr(values, "detail") <- detail
  values
}

## The points (residents and y*, on the trait scale) grouped by single
## linkage at `reach` of the bounds' width: an integer per point.
boundary_clusters <- function(z, width, reach) {
  n <- nrow(z)
  group <- seq_len(n)
  for (i in seq_len(n)) for (j in seq_len(n)) {
    if (i < j && max(abs(z[i, ] - z[j, ]) / width) < reach) {
      group[group == group[j]] <- group[i]
    }
  }
  match(group, unique(group))
}

## The highest fitness of a one-trait landscape away from the residents and
## y*: the maximum over the bounds outside 2 merge_tol of each of them, from
## a grid refined by optimise() around its best point. A maximum rather than
## the highest local peak, so it is continuous where a peak first appears:
## on the edge of a neighbourhood, where it rises past zero as the landscape
## beside a resident turns upward (a pitchfork turning subcritical), or
## inside the bounds.
boundary_other_peak <- function(comm, z, width, tf, control) {
  lb <- tf$fwd(comm$bounds)
  reach <- 2 * control$merge_tol * width
  away <- function(zz, slack = 0) vapply(zz, function(v) all(abs(v - z[, 1]) >= reach - slack), logical(1))
  grid <- seq(lb[1, 1], lb[1, 2], length.out = control$n_landscape)
  ## each neighbourhood's edges are kept unless another neighbourhood covers
  ## them; rounding can put an edge a hair inside the reach it was built at
  edges <- pmin(pmax(c(z[, 1] - reach, z[, 1] + reach), lb[1, 1]), lb[1, 2])
  grid <- sort(c(grid[away(grid)], edges[away(edges, 1e-9 * reach)]))
  if (length(grid) == 0L) {
    return(list(fitness = -Inf, at = NA_real_))
  }
  f <- function(zz) as.numeric(comm$fitness_function(trait_matrix(tf$inv(zz), comm$trait_names)))
  w <- f(grid)
  i <- which.max(w)
  best <- list(fitness = w[i], at = grid[i])
  ## refine only inside a stretch clear of every neighbourhood: removing a
  ## neighbourhood leaves its two edges adjacent on the grid
  if (i > 1L && i < length(grid) && !any(z[, 1] > grid[i - 1L] - reach & z[, 1] < grid[i + 1L] + reach)) {
    opt <- stats::optimise(f, grid[c(i - 1L, i + 1L)], maximum = TRUE, tol = 1e-10 * width)
    if (opt$objective > best$fitness) best <- list(fitness = opt$objective, at = opt$maximum)
  }
  best$at <- as.numeric(tf$inv(best$at))
  best
}

## The test functions that changed sign between two points.
boundary_crossed <- function(before, after) {
  nm <- intersect(names(before), names(after))
  nm[is.finite(before[nm]) & is.finite(after[nm]) & sign(before[nm]) != sign(after[nm]) &
       sign(before[nm]) != 0]
}

## Locate where the earliest of the `crossed` test functions changes sign on
## the step from u along t, whose far end `far` (corrected, with its tests)
## has them crossed: false position (Illinois) in the step length on the
## smallest of them, each signed positive before the step, until the change
## is bracketed to event_tol. A trial point is corrected without its Jacobian
## unless `bordered` is among them, which needs it; the point returned, at the
## crossed end of the bracket, has its Jacobian, tangent and every test.
## Where a corrector fails on the way the best point so far is returned.
boundary_locate_event <- function(residual, jacobian, u, t, J, tests, crossed, ds, far, info, control) {
  need_J <- "bordered" %in% crossed
  signs <- sign(tests[crossed])
  h <- function(pt) min(signs * pt$tests[crossed])
  at <- function(step) {
    out <- boundary_corrector(residual, jacobian, u, t, J, step, control, tangent = need_J)
    if (is.null(out$u)) return(NULL)
    out$tests <- boundary_test_functions(residual, out$u, out$J, out$t, info, control)
    out
  }
  lo <- 0
  h_lo <- h(list(tests = tests))
  hi <- ds
  h_hi <- h(far)
  best <- far
  kept <- 0L
  half <- control$event_tol / 2
  for (i in seq_len(100L)) {
    if (hi - lo <= control$event_tol) break
    ## a step within half of event_tol of an end crosses to close the bracket
    mid <- min(max((lo * h_hi - hi * h_lo) / (h_hi - h_lo), lo + half), hi - half)
    out <- at(mid)
    if (is.null(out)) break
    if (h(out) > 0) {
      lo <- mid
      h_lo <- h(out)
      if (kept == 1L) h_hi <- h_hi / 2
      kept <- 1L
    } else {
      hi <- mid
      h_hi <- h(out)
      best <- out
      if (kept == -1L) h_lo <- h_lo / 2
      kept <- -1L
    }
  }
  if (is.null(best$J)) {
    best$J <- jacobian(best$u)
    best$t <- boundary_tangent(best$J, t)
    best$tests <- boundary_test_functions(residual, best$u, best$J, best$t, info, control)
  }
  best
}

## Why u is off the branch's domain -- "left the box" (the free parameters
## outside `box`) or "left the trait bounds" (a resident or y* outside the
## community's) -- or NULL if it is inside.
boundary_outside <- function(u, info, box, start) {
  v <- boundary_unknowns(u, info)
  if (!is.null(box) && any(v$p[info$free] < box[, 1] | v$p[info$free] > box[, 2])) {
    return("left the box")
  }
  x <- rbind(matrix(v$x, info$m, info$k), if (!is.null(v$y)) matrix(v$y, 1L, info$k))
  b <- start$bounds
  if (!all(sweep(x, 2, b[, 1], ">=") & sweep(x, 2, b[, 2], "<="))) {
    return("left the trait bounds")
  }
  NULL
}

## The traces of both directions joined into one branch, in order of
## increasing arclength s, which is negative along the backward trace. Each
## trace's tangents point the way it travelled, so the backward ones are
## turned round; its first point is the start, which the forward trace also
## holds.
boundary_branch <- function(traces, info, seed, start, box, control) {
  forward <- traces$forward$points
  backward <- traces$backward$points
  flip <- function(pts) lapply(pts, function(pt) { pt$t <- -pt$t; pt$s <- -pt$s; pt })
  pts <- if (is.null(backward)) {
    forward
  } else if (is.null(forward)) {
    rev(flip(backward))
  } else {
    c(rev(flip(backward))[-length(backward)], forward)
  }
  start_row <- if (is.null(backward)) 1L else length(backward)
  free_names <- names(info$p)[info$free]
  labels <- resident_labels(info$trait_names, info$m)
  u <- do.call(rbind, lapply(pts, `[[`, "u"))
  tangent <- do.call(rbind, lapply(pts, `[[`, "t"))
  colnames(u) <- colnames(tangent) <- info$unknown_names
  v <- lapply(pts, function(pt) boundary_unknowns(pt$u, info))
  points <- data.frame(s = vapply(pts, `[[`, 0, "s"), check.names = FALSE)
  for (j in seq_along(free_names)) points[[free_names[j]]] <- vapply(v, function(x) x$p[[info$free[j]]], 0)
  for (j in seq_along(labels)) points[[labels[j]]] <- vapply(v, function(x) x$x[j], 0)
  if (identical(info$type, "invasion")) {
    for (j in seq_len(info$k)) points[[paste0("y:", info$trait_names[j])]] <- vapply(v, function(x) x$y[j], 0)
  }
  for (j in seq_len(info$m)) points[[sprintf("birth_rate[%d]", j)]] <- vapply(pts, function(pt) pt$birth_rate[j], 0)
  test_names <- pts[[1L]]$test_names
  for (nm in test_names) points[[nm]] <- vapply(pts, function(pt) pt$tests[match(nm, pt$test_names)], 0)
  points$iterations <- vapply(pts, `[[`, 0L, "iterations")

  events <- NULL
  for (end in names(traces)) {
    ev <- traces[[end]]$event
    if (is.null(ev)) next
    row <- if (end == "forward") nrow(points) else 1L
    rec <- data.frame(end = end, type = ev$type, s = points$s[row], stringsAsFactors = FALSE)
    for (nm in free_names) rec[[nm]] <- points[[nm]][row]
    rec$detail <- I(list(ev$detail))
    events <- rbind(events, rec)
  }
  if (is.null(events)) {
    events <- data.frame(end = character(0), type = character(0), s = numeric(0), stringsAsFactors = FALSE)
    for (nm in free_names) events[[nm]] <- numeric(0)
    events$detail <- I(list())
  }
  structure(list(
    points = points, u = u, tangent = tangent, start = start_row, events = events,
    termination = vapply(traces, `[[`, "", "termination"),
    rejected = sum(vapply(traces, `[[`, 0L, "rejected")),
    evaluations = sum(vapply(traces, `[[`, 0L, "evaluations")),
    type = info$type, free = free_names, which = info$which, p = info$p,
    trait_names = info$trait_names, m = info$m, k = info$k, box = box, info = info),
    class = "regnans_branch")
}

## Cold verification of every verify_every-th point and the ends, through
## the parallel map. An event's point is where two boundaries meet, so the
## labels a small offset either side of it cannot tell them apart: the point
## before it is verified instead.
boundary_verify <- function(branch, start, parameter, control) {
  n <- nrow(branch$u)
  ends <- c(backward = 1L, forward = n)
  for (end in branch$events$end) {
    ends[[end]] <- ends[[end]] + if (end == "forward") -1L else 1L
  }
  rows <- seq(1L, n, by = control$verify_every)
  rows <- sort(unique(c(rows[rows >= ends[["backward"]] & rows <= ends[["forward"]]], ends)))
  jobs <- lapply(rows, function(i) list(u = branch$u[i, ], t = branch$tangent[i, ]))
  out <- regnans_map(jobs, boundary_verify_point, start = start, parameter = parameter,
                     info = branch$info, control = control)
  failed <- vapply(out, inherits, logical(1), "error")
  if (any(failed)) {
    stop("boundary verification failed: ", conditionMessage(out[[which(failed)[1L]]]))
  }
  data.frame(
    s = branch$points$s[rows],
    residual = vapply(out, `[[`, 0, "residual"),
    n_minus = vapply(out, `[[`, 0L, "n_minus"), status_minus = vapply(out, `[[`, "", "status_minus"),
    n_plus = vapply(out, `[[`, 0L, "n_plus"), status_plus = vapply(out, `[[`, "", "status_plus"),
    ok = vapply(out, `[[`, TRUE, "ok"),
    evaluations = vapply(out, `[[`, 0L, "evaluations"),
    stringsAsFactors = FALSE)
}

## One verified point, on a worker: the residual re-solved from the model's
## initial densities (no warm start), and the labels either side of the curve
## across it in the plane of the free parameters.
boundary_verify_point <- function(job, start, parameter, info, control) {
  tryCatch({
    v <- boundary_unknowns(job$u, info)
    cold <- start
    cold$birth_rate <- rep(start$birth_rate_initial, info$m)
    residual <- boundary_residual_fn(cold, parameter, info$p, info$free, info$type, info$which)
    r <- residual(job$u)
    tp <- job$t[length(job$t) - 1:0]
    normal <- c(-tp[2], tp[1]) / sqrt(sum(tp^2))
    h <- control$verify_offset * pmax(1, abs(v$p[info$free]))
    side <- function(sgn) {
      q <- v$p
      q[info$free] <- q[info$free] + sgn * h * normal
      attr(community_label(parameter_community(parameter, q, info$trait_names), control), "label")
    }
    minus <- side(-1)
    plus <- side(1)
    ## one side is the boundary's own community, stable; across it at least
    ## one species more, or a count censored at the boundary's (DD99 beyond
    ## its branching locus fills the box)
    own <- function(l) l$n == info$m && l$status == "stable"
    more <- function(l) (l$status == "stable" && l$n == info$m + 1L) ||
      (l$status %in% c("capped", "edge") && l$n >= info$m)
    ok <- max(abs(r)) <= control$verify_tol && ((own(minus) && more(plus)) || (own(plus) && more(minus)))
    list(residual = max(abs(r)), n_minus = as.integer(minus$n), status_minus = minus$status,
         n_plus = as.integer(plus$n), status_plus = plus$status, ok = ok,
         evaluations = as.integer(attr(residual, "evaluations")() + minus$evaluations + plus$evaluations))
  }, error = function(e) e)
}

##' @param x A \code{regnans_branch}.
##' @param ... Ignored.
##' @rdname community_continue_boundary
##' @export
print.regnans_branch <- function(x, ...) {
  cat(sprintf("<regnans_branch: %s boundary of %d resident%s in %s, %d points, arclength %s>\n",
              x$type, x$m, if (x$m == 1L) "" else "s", paste(x$free, collapse = " x "),
              nrow(x$points), signif(diff(range(x$points$s)), 4)))
  for (end in names(x$termination)) {
    cat(sprintf("  %s: %s\n", end, x$termination[[end]]))
  }
  if (nrow(x$events) > 0L) {
    for (i in seq_len(nrow(x$events))) {
      cat(sprintf("  event (%s) %s at %s\n", x$events$end[i], x$events$type[i],
                  paste(x$free, signif(unlist(x$events[i, x$free]), 6), sep = " = ", collapse = ", ")))
    }
  }
  cat(sprintf("  verified: %d of %d points\n", sum(x$verification$ok), nrow(x$verification)))
  cat(sprintf("  equilibrium solves: %d (%d steps rejected)\n", x$evaluations, x$rejected))
  invisible(x)
}

##' @param row.names,optional Ignored.
##' @param n Points per step: with \code{n > 1} each step is interpolated by
##' the cubic Hermite polynomial through its ends and their tangents, so a
##' coarse branch plots as a smooth curve.
##' @rdname community_continue_boundary
##' @export
as.data.frame.regnans_branch <- function(x, row.names = NULL, optional = FALSE, n = 1L, ...) {
  n <- as.integer(n)
  if (n <= 1L) return(x$points)
  u <- x$u
  tg <- x$tangent
  s <- x$points$s
  rows <- list()
  for (i in seq_len(nrow(u) - 1L)) {
    h <- s[i + 1L] - s[i]
    tau <- seq(0, 1, length.out = n + 1L)[-(n + 1L)]
    h00 <- 2 * tau^3 - 3 * tau^2 + 1
    h10 <- tau^3 - 2 * tau^2 + tau
    h01 <- -2 * tau^3 + 3 * tau^2
    h11 <- tau^3 - tau^2
    seg <- outer(h00, u[i, ]) + h * outer(h10, tg[i, ]) + outer(h01, u[i + 1L, ]) + h * outer(h11, tg[i + 1L, ])
    rows[[i]] <- cbind(s = s[i] + tau * h, seg)
  }
  rows[[length(rows) + 1L]] <- cbind(s = s[nrow(u)], matrix(u[nrow(u), ], 1L, dimnames = list(NULL, colnames(u))))
  out <- as.data.frame(do.call(rbind, rows), check.names = FALSE)
  out
}

##' @param y Ignored.
##' @rdname community_continue_boundary
##' @export
plot.regnans_branch <- function(x, y, ...) {
  plot(boundary_diagram(x), ...)
}

##' The branch an event starts.
##'
##' Each event that ends a branch of \code{\link{community_continue_boundary}}
##' is a point where it meets another boundary, and this gives the start of
##' that boundary, to pass on to \code{community_continue_boundary()}:
##' \describe{
##'   \item{\code{separation}}{residents (or a resident and \eqn{y^*}) merge:
##'     each cluster (single linkage at 2.5 \code{merge_tol}, so that a pair
##'     straddling \eqn{y^*} joins it) becomes one resident at the
##'     density-weighted mean of its residents, and the boundary that
##'     continues is that resident's branching locus (the pitchfork the
##'     merging residents emerged from).}
##'   \item{\code{resident_curvature}}{a resident stops being an ESS: the
##'     boundary that starts is its branching locus, the same residents.}
##'   \item{\code{density}}{a resident goes extinct: the boundary that starts
##'     is the invasion boundary of the community without it, the resident's
##'     own trait the marginal mutant.}
##'   \item{\code{other_peak}}{another peak becomes marginal: the invasion
##'     boundary of the same residents at that peak.}
##' }
##' The annihilation of the marginal peak (\code{peak_curvature}) and a
##' crossing branch (\code{bordered}) start no boundary of these kinds, and
##' are an error.
##'
##' @title Start the branch an event seeds
##' @param branch A \code{regnans_branch} that ended at an event.
##' @param parameter The map the branch was traced with.
##' @param end Which end's event: \code{"forward"} or \code{"backward"}.
##' @return A list of arguments for \code{\link{community_continue_boundary}}:
##' \code{community}, \code{parameter}, \code{y}, \code{type}, \code{which},
##' \code{p} and \code{free}; pass it with \code{do.call()}.
##' @author Daniel Falster
##' @export
boundary_switch <- function(branch, parameter, end = c("forward", "backward")) {
  end <- match.arg(end)
  if (!inherits(branch, "regnans_branch")) stop("boundary_switch needs a regnans_branch")
  ev <- branch$events[branch$events$end == end, , drop = FALSE]
  if (nrow(ev) == 0L) {
    why <- if (end %in% names(branch$termination)) branch$termination[[end]] else "not traced"
    stop("the branch's ", end, " end is not an event (", why, ")")
  }
  row <- if (end == "forward") nrow(branch$u) else 1L
  info <- branch$info
  v <- boundary_unknowns(branch$u[row, ], info)
  m <- info$m
  k <- info$k
  x <- matrix(v$x, m, k)
  n <- as.numeric(unlist(branch$points[row, sprintf("birth_rate[%d]", seq_len(m))]))
  detail <- ev$detail[[1L]]
  type <- ev$type[1L]
  build <- function(traits, birth_rate) {
    community_clear_residents(parameter_community(parameter, v$p, info$trait_names)) |>
      community_add(trait_matrix(traits, info$trait_names), birth_rate = birth_rate)
  }
  args <- switch(type,
    separation = {
      group <- detail$separation[seq_len(m)]
      keep <- unique(group)
      traits <- t(vapply(keep, function(g) {
        i <- which(group == g)
        colSums(x[i, , drop = FALSE] * n[i]) / sum(n[i])
      }, numeric(k)))
      dim(traits) <- c(length(keep), k)
      merged <- match(detail$separation[length(detail$separation)], keep)
      if (identical(info$type, "branching") || is.na(merged)) {
        merged <- which(tabulate(group)[keep] > 1L)[1L]
      }
      list(community = build(traits, vapply(keep, function(g) sum(n[group == g]), 0)),
           y = NULL, type = "branching", which = merged)
    },
    resident_curvature = list(community = build(x, n), y = NULL, type = "branching",
                              which = detail$resident_curvature),
    density = {
      i <- detail$density
      list(community = build(x[-i, , drop = FALSE], n[-i]), y = x[i, ], type = "invasion", which = 1L)
    },
    other_peak = list(community = build(x, n), y = detail$other_peak, type = "invasion", which = 1L),
    stop("an event of type ", type, " starts no boundary that boundary_switch() can follow"))
  c(args, list(parameter = parameter, p = v$p, free = info$free))
}

##' Collect coexistence boundary branches into a diagram.
##'
##' @title A diagram of coexistence boundaries
##' @param ... \code{regnans_branch} objects, or lists of them.
##' @return An object of class \code{boundary_diagram}, a list of branches,
##' with \code{plot()} and \code{as.data.frame()} methods.
##' @author Daniel Falster
##' @export
boundary_diagram <- function(...) {
  branches <- list(...)
  branches <- unlist(lapply(branches, function(b) if (inherits(b, "regnans_branch")) list(b) else b),
                     recursive = FALSE)
  if (length(branches) == 0L || !all(vapply(branches, inherits, logical(1), "regnans_branch"))) {
    stop("boundary_diagram needs regnans_branch objects")
  }
  free <- branches[[1L]]$free
  if (!all(vapply(branches, function(b) identical(b$free, free), logical(1)))) {
    stop("the branches of a diagram must share their free parameters")
  }
  structure(branches, class = "boundary_diagram")
}

##' @param x A \code{boundary_diagram}.
##' @param row.names,optional Ignored.
##' @param n Points per step, as for a branch.
##' @param ... Ignored.
##' @rdname boundary_diagram
##' @export
as.data.frame.boundary_diagram <- function(x, row.names = NULL, optional = FALSE, n = 8L, ...) {
  do.call(rbind, lapply(seq_along(x), function(i) {
    d <- as.data.frame(x[[i]], n = n)
    data.frame(branch = i, type = x[[i]]$type, d[, c("s", x[[i]]$free)], check.names = FALSE)
  }))
}

##' @param y Ignored.
##' @rdname boundary_diagram
##' @export
plot.boundary_diagram <- function(x, y, n = 8L, ...) {
  free <- x[[1L]]$free
  curves <- as.data.frame(x, n = n)
  curves$branch <- factor(curves$branch)
  events <- do.call(rbind, lapply(x, function(b) {
    if (nrow(b$events) == 0L) return(NULL)
    data.frame(b$events[, free, drop = FALSE], event = b$events$type, check.names = FALSE)
  }))
  p <- ggplot2::ggplot(curves, ggplot2::aes(x = .data[[free[1]]], y = .data[[free[2]]],
                                            group = .data[["branch"]], linetype = .data[["type"]])) +
    ggplot2::geom_path() +
    ggplot2::labs(x = free[1], y = free[2], linetype = "boundary")
  if (!is.null(events)) {
    p <- p + ggplot2::geom_point(data = events, ggplot2::aes(x = .data[[free[1]]], y = .data[[free[2]]],
                                                             shape = .data[["event"]]),
                                 inherit.aes = FALSE)
  }
  p
}
