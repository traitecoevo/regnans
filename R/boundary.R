# Where in parameter space a model supports 0, 1, 2, ... coexisting species.
#
#   boundary_control()   the tunables of the label, the scouting and the rays
#   community_label()    the label of one community: the size of the
#                        evolutionarily stable community its parameters
#                        support, with a status saying whether the count is
#                        genuine or censored
#   boundary_design()    well-spread points in a parameter box
#   boundary_scout()     the label at each point of a design, through the
#                        parallel map
#   boundary_bisect()    where the label changes along rays from a point,
#                        each ray bisected on a worker
#   boundary_sweep()     the labels along a path, cold and carried from the
#                        previous point in each direction: where they
#                        disagree the model has alternative attractors
#   community_solve_boundary()
#                        a boundary point solved for: the residents at their
#                        singular coalition and a marginal condition, one
#                        parameter free, by Newton on the augmented system
#
# The label is the end state of assembly by immigration and evolution: the
# residents evolve to the singular coalition they form together, an immigrant
# arrives at the maximum of invasion fitness over the trait box if its fitness
# is positive there, and so on until nothing can invade. Every round is the
# assembler's (assembler_control(evolve = ), so there is one assembly loop),
# and this file only reads off what it ends on. Only a count the model
# decides is a label; a count the search decided (the box's edge reached, the
# cap on residents, a coalition that cannot be resolved) is censored and says
# so, since a boundary bisected on it would be the search's, not the model's.

##' Control how coexistence labels and boundaries are computed.
##'
##' Returns the settings read by \code{\link{community_label}} and
##' \code{\link{boundary_scout}}. Pass a list of values to override the
##' defaults.
##'
##' @title Coexistence boundary settings
##' @param control A list of values to modify from the defaults:
##' \describe{
##'   \item{\code{max_residents}}{the cap on the community's size; a label
##'     that would admit a resident beyond it is \code{"capped"}.}
##'   \item{\code{fitness_tol}}{an immigrant establishes when the maximum of
##'     invasion fitness over the box exceeds this: wherever fitness is
##'     positive, the default sitting just above the rounding in a resident's
##'     own fitness (zero at equilibrium); a model whose fitness is noisier
##'     (plant, about 1e-4) needs it raised.}
##'   \item{\code{edge_tol}}{an immigrant within this fraction of the bounds'
##'     width (on the trait scale) of the box's edge makes the label
##'     \code{"edge"}: the box, not the model, set the count.}
##'   \item{\code{branch_offset}}{how far from a resident at a branching point
##'     its daughter is placed, as a fraction of the bounds' width
##'     (\code{\link{assembler_control}}).}
##'   \item{\code{eps_too_close}}{a maximum within this fraction of the
##'     bounds' width of a resident is that resident at a branching point.}
##'   \item{\code{evolve}, \code{evolve_tol}}{how the residents evolve between
##'     immigrations (\code{\link{assembler_control}}): \code{"newton"}
##'     solves for the coalition and classifies it, integrating the canonical
##'     equation only when the root is not an attractor; \code{"canonical"}
##'     always integrates.}
##'   \item{\code{max_rounds}}{rounds of evolution and immigration before the
##'     label is given up as \code{"unresolved"}.}
##'   \item{\code{trait_tol}}{two communities of the same size whose sorted
##'     traits differ by less than this fraction of the bounds' width are the
##'     same community (a round that returns to one is a cycle, and a sweep
##'     whose routes end on different ones has found hysteresis).}
##'   \item{\code{bisect_tol}}{\code{\link{boundary_bisect}} halves each
##'     bracket until it is this short, as a fraction of the ray.}
##' }
##' For \code{\link{community_continue_boundary}}:
##' \describe{
##'   \item{\code{ds}, \code{ds_min}, \code{ds_max}}{the first, smallest
##'     and largest pseudo-arclength step, in the units of the unknowns
##'     (traits and the two free parameters together).}
##'   \item{\code{max_steps}}{steps in each direction before stopping.}
##'   \item{\code{tangent_tol}}{a step whose tangent turns so that its dot
##'     product with the last falls below this is rejected and the step
##'     halved: the corrector has jumped to another branch.}
##'   \item{\code{corrector_tol}, \code{corrector_maxit}}{the corrector's
##'     tolerance on every residual and its iteration limit.}
##'   \item{\code{density_tol}, \code{merge_tol}}{the test functions'
##'     thresholds: a resident whose share of the density falls below
##'     \code{density_tol} is going extinct, and two residents (or a resident
##'     and the marginal mutant) closer than \code{merge_tol} of the bounds'
##'     width have merged.}
##'   \item{\code{n_landscape}}{grid points on which the landscape is
##'     searched for another peak (one trait).}
##'   \item{\code{event_tol}}{an event is located by searching the step
##'     until it is bracketed this finely in arclength.}
##'   \item{\code{verify_every}, \code{verify_offset}, \code{verify_tol}}{every
##'     \code{verify_every} points (and the last) the branch is verified cold:
##'     the residual re-solved from the model's initial densities must be within
##'     \code{verify_tol}, and the labels a relative \code{verify_offset}
##'     either side of the curve, across it in the plane of the free
##'     parameters, must be \code{n} and \code{n + 1} species, both stable.}
##' }
##' @return A list with those elements.
##' @author Daniel Falster
##' @export
boundary_control <- function(control = NULL) {
  defaults <- list(
    max_residents = 8L,
    fitness_tol   = 1e-8,
    edge_tol      = 1e-3,
    branch_offset = 1e-2,
    eps_too_close = 1e-3,
    evolve        = "newton",
    evolve_tol    = 1e-6,
    max_rounds    = 20L,
    trait_tol     = 1e-3,
    bisect_tol    = 1e-6,
    ## continuation
    ds = 0.05, ds_min = 1e-5, ds_max = 0.2,
    max_steps = 200L,
    tangent_tol = 0.9,
    corrector_tol = 1e-8, corrector_maxit = 10L,
    density_tol = 1e-6, merge_tol = 1e-2,
    n_landscape = 101L,
    event_tol = 1e-6,
    verify_every = 15L, verify_offset = 1e-3, verify_tol = 1e-6
  )
  control <- as.list(control)
  extra <- setdiff(names(control), names(defaults))
  if (length(extra) > 0L) {
    stop("Unknown control parameters ", paste(extra, collapse = ", "))
  }
  ret <- modifyList(defaults, control)
  ret$evolve <- match.arg(ret$evolve, c("newton", "canonical"))
  for (nm in c("max_residents", "max_rounds", "max_steps", "corrector_maxit", "n_landscape", "verify_every")) {
    v <- ret[[nm]]
    if (!is.numeric(v) || length(v) != 1L || !is.finite(v) || v < 1 || v != round(v)) {
      stop(nm, " must be a positive whole number")
    }
    ret[[nm]] <- as.integer(v)
  }
  for (nm in c("edge_tol", "branch_offset", "eps_too_close", "evolve_tol", "trait_tol", "bisect_tol",
               "ds", "ds_min", "ds_max", "corrector_tol", "density_tol", "merge_tol", "event_tol",
               "verify_offset", "verify_tol")) {
    v <- ret[[nm]]
    if (!is.numeric(v) || length(v) != 1L || !is.finite(v) || v <= 0) {
      stop(nm, " must be a single positive number")
    }
  }
  if (!is.numeric(ret$fitness_tol) || length(ret$fitness_tol) != 1L || !is.finite(ret$fitness_tol)) {
    stop("fitness_tol must be a single finite number")
  }
  if (!(ret$ds_min <= ret$ds && ret$ds <= ret$ds_max)) {
    stop("the steps must satisfy ds_min <= ds <= ds_max")
  }
  if (!is.numeric(ret$tangent_tol) || length(ret$tangent_tol) != 1L || !(ret$tangent_tol > 0 && ret$tangent_tol < 1)) {
    stop("tangent_tol must be a single number between 0 and 1")
  }
  ret
}

label_statuses <- c("stable", "capped", "edge", "unresolved")

##' The number of species a community's parameters support.
##'
##' Assembles the community by immigration and evolution until it is
##' evolutionarily stable, and reports its size. Each round the residents
##' evolve to the singular coalition they form together and one immigrant
##' arrives at the maximum of invasion fitness over the trait box, if its
##' fitness there exceeds \code{fitness_tol}; the rounds are the assembler's
##' (\code{\link{assembler_run}} with \code{\link{assembler_control}}'s
##' \code{evolve}), started from \code{community} as given: empty, or seeded
##' with a resident for a model whose empty community cannot be invaded (soft
##' selection) or whose fundamental fitness has no interior maximum (DD99's is
##' flat). The box is the community's bounds, narrowed to the fundamental
##' niche where the harness can find it (plant, from an empty start).
##'
##' The status says whether the count is the model's:
##' \describe{
##'   \item{\code{"stable"}}{no immigrant can establish and every resident is
##'     at a singular coalition the dynamics reaches: an evolutionarily stable
##'     community, and the only genuine label.}
##'   \item{\code{"capped"}}{an immigrant would take the community beyond
##'     \code{max_residents}; it is refused. The count is at least \code{n}.}
##'   \item{\code{"edge"}}{the immigrant's maximum lies on the box's edge, so
##'     the box set the count; it is refused (a first immigrant into an empty
##'     community is exempt, since evolution moves it). At least \code{n}.}
##'   \item{\code{"unresolved"}}{the residents' evolution fails (the
##'     coalition cannot be solved, or the dynamics integrated from it do not
##'     settle), or the rounds cycle or run out: the count is the search's,
##'     with \code{reason} saying what happened.}
##' }
##' A censored label is never a boundary of the model; a bisection
##' (\code{\link{boundary_scout}}) treats it as "more than \code{n}" and an
##' unresolved one as unknown.
##'
##' @title Label a community by the species its parameters support
##' @param community A \code{community}, the start of the assembly.
##' @param control See \code{\link{boundary_control}}.
##' @return The community the assembly ends on (solved to demographic
##' equilibrium, at its singular coalition), with \code{attr(., "label")}: a
##' list of \code{n}, \code{status}, \code{reason}, \code{rounds} and
##' \code{evaluations}, the equilibrium solves it cost.
##' @examples
##' gk <- function(d) {
##'   community_start(bounds(x = c(-4, 4)), trait_scale = "linear",
##'                   harness = harness_gk98(d = d, sigma = 1)) |>
##'     community_add(trait_matrix(0.5, "x"), birth_rate = 1)
##' }
##' # three patches, farther apart: one generalist, then a pair of
##' # specialists, then one per patch
##' for (d in c(1, 1.5, 2.5)) print(attr(community_label(gk(d)), "label")[c("n", "status")])
##' @author Daniel Falster
##' @export
community_label <- function(community, control = boundary_control()) {
  control <- boundary_control(control)
  if (!inherits(community, "community")) {
    stop("community_label needs a community")
  }
  assembler <- assembler_control(list(
    birth_type = "maximum",
    evolve = control$evolve, evolve_tol = control$evolve_tol,
    branch_offset = control$branch_offset, eps_too_close = control$eps_too_close,
    eps_fitness_invasion = control$fitness_tol,
    max_residents = control$max_residents, edge_tol = control$edge_tol,
    dead_birth_rate = community$demography_control$equilibrium_extinct_birth_rate,
    ## the harness narrows the box to the fundamental niche of an empty
    ## community only
    compute_viable_fitness = length(community) == 0L))
  tf <- community_trait_transform(community)
  ## the box of a seeded start is its bounds; an empty one's is known once the
  ## assembler has narrowed it
  if (length(community) > 0L) label_check_bounds(community, tf)
  a <- assembler_start(community, assembler)
  rounds <- 0L
  outcome <- NULL
  if (isTRUE(a$done)) {
    outcome <- if (isFALSE(a$evolved$converged)) {
      list(status = "unresolved", reason = a$evolved$reason)
    } else {
      list(status = "stable", reason = NULL)
    }
  } else {
    label_check_bounds(a$community, tf)
    keys <- list(label_key(a$community, tf))
  }
  while (is.null(outcome)) {
    if (rounds >= control$max_rounds) {
      outcome <- list(status = "unresolved",
                      reason = sprintf("still assembling after %d rounds", rounds))
      break
    }
    a <- assembler_run(a, 1L)
    rounds <- rounds + 1L
    if (!is.null(a$refused)) {
      outcome <- a$refused
    } else if (isFALSE(a$evolved$converged)) {
      outcome <- list(status = "unresolved", reason = a$evolved$reason)
    } else if (isTRUE(a$done)) {
      outcome <- list(status = "stable", reason = NULL)
    } else {
      key <- label_key(a$community, tf)
      if (any(vapply(keys, label_same, logical(1), key, control$trait_tol))) {
        outcome <- list(status = "unresolved", reason = "the rounds returned to an earlier community")
      }
      keys <- c(keys, list(key))
    }
  }
  out <- a$community
  attr(out, "label") <- list(n = nrow(out$traits), status = outcome$status,
                             reason = outcome$reason, rounds = rounds,
                             evaluations = as.integer(a$evaluations))
  out
}

## The search, the edge rule and the keys all measure in fractions of the
## box's width on the trait scale, which must therefore be finite.
label_check_bounds <- function(community, tf) {
  b <- community$bounds
  if (!all(is.finite(tf$fwd(b)))) {
    stop("community_label needs bounds that are finite on the trait scale",
         if (tf$scale == "log") " (a positive lower bound for a log-scale trait)", "; these are ",
         paste(rownames(b), apply(signif(b, 6), 1, paste, collapse = " to "), sep = ": ", collapse = ", "),
         call. = FALSE)
  }
}

## A community's size and sorted traits on the trait scale, as a fraction of
## the bounds' width, for telling one round's community from another's.
label_key <- function(community, tf) {
  lb <- tf$fwd(community$bounds)
  z <- tf$fwd(community$traits)
  z <- sweep(sweep(matrix(z, ncol = ncol(community$traits)), 2, lb[, 1], "-"), 2, lb[, 2] - lb[, 1], "/")
  list(n = nrow(z), z = z[order(z[, 1]), , drop = FALSE])
}
label_same <- function(a, b, tol) {
  a$n == b$n && (a$n == 0L || max(abs(a$z - b$z)) < tol)
}

##' Well-spread points in a parameter box.
##'
##' A deterministic low-discrepancy design (the Kronecker sequence with the
##' generalised golden ratio) over a box of parameter values, for
##' \code{\link{boundary_scout}}: no random numbers, so the same design every
##' time, and each prefix of it is itself well spread.
##'
##' @title Scouting design in a parameter box
##' @param box A matrix with one row per parameter and columns \code{lower},
##' \code{upper}, as \code{\link{bounds}} builds.
##' @param n Number of points.
##' @return An \code{n} by \code{nrow(box)} matrix, columns named by parameter.
##' @examples
##' boundary_design(bounds(d = c(1, 3), kappa = c(0.5, 3)), 8)
##' @author Daniel Falster
##' @export
boundary_design <- function(box, n) {
  box <- check_bounds(box, finite = TRUE)
  if (!is.numeric(n) || length(n) != 1L || n < 1 || n != round(n)) {
    stop("n must be a positive whole number")
  }
  frac <- util_kronecker_points(as.integer(n), nrow(box))
  out <- sweep(sweep(frac, 2, box[, 2] - box[, 1], "*"), 2, box[, 1], "+")
  colnames(out) <- rownames(box)
  out
}

##' Label the points of a parameter design.
##'
##' The label (\code{\link{community_label}}) at each row of \code{design},
##' the model at each built by \code{parameter}: the scouting pass that
##' shows how many regions a parameter plane holds and roughly where. Every
##' point is labelled from the community \code{parameter} builds, apart from
##' the others, so under a \code{future::plan()} with several workers the
##' points run in parallel with the same answer (\code{parameter} runs on the
##' workers, and under a plan whose workers are separate processes it must
##' not hold a solved plant community).
##'
##' @title Scout a parameter design for coexistence
##' @param parameter A function \code{p -> community}; see
##' \code{\link{community_parameter_map}}. The community it returns is the
##' start of each assembly.
##' @param design A matrix of parameter values, one row per point, as
##' \code{\link{boundary_design}} builds; or a vector for one parameter.
##' @param control See \code{\link{boundary_control}}.
##' @return A data frame of class \code{boundary_labels}: the parameters, then
##' \code{n}, \code{status}, \code{reason}, \code{rounds}, \code{evaluations}
##' (equilibrium solves) and the list columns \code{traits} and
##' \code{birth_rate} of the community each point ended on. Total solves in
##' \code{attr(., "evaluations")}.
##' @examples
##' comm <- community_start(bounds(x = c(-2, 2)), trait_scale = "linear",
##'                         harness = harness_dd99()) |>
##'   community_add(trait_matrix(0.5, "x"), birth_rate = 100)
##' by <- community_parameter_map(comm, c("sigma_C", "sigma_K"))
##' # one species where competition is wider than the resource, else the box
##' # fills: no finite community resists invasion
##' boundary_scout(by, boundary_design(bounds(sigma_C = c(0.5, 1.5), sigma_K = c(0.8, 1.2)), 4))
##' @author Daniel Falster
##' @export
boundary_scout <- function(parameter, design, control = boundary_control()) {
  control <- boundary_control(control)
  if (!is.function(parameter)) {
    stop("parameter must be a function p -> community; ",
         "community_parameter_map() builds one for an explicit harness")
  }
  design <- boundary_check_design(design, parameter)
  points <- lapply(seq_len(nrow(design)), function(i) design[i, ])
  rows <- regnans_map(points, boundary_label_point, parameter = parameter, control = control)
  failed <- vapply(rows, inherits, logical(1), "error")
  if (any(failed)) {
    stop("community_label failed at ", paste(names(points[[1]]), signif(points[[which(failed)[1]]], 6),
                                             sep = " = ", collapse = ", "),
         ": ", conditionMessage(rows[[which(failed)[1]]]))
  }
  boundary_labels(design, rows)
}

## The design as a matrix with one column per parameter, named, its columns
## in the order the map takes them: a map applies p by position, so a design
## whose names are another order would be built with its parameters swapped
## and reported under their names. A map from community_parameter_map() knows
## its parameters' names; for a caller's function, `names` (a design already
## checked) sets the order, and otherwise the design's own names stand.
boundary_check_design <- function(design, parameter, names = NULL) {
  if (!is.matrix(design)) {
    design <- matrix(as.numeric(design), ncol = 1L)
  }
  if (!is.numeric(design) || nrow(design) < 1L || !all(is.finite(design))) {
    stop("design must be a matrix of finite parameter values, one row per point")
  }
  if (is.null(names) && !is.null(attr(parameter, "pars"))) {
    names <- names(harness_parameter_values(parameter(design[1, ])$harness$pars, attr(parameter, "pars")))
  }
  if (!is.null(names) && length(names) != ncol(design)) {
    stop("design must have one column per parameter, ", paste(names, collapse = ", "))
  }
  if (is.null(colnames(design))) {
    colnames(design) <- if (!is.null(names)) names
      else if (ncol(design) == 1L) "p" else sprintf("p[%d]", seq_len(ncol(design)))
  } else if (!is.null(names) && !identical(colnames(design), names)) {
    if (!setequal(colnames(design), names)) {
      stop("design's columns (", paste(colnames(design), collapse = ", "),
           ") must be the map's parameters, ", paste(names, collapse = ", "))
    }
    design <- design[, names, drop = FALSE]
  }
  design
}

## One point of a design, on a worker: the label of the community parameter(p)
## builds, returned as numbers (a solved plant community cannot cross back),
## or the error, so the caller can name the point.
boundary_label_point <- function(p, parameter, control) {
  tryCatch({
    comm <- parameter(p)
    if (!inherits(comm, "community")) {
      stop("parameter(p) must return a community")
    }
    out <- community_label(comm, control)
    c(attr(out, "label"), list(traits = unname(out$traits), birth_rate = as.numeric(out$birth_rate)))
  }, error = function(e) e)
}

## The data frame boundary_scout() returns.
boundary_labels <- function(design, rows) {
  out <- as.data.frame(design)
  out$n <- vapply(rows, function(r) as.integer(r$n), integer(1))
  out$status <- vapply(rows, function(r) r$status, character(1))
  out$reason <- vapply(rows, function(r) if (is.null(r$reason)) NA_character_ else r$reason, character(1))
  out$rounds <- vapply(rows, function(r) as.integer(r$rounds), integer(1))
  out$evaluations <- vapply(rows, function(r) as.integer(r$evaluations), integer(1))
  out$traits <- I(lapply(rows, `[[`, "traits"))
  out$birth_rate <- I(lapply(rows, `[[`, "birth_rate"))
  attr(out, "evaluations") <- sum(out$evaluations)
  class(out) <- c("boundary_labels", "data.frame")
  out
}

##' @param x A \code{boundary_labels} object.
##' @param ... Ignored.
##' @rdname boundary_scout
##' @export
print.boundary_labels <- function(x, ...) {
  pars <- setdiff(names(x), c("n", "status", "reason", "rounds", "evaluations", "traits", "birth_rate"))
  cat(sprintf("<boundary_labels: %d points over %s>\n", nrow(x), paste(pars, collapse = ", ")))
  tab <- table(label = sprintf("%d %s", x$n, x$status))
  for (nm in names(tab)) cat(sprintf("  %s: %d point%s\n", nm, tab[[nm]], if (tab[[nm]] == 1L) "" else "s"))
  cat(sprintf("  equilibrium solves: %d\n", attr(x, "evaluations")))
  invisible(x)
}

## The label as a string, for telling two labels apart: the count and its
## status. NA for an unresolved label, which is not known to differ from
## anything.
label_string <- function(row) {
  if (identical(row$status, "unresolved")) NA_character_ else sprintf("%d %s", as.integer(row$n), row$status)
}

##' Find where the coexistence label changes along rays.
##'
##' Each ray runs from \code{from} to one row of \code{to} in parameter
##' space, \eqn{p(t) = from + t (to - from)} for \eqn{t} in \eqn{[0, 1]}, and
##' is labelled at both ends (\code{\link{community_label}}). Where the labels
##' differ the ray is bisected: the label at the midpoint decides which half
##' holds a change, and a midpoint whose label differs from both ends puts a
##' change in each half, until a bracket is shorter than \code{bisect_tol}. A
##' label that only counts upwards along a ray, as nested regions make it
##' from an interior point, has one change per pair of labels, and bisection
##' cannot then pick the wrong one. A censored label is a label like any other
##' here (its count and status), so a change into or out of one is reported
##' with both labels and is a boundary of the search rather than of the model
##' unless both are \code{"stable"}; an unresolved label is not known to
##' differ from anything, so a bracket with one at an end, or meeting one at
##' its midpoint, is reported unrefined, with \code{resolved = FALSE}.
##'
##' Every label is assembled cold from the community \code{parameter} builds,
##' so the points are independent: each ray is one job of the parallel map,
##' with the same answer under any \code{future::plan()}. The rays' residents
##' either side of each change are what a solve for the boundary point starts
##' from.
##'
##' @title Bisect coexistence boundaries along rays
##' @param parameter A function \code{p -> community}; see
##' \code{\link{community_parameter_map}}.
##' @param from The parameter values every ray starts from: typically inside
##' the region of fewest species.
##' @param to The far end of each ray, one row per ray (a vector is one ray).
##' @param control See \code{\link{boundary_control}}.
##' @return A data frame of class \code{boundary_points}, one row per change
##' found: \code{ray}, the parameters at the bracket's midpoint, \code{t_lo}
##' and \code{t_hi} (the bracket along the ray), \code{n_lo}, \code{status_lo},
##' \code{n_hi}, \code{status_hi}, \code{resolved}, \code{reason}, and the list
##' columns \code{traits_lo}, \code{birth_rate_lo}, \code{traits_hi},
##' \code{birth_rate_hi}. \code{attr(., "rays")} has one row per ray with its
##' labels at the ends, its changes and its equilibrium solves, whose total is
##' \code{attr(., "evaluations")}.
##' @examples
##' comm <- community_start(bounds(x = c(-2, 2)), trait_scale = "linear",
##'                         harness = harness_dd99()) |>
##'   community_add(trait_matrix(0.5, "x"), birth_rate = 100)
##' by <- community_parameter_map(comm, c("sigma_C", "sigma_K"))
##' # from wide competition, where one species is stable, towards narrow
##' boundary_bisect(by, c(1.5, 1), rbind(c(0.5, 1.2), c(1.5, 2)),
##'                 boundary_control(list(bisect_tol = 1e-4)))
##' @author Daniel Falster
##' @export
boundary_bisect <- function(parameter, from, to, control = boundary_control()) {
  control <- boundary_control(control)
  if (!is.function(parameter)) {
    stop("parameter must be a function p -> community; ",
         "community_parameter_map() builds one for an explicit harness")
  }
  from_names <- if (is.matrix(from)) colnames(from) else names(from)
  from <- boundary_check_design(matrix(from, nrow = 1L, dimnames = list(NULL, from_names)), parameter)
  if (!is.matrix(to)) {
    to <- matrix(as.numeric(to), nrow = 1L, dimnames = list(NULL, names(to)))
  }
  if (ncol(to) != ncol(from)) {
    stop("to must have one column per parameter (", ncol(from), ")")
  }
  ## to's columns follow from's names wherever there are any to follow
  named <- !is.null(from_names) || !is.null(attr(parameter, "pars"))
  to <- boundary_check_design(to, parameter, names = if (named) colnames(from))
  colnames(to) <- colnames(from)
  if (any(rowSums(abs(sweep(to, 2, from[1, ]))) == 0)) {
    stop("a ray must end away from where it starts")
  }
  ends <- lapply(seq_len(nrow(to)), function(i) to[i, ])
  rays <- regnans_map(ends, boundary_ray, from = from[1, ], parameter = parameter, control = control)
  failed <- vapply(rays, inherits, logical(1), "error")
  if (any(failed)) {
    stop("boundary_bisect failed on ray ", which(failed)[1], ": ", conditionMessage(rays[[which(failed)[1]]]))
  }
  boundary_points(from[1, ], to, rays)
}

## One ray, on a worker: its labels bisected wherever they change. Returns
## numbers only, or the error.
boundary_ray <- function(to, from, parameter, control) {
  tryCatch({
    evaluations <- 0L
    label <- function(t) {
      r <- boundary_label_point(from + t * (to - from), parameter, control)
      if (inherits(r, "error")) stop(r)
      evaluations <<- evaluations + as.integer(r$evaluations)
      r
    }
    changes <- list()
    record <- function(t_lo, lo, t_hi, hi, resolved, reason = NA_character_) {
      changes[[length(changes) + 1L]] <<- list(t_lo = t_lo, t_hi = t_hi, lo = lo, hi = hi,
                                               resolved = resolved, reason = reason)
    }
    search <- function(t_lo, lo, t_hi, hi) {
      a <- label_string(lo)
      b <- label_string(hi)
      if (!is.na(a) && !is.na(b) && a == b) {
        return(invisible())
      }
      if (is.na(a) || is.na(b)) {
        return(record(t_lo, lo, t_hi, hi, FALSE, "an end of the bracket is unresolved"))
      }
      if (t_hi - t_lo <= control$bisect_tol) {
        return(record(t_lo, lo, t_hi, hi, TRUE))
      }
      t_mid <- (t_lo + t_hi) / 2
      mid <- label(t_mid)
      m <- label_string(mid)
      if (is.na(m)) {
        return(record(t_lo, lo, t_hi, hi, FALSE,
                      sprintf("the label is unresolved at t = %s inside the bracket: %s",
                              signif(t_mid, 6), mid$reason)))
      }
      search(t_lo, lo, t_mid, mid)
      search(t_mid, mid, t_hi, hi)
    }
    start <- label(0)
    end <- label(1)
    search(0, start, 1, end)
    list(start = start, end = end, changes = changes, evaluations = evaluations)
  }, error = function(e) e)
}

## The data frame boundary_bisect() returns.
boundary_points <- function(from, to, rays) {
  pars <- colnames(to)
  rows <- list()
  for (i in seq_along(rays)) {
    for (ch in rays[[i]]$changes) {
      t <- (ch$t_lo + ch$t_hi) / 2
      p <- from + t * (to[i, ] - from)
      rows[[length(rows) + 1L]] <- c(list(ray = i), as.list(stats::setNames(p, pars)),
        list(t_lo = ch$t_lo, t_hi = ch$t_hi,
             n_lo = as.integer(ch$lo$n), status_lo = ch$lo$status,
             n_hi = as.integer(ch$hi$n), status_hi = ch$hi$status,
             resolved = ch$resolved, reason = ch$reason,
             traits_lo = list(ch$lo$traits), birth_rate_lo = list(ch$lo$birth_rate),
             traits_hi = list(ch$hi$traits), birth_rate_hi = list(ch$hi$birth_rate)))
    }
  }
  scalar <- c("ray", pars, "t_lo", "t_hi", "n_lo", "status_lo", "n_hi", "status_hi", "resolved", "reason")
  listcols <- c("traits_lo", "birth_rate_lo", "traits_hi", "birth_rate_hi")
  out <- if (length(rows) == 0L) {
    df <- data.frame(ray = integer(0))
    for (nm in pars) df[[nm]] <- numeric(0)
    df$t_lo <- numeric(0); df$t_hi <- numeric(0)
    df$n_lo <- integer(0); df$status_lo <- character(0)
    df$n_hi <- integer(0); df$status_hi <- character(0)
    df$resolved <- logical(0); df$reason <- character(0)
    for (nm in listcols) df[[nm]] <- I(list())
    df
  } else {
    df <- as.data.frame(lapply(stats::setNames(scalar, scalar), function(nm) {
      v <- lapply(rows, `[[`, nm)
      unlist(v)
    }), stringsAsFactors = FALSE, check.names = FALSE)
    for (nm in listcols) df[[nm]] <- I(lapply(rows, function(r) r[[nm]][[1L]]))
    df
  }
  summary <- data.frame(
    ray = seq_along(rays),
    n_start = vapply(rays, function(r) as.integer(r$start$n), integer(1)),
    status_start = vapply(rays, function(r) r$start$status, character(1)),
    n_end = vapply(rays, function(r) as.integer(r$end$n), integer(1)),
    status_end = vapply(rays, function(r) r$end$status, character(1)),
    changes = vapply(rays, function(r) length(r$changes), integer(1)),
    evaluations = vapply(rays, function(r) as.integer(r$evaluations), integer(1)),
    stringsAsFactors = FALSE)
  attr(out, "rays") <- summary
  attr(out, "evaluations") <- sum(summary$evaluations)
  class(out) <- c("boundary_points", "data.frame")
  out
}

##' @param x A \code{boundary_points} object.
##' @param ... Ignored.
##' @rdname boundary_bisect
##' @export
print.boundary_points <- function(x, ...) {
  rays <- attr(x, "rays")
  cat(sprintf("<boundary_points: %d change%s on %d ray%s>\n", nrow(x), if (nrow(x) == 1L) "" else "s",
              nrow(rays), if (nrow(rays) == 1L) "" else "s"))
  pars <- setdiff(names(x), c("ray", "t_lo", "t_hi", "n_lo", "status_lo", "n_hi", "status_hi", "resolved",
                              "reason", "traits_lo", "birth_rate_lo", "traits_hi", "birth_rate_hi"))
  for (i in seq_len(nrow(x))) {
    cat(sprintf("  ray %d at %s: %d %s -> %d %s%s\n", x$ray[i],
                paste(pars, signif(unlist(x[i, pars]), 6), sep = " = ", collapse = ", "),
                x$n_lo[i], x$status_lo[i], x$n_hi[i], x$status_hi[i],
                if (x$resolved[i]) "" else " (unresolved)"))
  }
  cat(sprintf("  equilibrium solves: %d\n", attr(x, "evaluations")))
  invisible(x)
}

##' Look for alternative attractors along a path in parameter space.
##'
##' Where a model has alternative stable communities the label is not a
##' function of the parameters alone: which community is reached depends on
##' where the assembly started. A sweep labels each point of a path three
##' ways: cold, from the community \code{parameter} builds
##' (\code{\link{boundary_scout}}); and carried, from the community the
##' previous point ended on, along the path forward and back. Each carried
##' label starts from its predecessor's residents and densities, so it
##' follows an attractor for as long as it persists, as a slowly changing
##' environment would. A point whose three labels differ in count or status,
##' or whose communities differ in the residents' sorted traits by more than
##' \code{trait_tol} of the bounds' width, is flagged: the history decided. An
##' unresolved label is not known to differ from anything, so it is left out
##' of the comparison, and a point where the routes that did resolve agree is
##' flagged \code{NA}. Every route searches the same box: a carried label
##' starts from the box the cold one would search, narrowed to the
##' fundamental niche where the harness can find it (plant), and its
##' predecessor's residents outside that box are dropped.
##'
##' The cold labels are independent jobs of the parallel map; each carried
##' direction is one job (a chain, sequential by nature). The answer does not
##' depend on the plan.
##'
##' @title Sweep a path for hysteresis
##' @param parameter A function \code{p -> community}; see
##' \code{\link{community_parameter_map}}.
##' @param p The path, one row per point in order (a vector for one
##' parameter).
##' @param control See \code{\link{boundary_control}}.
##' @return A data frame of class \code{boundary_sweep}, one row per point:
##' the parameters, then \code{n} and \code{status} for each route
##' (\code{cold_n}, \code{cold_status}, \code{forward_n}, ...,
##' \code{backward_status}), \code{hysteresis} (\code{NA} where an
##' unresolved route leaves it unknown), the list columns
##' \code{cold_traits}, \code{forward_traits} and \code{backward_traits}, and
##' \code{evaluations} (equilibrium solves, all three routes), whose total is
##' \code{attr(., "evaluations")}.
##' @examples
##' comm <- community_start(bounds(x = c(-2, 2)), trait_scale = "linear",
##'                         harness = harness_dd99()) |>
##'   community_add(trait_matrix(0.5, "x"), birth_rate = 100)
##' # one stable community at every point: nothing to flag
##' boundary_sweep(community_parameter_map(comm, "sigma_C"), c(1.2, 1.4, 1.6))
##' @author Daniel Falster
##' @export
boundary_sweep <- function(parameter, p, control = boundary_control()) {
  control <- boundary_control(control)
  if (!is.function(parameter)) {
    stop("parameter must be a function p -> community; ",
         "community_parameter_map() builds one for an explicit harness")
  }
  p <- boundary_check_design(p, parameter)
  if (nrow(p) < 2L) {
    stop("a sweep needs a path of at least two points")
  }
  points <- lapply(seq_len(nrow(p)), function(i) p[i, ])
  cold <- regnans_map(points, boundary_label_point, parameter = parameter, control = control)
  chains <- regnans_map(list(points, rev(points)), boundary_chain, parameter = parameter, control = control)
  for (r in c(cold, chains)) {
    if (inherits(r, "error")) stop("boundary_sweep failed: ", conditionMessage(r))
  }
  forward <- chains[[1L]]
  backward <- rev(chains[[2L]])

  base <- parameter(p[1L, ])
  tf <- community_trait_transform(base)
  key <- function(r) {
    comm <- base
    comm$traits <- trait_matrix(r$traits, base$trait_names)
    label_key(comm, tf)
  }
  out <- as.data.frame(p)
  for (route in c("cold", "forward", "backward")) {
    rows <- get(route)
    out[[paste0(route, "_n")]] <- vapply(rows, function(r) as.integer(r$n), integer(1))
    out[[paste0(route, "_status")]] <- vapply(rows, function(r) r$status, character(1))
  }
  out$hysteresis <- vapply(seq_len(nrow(p)), function(i) {
    boundary_hysteresis(list(cold[[i]], forward[[i]], backward[[i]]), key, control$trait_tol)
  }, logical(1))
  out$cold_traits <- I(lapply(cold, `[[`, "traits"))
  out$forward_traits <- I(lapply(forward, `[[`, "traits"))
  out$backward_traits <- I(lapply(backward, `[[`, "traits"))
  out$evaluations <- vapply(seq_len(nrow(p)), function(i) {
    as.integer(cold[[i]]$evaluations + forward[[i]]$evaluations + backward[[i]]$evaluations)
  }, integer(1))
  attr(out, "evaluations") <- sum(out$evaluations)
  class(out) <- c("boundary_sweep", "data.frame")
  out
}

## Whether the routes' labels at one point of a sweep disagree. An unresolved
## route is not known to differ from anything: two resolved routes that
## disagree (in count and status, or in their communities' keys) are
## hysteresis whatever the others did, and resolved routes that agree beside
## an unresolved one leave it unknown, NA.
boundary_hysteresis <- function(routes, key, trait_tol) {
  labels <- vapply(routes, label_string, character(1))
  resolved <- routes[!is.na(labels)]
  if (length(unique(labels[!is.na(labels)])) > 1L) return(TRUE)
  if (length(resolved) >= 2L) {
    keys <- lapply(resolved, key)
    if (!all(vapply(keys[-1L], label_same, logical(1), keys[[1L]], trait_tol))) return(TRUE)
  }
  if (length(resolved) < length(routes)) NA else FALSE
}

## One direction of a sweep, on a worker: each point labelled from the
## community the previous one ended on (the first from the map's own start).
## The box is the one the cold label at that point searches: a label narrows
## the box to the fundamental niche only from an empty start, so a map that
## starts empty has its box narrowed here before the residents are carried in,
## and a carried resident outside it is dropped. Returns a list of label
## results in path order, or the error.
boundary_chain <- function(points, parameter, control) {
  tryCatch({
    prev <- NULL
    lapply(points, function(p) {
      comm <- parameter(p)
      if (!inherits(comm, "community")) {
        stop("parameter(p) must return a community")
      }
      if (!is.null(prev) && nrow(prev$traits) > 0L) {
        box <- if (length(comm) == 0L) community_viable_bounds(community_demography(comm))$bounds
               else comm$bounds
        ## no viable box: nothing to carry, and the label finds the same
        if (!is.null(box)) {
          comm$bounds <- box
          x <- trait_matrix(prev$traits, comm$trait_names)
          inside <- apply(sweep(x, 2, comm$bounds[, 1], ">=") & sweep(x, 2, comm$bounds[, 2], "<="), 1, all)
          comm <- community_clear_residents(comm) |>
            community_add(x[inside, , drop = FALSE], birth_rate = prev$birth_rate[inside])
        }
      }
      out <- community_label(comm, control)
      prev <<- list(traits = unname(out$traits), birth_rate = as.numeric(out$birth_rate))
      c(attr(out, "label"), prev)
    })
  }, error = function(e) e)
}

##' @param x A \code{boundary_sweep} object.
##' @param ... Ignored.
##' @rdname boundary_sweep
##' @export
print.boundary_sweep <- function(x, ...) {
  cat(sprintf("<boundary_sweep: %d points, %d with hysteresis%s>\n", nrow(x), sum(x$hysteresis, na.rm = TRUE),
              if (anyNA(x$hysteresis)) sprintf(", %d unknown", sum(is.na(x$hysteresis))) else ""))
  cat(sprintf("  equilibrium solves: %d\n", attr(x, "evaluations")))
  invisible(x)
}

##' Solve for a point on a coexistence boundary.
##'
##' A boundary between regions supporting \eqn{n} and \eqn{n + 1} species is
##' the solution set of an augmented system: the \eqn{n} residents at their
##' singular coalition, each at demographic equilibrium, and a condition that
##' makes the community marginal. Two kinds of boundary:
##' \describe{
##'   \item{\code{"invasion"}}{a new peak of the fitness landscape, away from
##'     the residents, crossing zero: a mutant \eqn{y^*} at a critical point of
##'     invasion fitness, \eqn{\partial s / \partial y (y^*) = 0}, with
##'     \eqn{s(y^*) = 0}. Unknowns: the residents' traits, \eqn{y^*} and one
##'     parameter; \eqn{nk + k + 1} equations in as many unknowns.}
##'   \item{\code{"branching"}}{a resident becoming a branching point, where
##'     the pair it branches into emerges from it (a pitchfork): the
##'     determinant of resident \code{which}'s Hessian in the mutant direction
##'     is zero. Unknowns: the residents' traits and one parameter.}
##' }
##' With the other parameters fixed the system is square and is solved by
##' the package's Newton iteration (\code{\link{util_nlsolve}}, \code{solver
##' = "newton"}), the densities solved to equilibrium
##' inside each evaluation (\code{\link{community_demography}}). Its Jacobian
##' comes from the derivative layer: the residents' columns by finite
##' differences across the residents, each point solved to the cube of its
##' step and the points through the parallel map, the marginal mutant's from
##' the exact mutant-direction derivatives, the free parameter's by finite
##' differences across the parameter; so \code{2nk + 2} equilibrium solves per
##' Jacobian until the model supplies equilibrium sensitivities.
##'
##' The start is a boundary point found by \code{\link{boundary_bisect}}: the
##' community at the side with fewer species, \code{y} the newcomer's trait
##' on the other side, and the bracket's midpoint as \code{p}. \eqn{y^*} must
##' stay away from the residents: at a resident the conditions hold for any
##' parameter (its own fitness is zero, and its gradient too at a singular
##' coalition), so a start or an answer there is refused. An answer at which
##' \eqn{y^*} is not a strict maximum of fitness warns: there the boundary is
##' not where a new peak crosses zero.
##'
##' Traits and parameters are solved in raw units, whatever the community's
##' trait scale.
##'
##' @title Solve for a coexistence boundary point
##' @param community A \code{community} whose residents start the solve (their
##' densities, if positive, seed the first equilibrium), built by
##' \code{parameter(p)}.
##' @param parameter A function \code{p -> community}; see
##' \code{\link{community_parameter_map}}.
##' @param y For \code{type = "invasion"}, the marginal mutant's starting
##' traits, one value per trait.
##' @param type \code{"invasion"} or \code{"branching"}.
##' @param free The parameter solved for, by name or position in \code{p};
##' the others are held.
##' @param which For \code{type = "branching"}, the resident becoming a
##' branching point.
##' @param p The parameter values \code{community} was built at, the free
##' one's starting value among them. Read from the community for a map from
##' \code{\link{community_parameter_map}}.
##' @param tol,maxit Passed to \code{\link{util_nlsolve}}: the tolerance on
##' every residual and the iteration limit.
##' @return The community at the boundary point, solved to demographic
##' equilibrium at the parameters found, with \code{attr(., "boundary")}: a list
##' of \code{type}, \code{p} (all the parameters, named), \code{free},
##' \code{y} (\code{NULL} for a branching boundary), \code{residual},
##' \code{converged}, \code{message}, \code{iterations} and
##' \code{evaluations} (equilibrium solves).
##' @examples
##' gk <- function(p) {
##'   community_start(bounds(x = c(-4, 4)), trait_scale = "linear",
##'                   harness = harness_gk98(mu = c(-p[1], 0, p[1]))) |>
##'     community_add(trait_matrix(0.5, "x"), birth_rate = 1)
##' }
##' # the pair of specialists on the outer patches, and the middle patch's
##' # peak, near the spacing at which a third species invades
##' pair <- community_start(bounds(x = c(-4, 4)), trait_scale = "linear",
##'                         harness = harness_gk98(mu = c(-1.95, 0, 1.95))) |>
##'   community_add(trait_matrix(c(-1.3, 1.3), "x"), birth_rate = c(1.5, 1.5))
##' at <- community_solve_boundary(pair, gk, y = 0, p = 1.95)
##' attr(at, "boundary")$p
##' @author Daniel Falster
##' @export
community_solve_boundary <- function(community, parameter, y = NULL,
                                     type = c("invasion", "branching"), free = 1L,
                                     which = 1L, p = NULL, tol = 1e-8, maxit = 50) {
  type <- match.arg(type)
  p <- check_parameter_values(parameter, p, community)
  trait_names <- community$trait_names
  k <- length(trait_names)
  m <- nrow(community$traits)
  if (m < 1L) {
    stop("community_solve_boundary needs the residents of the side with fewer species, at least one")
  }
  if (length(free) != 1L || is.na(free)) {
    stop("free names one parameter to solve for")
  }
  free <- if (is.character(free)) match(free, names(p)) else as.integer(free)
  if (is.na(free) || free < 1L || free > length(p)) {
    stop("free must be one of the parameters, ", paste(names(p), collapse = ", "))
  }
  if (type == "invasion") {
    if (is.null(y) || length(y) != k || !all(is.finite(y))) {
      stop("an invasion boundary needs y, the marginal mutant's starting traits, one value per trait (", k, ")")
    }
    y <- as.numeric(y)
    boundary_check_away(community, y, "y starts")
  } else {
    if (!is.null(y)) stop("a branching boundary takes no y")
    if (length(which) != 1L || !which %in% seq_len(m)) {
      stop("which must name one of the ", m, " residents")
    }
  }
  merged <- residents_merged(as.numeric(community$traits), m, k, community_derivative_control(community))
  if (nrow(merged) > 0L) {
    stop("residents ", merged[1, 1], " and ", merged[1, 2], " start within the finite-difference step of ",
         "each other: two copies of one strategy are not a coalition")
  }

  residual <- boundary_residual_fn(community, parameter, p, free, type, which)
  info <- attr(residual, "info")
  u0 <- c(as.numeric(community$traits), y, p[[free]])
  jacobian_solves <- 0L
  jac <- function(u) {
    J <- community_boundary_jacobian(residual, u)
    jacobian_solves <<- jacobian_solves + attr(J, "evaluations")
    unname(J)
  }
  ## the package's Newton, which recomputes the Jacobian when a Broyden-updated
  ## one fails its line search; nleqslv keeps updating it, and on the GK98
  ## oracle its update turns singular two steps from the root
  sol <- util_nlsolve(u0, residual, tol = tol, maxit = maxit, solver = "newton",
                      require_converged = FALSE, jac = jac)
  r <- residual(as.numeric(sol))
  out <- attr(residual, "last")()
  v <- boundary_unknowns(as.numeric(sol), info)
  converged <- isTRUE(attr(sol, "converged"))
  message <- attr(sol, "message")
  if (type == "invasion") {
    if (converged && !boundary_check_away(out, v$y, "y ended", stop = FALSE)) {
      converged <- FALSE
      message <- "the marginal mutant converged onto a resident, where the conditions hold trivially"
    }
    H <- community_fitness_hessian(out, matrix(v$y, nrow = 1L))
    if (converged && !all(eigen(H, symmetric = TRUE, only.values = TRUE)$values < 0)) {
      warning("community_solve_boundary: y* is not a strict maximum of invasion fitness, ",
              "so this is not where a new peak crosses zero")
    }
  }
  if (!converged) {
    warning(sprintf("community_solve_boundary did not converge (%s)", message))
  }
  attr(out, "boundary") <- list(
    type = type, p = v$p, free = names(p)[free], y = v$y,
    residual = stats::setNames(r, info$residual_names), converged = converged, message = message,
    iterations = attr(sol, "iter"),
    evaluations = attr(residual, "evaluations")() + jacobian_solves)
  out
}

## Is y more than two finite-difference steps from every resident? Errors,
## naming the fix, unless stop = FALSE.
boundary_check_away <- function(community, y, what, stop = TRUE) {
  m <- nrow(community$traits)
  k <- length(community$trait_names)
  x <- rbind(community$traits, y, deparse.level = 0)
  near <- residents_merged(as.numeric(x), m + 1L, k, community_derivative_control(community), reach = 2)
  near <- near[near[, 2] == m + 1L, , drop = FALSE]
  if (nrow(near) == 0L) return(TRUE)
  if (!stop) return(FALSE)
  stop(what, " at resident ", near[1, 1], ": the marginal peak must be away from the residents. ",
       "Give the residents of the side with fewer species, and as y the newcomer's traits from the other side")
}

## The unknowns of a boundary solve, split: residents' traits x, the
## marginal mutant y (invasion only) and the parameters with the free ones
## set (one for a point, two for a continuation), last in u.
boundary_unknowns <- function(u, info) {
  mk <- info$m * info$k
  y <- if (identical(info$type, "invasion")) u[mk + seq_len(info$k)] else NULL
  p <- info$p
  nf <- length(info$free)
  p[info$free] <- u[length(u) - nf + seq_len(nf)]
  list(x = u[seq_len(mk)], y = y, p = p)
}

## What boundary_probe() evaluates at a solved community, and how many values
## it returns.
boundary_probe_spec <- function(info, y) {
  list(type = info$type, y = y, which = info$which,
       n = if (identical(info$type, "invasion")) info$k + 1L else 1L)
}

## The residual of a boundary solve: unknowns u -> the residents' selection
## gradients and the marginal conditions (free: the parameters in u, one or
## two), the community solved to
## equilibrium at the residents under parameter(p), warm-started from the
## previous evaluation's densities and solver state. A call at the previous
## unknowns returns that answer without solving. Residents that have merged
## (residents_merged()) are outside the domain, as for
## singularity_gradient_fn(): the residual is NA and nothing is solved, so a
## Newton step or corrector backtracks.
boundary_residual_fn <- function(community, parameter, p, free, type, which) {
  trait_names <- community$trait_names
  k <- length(trait_names)
  m <- nrow(community$traits)
  labels <- resident_labels(trait_names, m)
  info <- list(m = m, k = k, type = type, which = which, p = p, free = free,
               parameter = parameter, trait_names = trait_names,
               residual_names = c(paste0("g:", labels),
                                  if (type == "invasion") c(paste0("ds/d", trait_names, "(y)"), "s(y)")
                                  else sprintf("det(H[%d])", which)),
               unknown_names = c(labels, if (type == "invasion") paste0("y:", trait_names), names(p)[free]))
  seed <- singularity_seed_birth_rate(community, m)
  state <- NULL
  last_u <- NULL
  last_r <- NULL
  last_community <- NULL
  evaluations <- 0L
  ctrl <- community_derivative_control(community)
  fn <- function(u) {
    u <- as.numeric(u)
    if (identical(u, last_u)) {
      return(last_r)
    }
    v <- boundary_unknowns(u, info)
    if (nrow(residents_merged(v$x, m, k, ctrl)) > 0L) {
      return(rep(NA_real_, length(info$residual_names)))
    }
    base <- community_clear_residents(parameter_community(parameter, v$p, trait_names))
    comm <- singularity_solve_point(v$x, base, m, seed, state)
    evaluations <<- evaluations + 1L
    r <- c(as.numeric(comm$selection_gradient), boundary_probe(comm, boundary_probe_spec(info, v$y)))
    br <- as.numeric(comm$birth_rate)
    if (length(br) == m && all(is.finite(br) & br > 0)) seed <<- br
    state <<- comm$demography_state
    last_u <<- u
    last_r <<- r
    last_community <<- comm
    r
  }
  attr(fn, "info") <- info
  attr(fn, "last") <- function() last_community
  attr(fn, "evaluations") <- function() evaluations
  fn
}
