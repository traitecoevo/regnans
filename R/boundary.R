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
##'     same community (a round that returns to one is a cycle).}
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
    trait_tol     = 1e-3
  )
  control <- as.list(control)
  extra <- setdiff(names(control), names(defaults))
  if (length(extra) > 0L) {
    stop("Unknown control parameters ", paste(extra, collapse = ", "))
  }
  ret <- modifyList(defaults, control)
  ret$evolve <- match.arg(ret$evolve, c("newton", "canonical"))
  for (nm in c("max_residents", "max_rounds")) {
    v <- ret[[nm]]
    if (!is.numeric(v) || length(v) != 1L || !is.finite(v) || v < 1 || v != round(v)) {
      stop(nm, " must be a positive whole number")
    }
    ret[[nm]] <- as.integer(v)
  }
  for (nm in c("edge_tol", "branch_offset", "eps_too_close", "evolve_tol", "trait_tol")) {
    v <- ret[[nm]]
    if (!is.numeric(v) || length(v) != 1L || !is.finite(v) || v <= 0) {
      stop(nm, " must be a single positive number")
    }
  }
  if (!is.numeric(ret$fitness_tol) || length(ret$fitness_tol) != 1L || !is.finite(ret$fitness_tol)) {
    stop("fitness_tol must be a single finite number")
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
  a <- assembler_start(community, assembler)
  tf <- community_trait_transform(community)
  keys <- list(label_key(a$community, tf))
  rounds <- 0L
  outcome <- NULL
  if (isTRUE(a$done)) {
    outcome <- if (isFALSE(a$evolved$converged)) {
      list(status = "unresolved", reason = a$evolved$reason)
    } else {
      list(status = "stable", reason = NULL)
    }
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

## The design as a matrix with one column per parameter, named.
boundary_check_design <- function(design, parameter) {
  if (!is.matrix(design)) {
    design <- matrix(as.numeric(design), ncol = 1L)
  }
  if (!is.numeric(design) || nrow(design) < 1L || !all(is.finite(design))) {
    stop("design must be a matrix of finite parameter values, one row per point")
  }
  if (is.null(colnames(design))) {
    colnames(design) <- if (!is.null(attr(parameter, "pars"))) {
      names(harness_parameter_values(parameter(design[1, ])$harness$pars, attr(parameter, "pars")))
    } else if (ncol(design) == 1L) "p" else sprintf("p[%d]", seq_len(ncol(design)))
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
