# How singular strategies and coalitions move with model parameters.
#
#   community_parameter_map()          p -> community for an explicit harness
#   community_parameter_sensitivity()  dx*/dp = -(dg/dx)^-1 dg/dp at a singularity
#   community_continue_singularity()   x* followed along one parameter: the
#                                      sensitivity predicts, the singularity
#                                      solver corrects, the classifier marks
#                                      where the verdict changes
#
# A parameter enters through a function p -> community rather than a harness
# connector: an explicit harness closes its pars into the model primitives,
# plant's parameters live in its own Parameters object, and a parameter worth
# following may be an environmental input that no harness owns. Derivatives
# come from R/derivatives.R: dg/dx is community_selection_gradient_jacobian()
# (any number of residents) and dg/dp
# community_selection_gradient_parameter_jacobian().

##' Vary the parameters of a model with an explicit harness.
##'
##' Builds the function \code{p -> community} that
##' \code{\link{community_parameter_sensitivity}},
##' \code{\link{community_continue_singularity}} and
##' \code{\link{community_selection_gradient_parameter_jacobian}} take, for a
##' community whose harness was built by \code{\link{harness_explicit}} (the
##' reference models \code{harness_dd99()}, \code{harness_jj12()}, ...). Editing
##' \code{community$harness$pars} does nothing, since the parameters are closed
##' into the model's functions when the harness is built; the map rebuilds the
##' harness instead. For any other model, or a parameter that is not one of the
##' harness's, write the function directly: it takes the parameter values and
##' returns a community of the same traits.
##'
##' @title Parameter map for an explicit harness
##' @param community A \code{community} with an explicit harness.
##' @param names Names of the harness parameters (\code{names(harness$pars)}) to
##' vary.
##' @return A function taking a numeric vector \code{p}, the named parameters
##' concatenated in order (a vector-valued parameter gives one entry per element),
##' and returning \code{community} with those parameters set to \code{p}: same
##' residents, bounds and controls, reset so that nothing solved under the old
##' parameters is carried over. The current values, named (\code{"sigma_K"}, or
##' \code{"sigma_K[1]"}, \code{"sigma_K[2]"} for a vector), are in
##' \code{attr(., "value")}.
##' @examples
##' comm <- community_start(bounds(x = c(-2, 2)), trait_scale = "linear",
##'                         harness = harness_dd99())
##' by_x0 <- community_parameter_map(comm, "x0")
##' attr(by_x0, "value")
##' community_solve_singularity(by_x0(0.5))$traits
##' @author Daniel Falster
##' @export
community_parameter_map <- function(community, names) {
  h <- community$harness
  if (!is.function(h$rebuild)) {
    stop("community_parameter_map() needs a harness built by harness_explicit(); ",
         "for another model write the function p -> community yourself")
  }
  if (!is.character(names) || length(names) < 1L || anyDuplicated(names)) {
    stop("names must name distinct parameters of the harness")
  }
  unknown <- setdiff(names, names(h$pars))
  if (length(unknown) > 0L) {
    stop("The harness has no parameter ", paste(unknown, collapse = ", "),
         "; it has ", paste(names(h$pars), collapse = ", "))
  }
  pars <- h$pars
  value <- lapply(pars[names], function(v) {
    if (!is.numeric(v) || length(v) < 1L) {
      stop("community_parameter_map() varies numeric parameters only")
    }
    as.numeric(v)
  })
  sizes <- lengths(value)
  labels <- unlist(Map(function(nm, s) if (s == 1L) nm else sprintf("%s[%d]", nm, seq_len(s)),
                       names, sizes), use.names = FALSE)
  value <- stats::setNames(unlist(value, use.names = FALSE), labels)
  last <- cumsum(sizes)

  map <- function(p) {
    p <- as.numeric(p)
    if (length(p) != length(value)) {
      stop("p must have ", length(value), " value(s), for ", paste(labels, collapse = ", "))
    }
    for (i in seq_along(names)) {
      pars[[names[i]]] <- p[(last[i] - sizes[i] + 1L):last[i]]
    }
    out <- community
    out$harness <- harness_fd(h$rebuild(pars))
    community_reset(out)
  }
  attr(map, "value") <- value
  map
}

##' Sensitivity of a singular strategy or coalition to model parameters.
##'
##' A singular point \eqn{x^*} is a root of the residents' selection gradients
##' \eqn{g(x, p) = 0}, so where the resident Jacobian
##' \eqn{J = \partial g / \partial x} is non-singular it moves smoothly with the
##' parameters, at
##' \deqn{dx^* / dp = -J^{-1} \partial g / \partial p.}
##' \eqn{J} is \code{\link{community_selection_gradient_jacobian}} (the
##' \code{mk x mk} Jacobian for a coalition of \code{m} residents with \code{k}
##' traits, through the environment they share, so a coalition's residents move
##' together) and \eqn{\partial g / \partial p} is
##' \code{\link{community_selection_gradient_parameter_jacobian}}. A singular
##' \eqn{J} is a fold or bifurcation of the singularity, where it does not move
##' smoothly with \eqn{p}, and is an error.
##'
##' Cost: \code{2mk} equilibrium solves for \eqn{J} and two per parameter for
##' \eqn{\partial g / \partial p}, all warm-started from the community.
##'
##' The sensitivity is in raw trait units whatever the community's trait scale
##' (on a log scale, \eqn{d \log x^* / dp} is it divided by \eqn{x^*}).
##' \code{community} is taken to be the model at \code{p}, as
##' \code{parameter(p)} builds it; it is not re-solved to check. As in
##' \code{\link{community_classify_singularity}}, a community marked as not
##' converged, or whose residents a Newton step on their gradients would move
##' by more than a thousandth of the bounds' width, is answered with a warning.
##'
##' @title Parameter sensitivity of a singular strategy or coalition
##' @param community A \code{community} at a singular strategy or coalition,
##' typically from \code{\link{community_solve_singularity}} on
##' \code{parameter(p)}.
##' @param parameter A function \code{p -> community}; see
##' \code{\link{community_parameter_map}}.
##' @param p The parameter values \code{community} was built at; defaults to
##' those \code{community_parameter_map()} records.
##' @param birth_rate Birth rates to start each equilibrium solve from, one per
##' resident; defaults to the residents' own.
##' @return An \code{mk} by \code{np} matrix of \eqn{dx^*/dp}, rows named as in
##' \code{\link{community_selection_gradient_jacobian}} and columns by
##' parameter, with \code{attr(., "jacobian")}, \code{attr(.,
##' "parameter_jacobian")} and the equilibrium solves in
##' \code{attr(., "evaluations")}.
##' @examples
##' comm <- community_start(bounds(x = c(-2, 2)), trait_scale = "linear",
##'                         harness = harness_dd99(sigma_C = 0.4))
##' by <- community_parameter_map(comm, c("sigma_C", "sigma_K"))
##' pair <- community_solve_singularity(comm, x0 = trait_matrix(c(-0.3, 0.6), "x"))
##' # the dimorphic pair spreads apart as either kernel widens
##' community_parameter_sensitivity(pair, by)
##' @author Daniel Falster
##' @export
community_parameter_sensitivity <- function(community, parameter,
                                            p = attr(parameter, "value"),
                                            birth_rate = NULL) {
  p <- check_parameter_values(parameter, p)
  if (nrow(community$traits) < 1L) {
    stop("community_parameter_sensitivity needs at least one resident (the ",
         "singular strategy or coalition)")
  }
  unconverged <- isFALSE(attr(community, "converged"))
  if (unconverged) {
    warning("community_parameter_sensitivity: the community is marked as not converged, ",
            "so this is the sensitivity of the point it stopped at, which need not be ",
            "a singularity")
  }
  J <- community_selection_gradient_jacobian(community, birth_rate = birth_rate)
  evaluations <- attr(J, "evaluations")
  reach <- singularity_newton_reach(community, community$traits,
                                    attr(J, "selection_gradient"), J)
  if (!unconverged && isTRUE(reach > singularity_reach_tol)) {
    warning(sprintf(paste0(
      "community_parameter_sensitivity: the residents are not at a singular point; ",
      "a Newton step on the selection gradients moves a trait by %s of the bounds' width"),
      signif(reach, 2)))
  }
  S <- singularity_sensitivity(community, J, parameter, p, birth_rate)
  attr(S, "evaluations") <- attr(S, "evaluations") + evaluations
  S
}

## dx*/dp from a resident Jacobian J already in hand (the classifier's, in a
## continuation), at the cost of the parameter Jacobian only.
singularity_sensitivity <- function(community, J, parameter, p, birth_rate = NULL) {
  G <- community_selection_gradient_parameter_jacobian(community, parameter, p,
                                                       birth_rate = birth_rate)
  evaluations <- attr(G, "evaluations")
  attr(G, "evaluations") <- NULL
  attr(J, "selection_gradient") <- NULL
  attr(J, "evaluations") <- NULL
  S <- tryCatch(-solve(J, G), error = function(e) {
    stop("community_parameter_sensitivity: the resident Jacobian is singular here, ",
         "a fold or bifurcation of the singularity, so it does not move smoothly ",
         "with the parameters (", conditionMessage(e), ")", call. = FALSE)
  })
  dimnames(S) <- dimnames(G)
  attr(S, "jacobian") <- J
  attr(S, "parameter_jacobian") <- G
  attr(S, "evaluations") <- evaluations
  S
}

##' Follow a singular strategy or coalition along one parameter.
##'
##' Natural-parameter continuation of a singular point through the values
##' \code{p}, in the order given. At each value the point is predicted from the
##' last by its sensitivity (\code{\link{community_parameter_sensitivity}}),
##' \eqn{x^*_i \approx x^*_{i-1} + (dx^*/dp)(p_i - p_{i-1})}, corrected by
##' \code{\link{community_solve_singularity}} started there with the last
##' equilibrium densities, and classified by
##' \code{\link{community_classify_singularity}}. The classifier's resident
##' Jacobian is also the one the next prediction uses and the one the next
##' correction starts from, so each point costs the corrector's residual
##' solves, the classification's and two more for
##' \eqn{\partial g / \partial p}. Where consecutive classifications differ the
##' path records the change (\code{changes}): a branching point becoming a CSS,
##' say, as a kernel widens.
##'
##' The path stops at the first value the corrector cannot reach --- it does not
##' converge, a coalition loses or merges residents, or the point leaves the
##' bounds --- with one warning naming the value and the reason, and keeps the
##' points before it. Such a stop usually brackets a fold or a bifurcation; a
##' finer \code{p} near it, or a path restarted from \code{community} at the last
##' point, says more. Failing at \code{p[1]} is an error.
##'
##' @title Continuation of a singular strategy or coalition along a parameter
##' @param community A \code{community} whose residents (and their densities, if
##' solved) start the path: the singular point at \code{p[1]}, or near it.
##' @param parameter A function from one parameter value to a \code{community};
##' see \code{\link{community_parameter_map}} (built with one name).
##' @param p The parameter values to visit, \code{p[1]} the starting one.
##' @param solver,tol,maxit Passed to \code{\link{community_solve_singularity}}
##' for each correction. The \code{"bracket"} solver ignores the prediction and
##' is not offered.
##' @param ... Passed to \code{\link{community_classify_singularity}}
##' (\code{tol}, \code{speeds}, \code{invasion_tol}).
##' @return An object of class \code{singularity_path}: a list with \code{p}
##' (the values reached), \code{traits}, \code{predicted} (the corrector's
##' starting points, the community's residents at \code{p[1]}) and
##' \code{sensitivity} (\eqn{dx^*/dp}), each a matrix with one row per value and
##' one column per resident trait, stacked as in
##' \code{\link{community_selection_gradient_jacobian}}; \code{birth_rate} (one
##' row per value, one column per resident); \code{classification} and the full
##' \code{classifications}; \code{changes}, a data frame of \code{p_before},
##' \code{p_after}, \code{from}, \code{to} for each change of classification;
##' \code{evaluations}, the equilibrium solves at each value (their sum in
##' \code{attr(., "evaluations")}); \code{stopped}, \code{NULL} or the value and
##' reason the path stopped at; and \code{community}, the community at the last
##' point, from which the path can be continued.
##' @examples
##' comm <- community_start(bounds(x = c(-2, 2)), trait_scale = "linear",
##'                         harness = harness_dd99(sigma_C = 0.7))
##' by_sigma_C <- community_parameter_map(comm, "sigma_C")
##' # a branching point becomes a CSS as the competition kernel outgrows the
##' # resource kernel (sigma_K = 1)
##' path <- community_continue_singularity(community_solve_singularity(comm),
##'                                        by_sigma_C, seq(0.7, 1.3, by = 0.2))
##' path
##' path$changes
##' @author Daniel Falster
##' @export
community_continue_singularity <- function(community, parameter, p,
                                           solver = c("nleqslv", "newton", "dfsane"),
                                           tol = 1e-6, maxit = 100, ...) {
  solver <- match.arg(solver)
  if (!is.function(parameter)) {
    stop("parameter must be a function p -> community; ",
         "community_parameter_map() builds one for an explicit harness")
  }
  if (!is.numeric(p) || length(p) < 1L || !all(is.finite(p))) {
    stop("p must be a vector of finite parameter values, the first the one ",
         "the path starts at")
  }
  trait_names <- community$trait_names
  k <- length(trait_names)
  m <- nrow(community$traits)
  if (m < 1L) {
    stop("community_continue_singularity needs at least one resident to start from")
  }
  labels <- resident_labels(trait_names, m)

  x <- as.numeric(community$traits)
  n <- as.numeric(community$birth_rate)
  if (length(n) != m || !all(is.finite(n) & n > 0)) n <- NULL
  S <- NULL
  J <- NULL
  points <- list()
  classifications <- list()
  stopped <- NULL
  last <- NULL

  for (i in seq_along(p)) {
    predicted <- if (i == 1L) x else x + as.numeric(S) * (p[i] - p[i - 1L])
    start <- parameter_community(parameter, p[i], community)
    ## the corrector's warnings are why the path stops, if it does; they are
    ## raised again if it does not
    caught <- list()
    sol <- withCallingHandlers(
      tryCatch(
        community_solve_singularity(start, x0 = trait_matrix(matrix(predicted, m, k), trait_names),
                                    solver = solver, tol = tol, maxit = maxit,
                                    birth_rate = n, edge_ok = FALSE,
                                    jacobian = J),
        error = function(e) e),
      warning = function(w) {
        caught[[length(caught) + 1L]] <<- w
        invokeRestart("muffleWarning")
      })
    if (inherits(sol, "error") || !isTRUE(attr(sol, "converged"))) {
      reason <- if (inherits(sol, "error")) conditionMessage(sol)
                else if (length(caught) > 0L) paste(unique(vapply(caught, conditionMessage, "")), collapse = "; ")
                else "the corrector did not converge"
      if (i == 1L) {
        stop("community_continue_singularity could not solve for the starting point at p = ",
             signif(p[1L], 6), ": ", reason)
      }
      stopped <- list(p = p[i], reason = reason)
      break
    }
    for (w in caught) warning(w)

    cl <- community_classify_singularity(sol, ...)
    J <- cl$jacobian
    S <- singularity_sensitivity(sol, J, parameter, p[i])
    x <- as.numeric(sol$traits)
    n <- as.numeric(sol$birth_rate)
    points[[i]] <- list(traits = x, predicted = predicted, sensitivity = as.numeric(S),
                        birth_rate = n,
                        evaluations = attr(sol, "evaluations") + cl$evaluations +
                          attr(S, "evaluations"))
    classifications[[i]] <- cl
    last <- sol
  }

  rows <- function(what, cols) {
    matrix(unlist(lapply(points, `[[`, what)), ncol = length(cols), byrow = TRUE,
           dimnames = list(NULL, cols))
  }
  reached <- p[seq_along(points)]
  classification <- vapply(classifications, `[[`, "", "classification")
  at <- which(classification[-1L] != classification[-length(classification)])
  evaluations <- vapply(points, `[[`, 0L, "evaluations")

  if (!is.null(stopped)) {
    warning(sprintf("community_continue_singularity stopped at p = %s: %s",
                    signif(stopped$p, 6), stopped$reason))
  }

  structure(list(
    p = reached,
    traits = rows("traits", labels),
    predicted = rows("predicted", labels),
    sensitivity = rows("sensitivity", labels),
    birth_rate = rows("birth_rate", if (m == 1L) "birth_rate" else sprintf("birth_rate[%d]", seq_len(m))),
    classification = classification,
    classifications = classifications,
    changes = data.frame(p_before = reached[at], p_after = reached[at + 1L],
                         from = classification[at], to = classification[at + 1L],
                         stringsAsFactors = FALSE),
    evaluations = evaluations,
    stopped = stopped,
    community = last),
    evaluations = sum(evaluations),
    class = "singularity_path")
}

##' @param x A \code{singularity_path} object.
##' @rdname community_continue_singularity
##' @export
print.singularity_path <- function(x, ...) {
  m <- length(x$community$birth_rate)
  cat(sprintf("<singularity_path: %d point%s, p from %s to %s%s>\n",
              length(x$p), if (length(x$p) == 1L) "" else "s",
              signif(x$p[1L], 6), signif(x$p[length(x$p)], 6),
              if (m > 1L) sprintf(" (coalition of %d)", m) else ""))
  runs <- rle(x$classification)
  end <- cumsum(runs$lengths)
  start <- end - runs$lengths + 1L
  cat("  classification:",
      paste(sprintf("%s (p %s to %s)", runs$values, signif(x$p[start], 6), signif(x$p[end], 6)),
            collapse = ", "), "\n")
  if (!is.null(x$stopped)) {
    cat(sprintf("  stopped at p = %s: %s\n", signif(x$stopped$p, 6), x$stopped$reason))
  }
  cat(sprintf("  equilibrium solves: %d\n", attr(x, "evaluations")))
  invisible(x)
}
