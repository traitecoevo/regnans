# Derivatives of invasion fitness and of the selection gradient.
#
# Every method that needs a derivative -- the selection gradient, the mutant
# Hessian behind the ESS test, the resident Jacobian behind convergence
# stability -- obtains it here and never from a stencil of its own. This is the
# one place that decides where a derivative comes from, so that a model able to
# supply exact derivatives (automatic differentiation in plant, closed forms in
# the reference models) can do so without any consumer changing.
#
# A harness that can differentiate its fitness sets closures in
# community$fitness_derivatives when it builds community$fitness_function;
# harness_fd() fills whatever is missing with the finite-difference machinery
# in util_gradient.R and records the source of each. Consumers call
# community_fitness_gradient() and friends and never branch on the source.
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

## The derivative quantities a harness can provide in the mutant direction.
## Each names a closure `function(y, control)` in community$fitness_derivatives.
fitness_derivative_names <- c("fitness_gradient", "fitness_hessian")

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

community_fitness_derivative <- function(community, what) {
  fn <- community$fitness_derivatives[[sub("^fitness_", "", what)]]
  if (!is.function(fn)) {
    community_fitness_function(community)
    stop("The community has no ", what, "; its harness was not passed through ",
         "harness_fd() -- see community_start()")
  }
  fn
}

## Fill in finite-difference closures for any derivative the harness did not
## supply, and record where each one comes from. Called by harness_fd() right
## after the harness has built the fitness function.
fd_fill_fitness_derivatives <- function(community) {
  f <- community$fitness_function
  if (!is.function(f)) {
    return(community)
  }
  der <- community$fitness_derivatives
  if (is.null(der)) der <- list()
  source <- character(0)
  if (!is.function(der$gradient)) {
    der$gradient <- function(y, control) fd_fitness_gradient(f, y, control)
    source["fitness_gradient"] <- "finite difference"
  } else {
    source["fitness_gradient"] <- "model"
  }
  if (!is.function(der$hessian)) {
    der$hessian <- function(y, control) fd_fitness_hessian(f, y, control)
    source["fitness_hessian"] <- "finite difference"
  } else {
    source["fitness_hessian"] <- "model"
  }
  der$source <- source
  community$fitness_derivatives <- der
  community
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
  trait_names <- community$trait_names
  k <- length(trait_names)
  y <- derivative_points(y, trait_names)
  ctrl <- community_derivative_control(community)

  g <- community_fitness_derivative(community, "fitness_gradient")(y, ctrl)
  value <- attr(g, "value")
  if (!is.matrix(g)) {
    g <- matrix(as.numeric(g), ncol = k)
  }
  if (!is.numeric(g) || !identical(dim(g), c(nrow(y), k))) {
    stop("fitness_gradient must return a ", nrow(y), " by ", k, " matrix")
  }
  if (is.null(value)) {
    value <- as.numeric(community_fitness_function(community)(y))
  }
  g <- matrix(as.numeric(g), nrow(y), k)
  dimnames(g) <- list(NULL, trait_names)
  attr(g, "value") <- value
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
  k <- length(community$trait_names)
  ctrl <- community_derivative_control(community)

  H <- community_fitness_derivative(community, "fitness_hessian")(y, ctrl)
  if (!is.numeric(H) || length(H) != k * k) {
    stop("fitness_hessian must return a ", k, " by ", k, " matrix")
  }
  H <- matrix(as.numeric(H), k, k)
  dimnames(H) <- list(community$trait_names, community$trait_names)
  H
}

fd_fitness_hessian <- function(f, y, ctrl) {
  util_hessian(f, y[1, ], d = ctrl$d_second, eps = ctrl$eps_second,
               r = ctrl$r_hessian)
}

##' Jacobian of the selection gradient with respect to the resident traits.
##'
##' How the selection gradients of the residents change as the residents move,
##' each re-solved to demographic equilibrium: for one resident the \code{k}
##' by \code{k} matrix \eqn{J_{ij} = d g_i / d x_j}; for \code{m} residents the
##' \code{mk} by \code{mk} matrix of every resident's gradient with respect to
##' every resident's traits, its gradients and traits stacked as the trait
##' matrix is stored (trait by trait, residents within each). The
##' off-diagonal blocks are how a resident's selection responds to the others
##' through the environment they make. Eigenvalues with negative real parts
##' mean selection carries nearby residents towards the point ---
##' convergence stability. Each stencil point solves the perturbed community to
##' demographic equilibrium, so this is the expensive derivative: \code{2mk}
##' equilibrium solves per Richardson level, finite differences across the
##' residents until the model supplies equilibrium sensitivities.
##'
##' @title Resident Jacobian of the selection gradient
##' @param community A \code{community} with one or more residents, no two of
##' them within the finite-difference step of each other (an error: the
##' stencil cannot tell them apart).
##' @param birth_rate Birth rates to start each equilibrium solve from, one per
##' resident; defaults to the residents' own, with each solve then
##' warm-starting from the last.
##' @return An \code{mk} by \code{mk} matrix, rows and columns named by trait
##' (and \code{[i]} for resident \code{i} when there are several), with the
##' stacked selection gradient at the residents in
##' \code{attr(., "selection_gradient")} and the equilibrium solves it cost in
##' \code{attr(., "evaluations")}.
##' @author Daniel Falster
##' @export
community_selection_gradient_jacobian <- function(community, birth_rate = NULL) {
  trait_names <- community$trait_names
  m <- nrow(community$traits)
  if (m < 1L) {
    stop("community_selection_gradient_jacobian needs at least one resident")
  }
  x <- as.numeric(community$traits)
  gradient <- singularity_gradient_fn(community, m, birth_rate = birth_rate)
  g0 <- gradient(x)
  ctrl <- community_derivative_control(community)
  J <- resident_jacobian(gradient, x, ctrl)
  ## a stencil point that merged two residents was refused, leaving the
  ## Jacobian across them undefined (the solvers back off from such a point;
  ## a caller asking for the Jacobian there has residents too close to tell
  ## apart)
  if (attr(gradient, "refused")() > 0L) {
    near <- residents_merged(x, m, length(trait_names), ctrl, reach = 2)
    stop("community_selection_gradient_jacobian: ",
         if (nrow(near) > 0L)
           paste(sprintf("residents %d and %d", near[, 1], near[, 2]), collapse = "; ")
         else "two residents",
         " are within the finite-difference step of each other, so the Jacobian across them is undefined")
  }
  labels <- resident_labels(trait_names, m)
  dimnames(J) <- list(labels, labels)
  names(g0) <- labels
  attr(J, "selection_gradient") <- g0
  attr(J, "evaluations") <- attr(gradient, "evaluations")()
  J
}

## The resident Jacobian of a stacked gradient closure (singularity_gradient_fn)
## at x, raw trait units. Solvers that already hold the closure call this
## directly, so the stencil shares its warm starts and its count of solves.
resident_jacobian <- function(gradient, x, ctrl) {
  util_jacobian(gradient, x, d = ctrl$d_second, eps = ctrl$eps_second,
                r = ctrl$r_jacobian)
}

##' Check a harness's derivatives against finite differences.
##'
##' Every derivative a model supplies must agree with a finite difference of
##' its own fitness function; this is the contract a harness author signs by
##' advertising a derivative (\code{harness$provides}). The check evaluates each
##' advertised quantity at \code{n_points} well-spread points inside the
##' community's bounds (on its trait scale) and compares every entry with a
##' three-level Richardson finite difference. A harness that advertises nothing
##' has nothing to check and returns an empty table.
##'
##' Call it from a package's own tests with a tolerance suited to the model:
##' the reference models agree to machine precision; a model whose fitness
##' comes from an adaptive ODE integration needs something like
##' \code{tol_rel = 1e-2}. The absolute tolerance exists for entries that are
##' close to zero, where the finite-difference reference itself is only good to
##' its roundoff floor --- about \code{1e-7} for a second difference of an
##' order-one fitness with the default steps.
##'
##' @title Verify model-supplied derivatives
##' @param community A \code{community}; solved to demographic equilibrium if it
##' has residents (an empty community is solved here).
##' @param n_points Number of trait points to test at.
##' @param tol_rel,tol_abs An entry passes when
##' \code{abs(model - fd) <= tol_abs + tol_rel * abs(fd)}.
##' @return Invisibly, a data frame with one row per entry checked:
##' \code{quantity}, \code{point}, \code{entry}, \code{model}, \code{fd},
##' \code{abs_err}, \code{rel_err}, \code{pass}. Any failure is an error of
##' class \code{derivative_check_error} carrying the table as \code{$results}.
##' @author Daniel Falster
##' @export
harness_check_derivatives <- function(community, n_points = 5, tol_rel = 1e-4,
                                      tol_abs = 1e-6) {
  provides <- intersect(community$harness$provides, fitness_derivative_names)
  trait_names <- community$trait_names
  k <- length(trait_names)
  empty <- data.frame(quantity = character(0), point = integer(0),
                      entry = character(0), model = numeric(0), fd = numeric(0),
                      abs_err = numeric(0), rel_err = numeric(0),
                      pass = logical(0), stringsAsFactors = FALSE)
  if (length(provides) == 0L) {
    return(invisible(empty))
  }
  if (!is.function(community$fitness_function)) {
    community <- community_demography(community)
  }
  f <- community_fitness_function(community)
  der <- community$fitness_derivatives
  if (!identical(unname(der$source[provides]), rep("model", length(provides)))) {
    stop("The harness advertises ", paste(provides, collapse = ", "),
         " but the community carries no model-supplied closure for it")
  }
  ctrl <- community_derivative_control(community)
  reference <- modifyList(ctrl, list(r_gradient = 3L, r_hessian = 3L))
  y <- derivative_check_points(community, n_points)

  rows <- list()
  add <- function(quantity, point, entry, model, fd) {
    rows[[length(rows) + 1L]] <<- data.frame(
      quantity = quantity, point = point, entry = entry,
      model = model, fd = fd, stringsAsFactors = FALSE)
  }
  if ("fitness_gradient" %in% provides) {
    g_model <- community_fitness_gradient(community, y)
    g_fd <- fd_fitness_gradient(f, y, reference)
    for (i in seq_len(nrow(y))) {
      for (d in seq_len(k)) {
        add("fitness_gradient", i, trait_names[d], g_model[i, d], g_fd[i, d])
      }
    }
  }
  if ("fitness_hessian" %in% provides) {
    for (i in seq_len(nrow(y))) {
      H_model <- community_fitness_hessian(community, y[i, , drop = FALSE])
      H_fd <- fd_fitness_hessian(f, y[i, , drop = FALSE], reference)
      for (d in seq_len(k)) {
        for (e in seq_len(k)) {
          add("fitness_hessian", i, paste(trait_names[d], trait_names[e], sep = ":"),
              H_model[d, e], H_fd[d, e])
        }
      }
    }
  }
  results <- do.call(rbind, rows)
  results$abs_err <- abs(results$model - results$fd)
  results$rel_err <- results$abs_err / pmax(abs(results$fd), .Machine$double.eps)
  results$pass <- results$abs_err <= tol_abs + tol_rel * abs(results$fd)

  if (!all(results$pass)) {
    worst <- results[which.max(results$abs_err / (tol_abs + tol_rel * abs(results$fd))), ]
    stop(structure(class = c("derivative_check_error", "error", "condition"), list(
      message = sprintf(
        "%d of %d derivative entries disagree with finite differences; worst is %s[%s] at point %d: model %g, finite difference %g",
        sum(!results$pass), nrow(results), worst$quantity, worst$entry,
        worst$point, worst$model, worst$fd),
      call = NULL, results = results)))
  }
  invisible(results)
}

## Well-spread, deterministic points inside the bounds on the trait scale (a
## Kronecker sequence with the generalised golden ratio), kept off the edges.
derivative_check_points <- function(community, n) {
  trait_names <- community$trait_names
  k <- length(trait_names)
  tf <- community_trait_transform(community)
  lo <- tf$fwd(community$bounds[, 1])
  hi <- tf$fwd(community$bounds[, 2])
  if (any(!is.finite(lo) | !is.finite(hi))) {
    stop("harness_check_derivatives needs finite bounds on the trait scale")
  }
  phi <- 2
  for (i in seq_len(30)) phi <- (1 + phi)^(1 / (k + 1))
  frac <- (0.5 + outer(seq_len(n), phi^-seq_len(k))) %% 1
  z <- sweep(sweep(0.05 + 0.9 * frac, 2, hi - lo, "*"), 2, lo, "+")
  y <- tf$inv(z)
  dim(y) <- c(n, k)
  colnames(y) <- trait_names
  y
}
