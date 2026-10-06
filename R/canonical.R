# The canonical equation of adaptive dynamics: trait values of every resident
# move up their selection gradient at a speed set by mutation supply, the
# community being re-solved to demographic equilibrium at every step. When
# selection stops, each resident is tested for evolutionary stability; one that
# sits at a fitness minimum is split in two and the dimorphic equation
# continues, so branching, the build-up of a polymorphic coalition and the
# extinction of residents it drives out all fall out of one integration.
#
# Each right-hand side evaluation is one equilibrium solve (warm-started from
# the last) plus the fitness gradient of every resident, so the stepper is
# chosen for how few evaluations it needs, not for accuracy in the trajectory:
# an explicit adaptive Bogacki-Shampine 3(2) by default, or a linearly implicit
# Rosenbrock step for the stiff case of many traits or species (#43), which
# takes the Jacobian of the right-hand side -- by finite differences over the
# gradient today, exactly once the model supplies equilibrium sensitivities.

##' Control the canonical equation.
##'
##' @title Canonical-equation settings
##' @param control A list of values to modify from the defaults:
##' \describe{
##'   \item{\code{rate}}{scalar on the right-hand side: mutation rate times
##'     mutational variance, in units of evolutionary time.}
##'   \item{\code{vcv}}{mutational covariance on the trait scale, a k by k
##'     matrix; \code{NULL} is the identity.}
##'   \item{\code{density}}{multiply each resident's speed by its equilibrium
##'     density relative to the initial total (mutation supply scales with
##'     population size).}
##'   \item{\code{stepper}}{\code{"rk23"} (explicit, adaptive) or
##'     \code{"rosenbrock"} (linearly implicit ROS2, adaptive, for stiff
##'     dynamics).}
##'   \item{\code{t_max}, \code{max_steps}}{limits on evolutionary time and
##'     accepted steps.}
##'   \item{\code{rtol}, \code{atol}}{step-size control on the trait scale.}
##'   \item{\code{dt0}, \code{dt_max}}{initial and largest step; \code{NULL}
##'     chooses \code{dt0} from the initial speed.}
##'   \item{\code{gradient_tol}}{selection has stopped when every resident's
##'     gradient on the trait scale is below this.}
##'   \item{\code{polish}, \code{polish_tol}}{once every gradient is below
##'     \code{polish_tol}, finish the approach to the stationary coalition by
##'     Newton on the residents' gradients instead of integrating the slow tail
##'     (the approach to a stable point is exponential, and its rate can be
##'     small).}
##'   \item{\code{branch}}{test stationary residents for evolutionary stability
##'     and split the invadable ones.}
##'   \item{\code{branch_distance}}{how far, as a fraction of the trait range on
##'     the trait scale, the two daughters of a split sit from the parent.}
##'   \item{\code{max_residents}}{stop rather than branch beyond this many.}
##'   \item{\code{extinct_fraction}}{a resident whose equilibrium density falls
##'     below this fraction of the total is removed.}
##'   \item{\code{classify_tol}}{eigenvalue magnitude below which the stability
##'     test is undecided and the resident is left alone.}
##' }
##' @return A list of the settings above.
##' @author Daniel Falster
##' @export
canonical_control <- function(control = NULL) {
  defaults <- list(
    rate = 1, vcv = NULL, density = TRUE,
    stepper = "rk23",
    t_max = 100, max_steps = 2000L,
    rtol = 1e-3, atol = 1e-5, dt0 = NULL, dt_max = Inf,
    gradient_tol = 1e-4,
    polish = TRUE, polish_tol = 1e-2,
    branch = TRUE, branch_distance = 0.02, max_residents = 8L,
    extinct_fraction = 1e-3,
    classify_tol = 1e-8
  )
  control <- as.list(control)
  extra <- setdiff(names(control), names(defaults))
  if (length(extra) > 0L) {
    stop("Unknown control parameters ", paste(extra, collapse = ", "))
  }
  ret <- modifyList(defaults, control, keep.null = TRUE)
  ret$stepper <- match.arg(ret$stepper, c("rk23", "rosenbrock"))
  for (nm in c("rate", "t_max", "rtol", "atol", "gradient_tol", "polish_tol", "branch_distance",
               "extinct_fraction", "classify_tol")) {
    v <- ret[[nm]]
    if (!is.numeric(v) || length(v) != 1L || !is.finite(v) || v <= 0) {
      stop(nm, " must be a single positive number")
    }
  }
  if (!is.null(ret$dt0) && (!is.numeric(ret$dt0) || ret$dt0 <= 0)) stop("dt0 must be positive or NULL")
  if (!is.numeric(ret$dt_max) || ret$dt_max <= 0) stop("dt_max must be positive")
  for (nm in c("max_steps", "max_residents")) {
    v <- ret[[nm]]
    if (!is.numeric(v) || length(v) != 1L || v < 1 || v != round(v)) stop(nm, " must be a positive whole number")
    ret[[nm]] <- as.integer(v)
  }
  if (!is.logical(ret$branch) || length(ret$branch) != 1L || is.na(ret$branch)) stop("branch must be TRUE or FALSE")
  if (!is.logical(ret$polish) || length(ret$polish) != 1L || is.na(ret$polish)) stop("polish must be TRUE or FALSE")
  if (!is.logical(ret$density) || length(ret$density) != 1L || is.na(ret$density)) stop("density must be TRUE or FALSE")
  if (!is.null(ret$vcv) && (!is.matrix(ret$vcv) || nrow(ret$vcv) != ncol(ret$vcv))) stop("vcv must be a square matrix")
  ret
}

## The right-hand side as a closure over the community: residents z (m x k on
## the trait scale) -> their speeds, with the solved community kept for warm
## starts and for the caller.
canonical_rhs <- function(base, tf, control, k, n_total0) {
  birth_rate <- NULL
  state <- NULL
  last <- NULL
  V <- if (is.null(control$vcv)) diag(1, k) else control$vcv
  evaluations <- 0L
  f <- function(z) {
    m <- length(z) / k
    x <- tf$inv(matrix(z, m, k))
    comm <- base |>
      community_add(trait_matrix(x, base$trait_names), birth_rate = birth_rate)
    comm$demography_state <- state
    comm <- community_demography(comm)
    evaluations <<- evaluations + 1L
    n <- as.numeric(comm$birth_rate)
    if (all(is.finite(n)) && all(n > 0)) birth_rate <<- n
    state <<- comm$demography_state
    last <<- comm
    g <- community_fitness_gradient(comm)                       # m x k, raw units
    g_z <- if (identical(tf$scale, "log")) g * x else g           # w.r.t. the trait scale
    weight <- if (control$density) n / n_total0 else rep(1, m)
    speed <- control$rate * weight * (g_z %*% V)
    list(f = as.numeric(speed), g_z = g_z, n = n, community = comm)
  }
  list(f = f, last = function() last, evaluations = function() evaluations,
       reset = function() { birth_rate <<- NULL; state <<- NULL })
}

## Finish the approach to a stationary coalition: Newton on the residents'
## gradients, g(z) = 0, with the equilibrium re-solved inside each evaluation.
## The Jacobian is finite differences across the residents -- one equilibrium
## solve per unknown -- which is the irreducible part until the model supplies
## equilibrium sensitivities (the gradient inside each evaluation is already
## exact where the harness provides it). Returns the new z, or NULL if Newton
## did not converge, in which case the integration simply carries on.
canonical_polish <- function(rhs, z, control) {
  residual <- function(zz) as.numeric(rhs(zz)$g_z)
  sol <- tryCatch(util_nlsolve(z, residual, tol = control$gradient_tol / 10, maxit = 50,
                               solver = "newton", require_converged = FALSE,
                               max_step = 0.1 * max(abs(z), 1)),
                  error = function(e) NULL)
  if (is.null(sol) || !isTRUE(attr(sol, "converged"))) return(NULL)
  as.numeric(sol)
}

## Step-size control shared by both steppers: error norm against the tolerance,
## accept at <= 1, and the next step from the usual power law.
canonical_step_control <- function(err, y, h, order, control) {
  scale <- control$atol + control$rtol * abs(y)
  norm <- max(abs(err) / scale)
  factor <- if (norm == 0) 5 else min(5, max(0.2, 0.9 * norm^(-1 / (order + 1))))
  list(accept = norm <= 1, h_next = min(h * factor, control$dt_max))
}

## Bogacki-Shampine 3(2), first-same-as-last: three right-hand sides per step.
step_rk23 <- function(rhs, y, h, k1, control) {
  k2 <- rhs(y + h / 2 * k1$f)
  k3 <- rhs(y + 3 * h / 4 * k2$f)
  y_new <- y + h * (2 / 9 * k1$f + 1 / 3 * k2$f + 4 / 9 * k3$f)
  k4 <- rhs(y_new)
  err <- h * (-5 / 72 * k1$f + 1 / 12 * k2$f + 1 / 9 * k3$f - 1 / 8 * k4$f)
  c(canonical_step_control(err, y_new, h, 2, control), list(y = y_new, k_end = k4))
}

## ROS2 (Verwer, Spee, Blom & Hundsdorfer 1999): two stages, order two,
## L-stable at gamma = 1 + 1/sqrt(2), with the linearly implicit Euler step as
## the embedded first-order estimate. One Jacobian and two right-hand sides per
## step; the Jacobian is of the right-hand side itself, by finite differences
## over the residents (one evaluation per unknown), until the model supplies
## equilibrium sensitivities.
step_rosenbrock <- function(rhs, y, h, k1, control) {
  n <- length(y)
  J <- util_fd_jacobian(function(z) rhs(z)$f, y, k1$f, h = 1e-5 * pmax(abs(y), 1))
  gamma <- 1 + 1 / sqrt(2)
  W <- diag(1, n) - gamma * h * J
  s1 <- solve(W, k1$f)
  f2 <- rhs(y + h * s1)
  s2 <- solve(W, f2$f - 2 * s1)
  y_new <- y + h * (1.5 * s1 + 0.5 * s2)
  err <- h * (0.5 * s1 + 0.5 * s2)
  k_end <- rhs(y_new)
  c(canonical_step_control(err, y_new, h, 1, control), list(y = y_new, k_end = k_end))
}

##' Integrate the canonical equation of adaptive dynamics.
##'
##' Every resident's trait moves up its selection gradient,
##' \deqn{\frac{dz_i}{dt} = \mathrm{rate}\; \tilde n_i\; V\, \nabla_z s(z_i),}
##' on the community's trait scale \eqn{z}, with the community re-solved to
##' demographic equilibrium at every evaluation; \eqn{\tilde n_i} is the
##' resident's equilibrium density relative to the initial total and \eqn{V}
##' the mutational covariance (see \code{\link{canonical_control}}). When
##' every gradient has vanished the residents are classified: one whose
##' fitness Hessian in its own direction has a positive eigenvalue sits at a
##' fitness minimum and is split into two daughters along that eigenvector,
##' after which the polymorphic equation continues; when all are at fitness
##' maxima the coalition is evolutionarily stable and the integration stops.
##' A resident whose density falls below \code{extinct_fraction} of the total
##' is removed.
##'
##' @title Canonical equation with branching
##' @param community A \code{community}; its residents (or \code{x0}) are the
##' starting strategies.
##' @param x0 Starting trait values, one row per resident, if the community
##' carries none.
##' @param control See \code{\link{canonical_control}}.
##' @return A \code{canonical_equation} object: \code{trajectory} (a tibble of
##' \code{time}, \code{lineage}, the trait columns, \code{density} and the
##' gradient columns \code{gradient_<trait>} on the trait scale, one row per
##' resident per accepted step), \code{events} (\code{time}, \code{event} ---
##' \code{"branch"}, \code{"extinct"}, \code{"polish"}, \code{"stable"},
##' \code{"t_max"}, \code{"max_steps"}, \code{"max_residents"} --- and
##' \code{lineage}),
##' \code{community} (the final solved community), \code{evaluations} (how
##' many equilibrium solves it cost), \code{outcome} and the control.
##' @author Daniel Falster
##' @export
community_canonical_equation <- function(community, x0 = NULL, control = canonical_control()) {
  control <- canonical_control(control)
  trait_names <- community$trait_names
  k <- length(trait_names)
  tf <- community_trait_transform(community)
  range_z <- as.numeric(diff(t(tf$fwd(community$bounds))))
  if (any(!is.finite(range_z))) stop("community_canonical_equation needs finite bounds")
  if (!is.null(control$vcv) && nrow(control$vcv) != k) stop("vcv must be ", k, " by ", k)

  if (is.null(x0)) {
    if (nrow(community$traits) == 0L) stop("Give starting traits in x0 or as residents on the community")
    x0 <- community$traits
  }
  x0 <- if (is.matrix(x0)) x0 else matrix(x0, ncol = k)
  if (ncol(x0) != k) stop("x0 must have ", k, " column(s)")
  m <- nrow(x0)
  base <- community_clear_residents(community)

  ## initial solve sets the density scale and the first right-hand side
  first <- base |> community_add(trait_matrix(x0, trait_names)) |> community_demography()
  n_total0 <- sum(as.numeric(first$birth_rate))
  if (!is.finite(n_total0) || n_total0 <= 0) stop("The starting community has no positive equilibrium density")
  rhs <- canonical_rhs(base, tf, control, k, n_total0)
  step <- switch(control$stepper, rk23 = step_rk23, rosenbrock = step_rosenbrock)

  z <- as.numeric(tf$fwd(x0))
  lineage <- seq_len(m)
  next_lineage <- m + 1L
  t <- 0
  k1 <- rhs$f(z)
  h <- if (is.null(control$dt0)) min(control$dt_max, 0.01 * max(range_z) / max(max(abs(k1$f)), 1e-12)) else control$dt0

  rows <- list()
  events <- list()
  ## one record per time: an event at the same time as the step just taken
  ## (a polish, a split, an extinction) replaces that step's rows
  record <- function(t, z, k) {
    x <- tf$inv(matrix(z, length(lineage), length(trait_names)))
    tr <- as.data.frame(x); names(tr) <- trait_names
    gz <- as.data.frame(k$g_z); names(gz) <- paste0("gradient_", trait_names)
    rows <<- Filter(function(r) r$time[1] != t, rows)
    rows[[length(rows) + 1L]] <<- cbind(data.frame(time = t, lineage = lineage), tr, density = k$n, gz)
  }
  event <- function(t, what, who = NA_integer_) {
    events[[length(events) + 1L]] <<- data.frame(time = t, event = what, lineage = who)
  }
  record(t, z, k1)
  outcome <- NULL
  steps <- 0L

  while (is.null(outcome)) {
    if (t >= control$t_max) { outcome <- "t_max"; break }
    if (steps >= control$max_steps) { outcome <- "max_steps"; break }
    h <- min(h, control$t_max - t)
    attempt <- step(rhs$f, z, h, k1, control)
    if (!attempt$accept) { h <- attempt$h_next; next }
    steps <- steps + 1L
    t <- t + h
    z <- attempt$y
    k1 <- attempt$k_end
    h <- attempt$h_next

    ## extinction: a resident the others have driven out stops counting
    mm <- length(lineage)
    frac <- k1$n / sum(k1$n)
    gone <- which(frac < control$extinct_fraction)
    if (length(gone) > 0L && length(gone) < mm) {
      for (i in gone) event(t, "extinct", lineage[i])
      keep <- setdiff(seq_len(mm), gone)
      zm <- matrix(z, mm, k)[keep, , drop = FALSE]
      z <- as.numeric(zm)
      lineage <- lineage[keep]
      rhs$reset()
      k1 <- rhs$f(z)
    }
    record(t, z, k1)

    ## near stationarity, finish with Newton on the coalition's gradients rather
    ## than integrating the slow exponential tail
    if (control$polish && max(abs(k1$g_z)) < control$polish_tol &&
        max(abs(k1$g_z)) >= control$gradient_tol) {
      polished <- canonical_polish(rhs$f, z, control)
      if (!is.null(polished)) {
        z <- polished
        k1 <- rhs$f(z)
        event(t, "polish")
        record(t, z, k1)
      }
    }

    ## stationary: selection has stopped on every resident
    if (max(abs(k1$g_z)) < control$gradient_tol) {
      if (!control$branch) { outcome <- "stable"; break }
      comm <- k1$community
      mm <- length(lineage)
      zm <- matrix(z, mm, k)
      xm <- tf$inv(zm)
      split <- list()
      for (i in seq_len(mm)) {
        H <- community_fitness_hessian(comm, xm[i, , drop = FALSE])
        if (identical(tf$scale, "log")) {
          ## curvature on the trait scale: x^2 H + x g, the gradient term ~0 here
          H <- H * outer(xm[i, ], xm[i, ])
        }
        e <- eigen((H + t(H)) / 2, symmetric = TRUE)
        if (max(e$values) > control$classify_tol) {
          split[[length(split) + 1L]] <- list(i = i, direction = e$vectors[, which.max(e$values)])
        }
      }
      if (length(split) == 0L) { outcome <- "stable"; break }
      if (mm + length(split) > control$max_residents) { outcome <- "max_residents"; break }
      new_z <- zm
      new_lineage <- lineage
      for (s in split) {
        d <- s$direction / sqrt(sum(s$direction^2)) * control$branch_distance * range_z
        new_z[s$i, ] <- zm[s$i, ] - d
        new_z <- rbind(new_z, zm[s$i, ] + d)
        new_lineage <- c(new_lineage, next_lineage)
        event(t, "branch", lineage[s$i])
        next_lineage <- next_lineage + 1L
      }
      z <- as.numeric(new_z)
      lineage <- new_lineage
      rhs$reset()
      k1 <- rhs$f(z)
      record(t, z, k1)
      h <- if (is.null(control$dt0)) 0.01 * max(range_z) / max(max(abs(k1$f)), 1e-12) else control$dt0
    }
  }
  event(t, outcome)

  structure(list(
    trajectory = tibble::as_tibble(do.call(rbind, rows)),
    events = tibble::as_tibble(do.call(rbind, events)),
    community = rhs$last(),
    evaluations = rhs$evaluations(),
    steps = steps,
    outcome = outcome,
    trait_names = trait_names,
    trait_scale = tf$scale,
    control = control
  ), class = "canonical_equation")
}

##' @export
print.canonical_equation <- function(x, ...) {
  n_branch <- sum(x$events$event == "branch")
  final <- x$trajectory[x$trajectory$time == max(x$trajectory$time), ]
  cat(sprintf("<canonical_equation: %s after %d steps (%d equilibrium solves), t = %.3g; %d branching event%s, %d resident%s at the end>\n",
              x$outcome, x$steps, x$evaluations, max(x$trajectory$time),
              n_branch, if (n_branch == 1L) "" else "s", nrow(final), if (nrow(final) == 1L) "" else "s"))
  invisible(x)
}

##' Plot trait trajectories from the canonical equation.
##'
##' Each lineage is a line of trait value against evolutionary time; a
##' branching event is where one line becomes two. With several traits, one
##' panel per trait.
##'
##' @title Plot a canonical-equation trajectory
##' @param x A \code{canonical_equation}.
##' @param object A \code{canonical_equation}, for \code{autoplot}.
##' @param ... Ignored.
##' @return A \code{ggplot} object.
##' @author Daniel Falster
##' @export
plot.canonical_equation <- function(x, ...) {
  tr <- x$trajectory
  long <- do.call(rbind, lapply(x$trait_names, function(nm) {
    data.frame(time = tr$time, lineage = tr$lineage, trait = nm, value = tr[[nm]])
  }))
  branches <- x$events[x$events$event == "branch", ]
  p <- ggplot2::ggplot(long, ggplot2::aes(x = .data[["time"]], y = .data[["value"]],
                                          group = .data[["lineage"]])) +
    ggplot2::geom_line(linewidth = 0.6)
  if (nrow(branches) > 0L) {
    p <- p + ggplot2::geom_vline(xintercept = branches$time, linetype = "dotted", colour = "grey50")
  }
  if (length(x$trait_names) > 1L) {
    p <- p + ggplot2::facet_wrap(~ trait, scales = "free_y")
  }
  p <- p + ggplot2::xlab("evolutionary time") +
    ggplot2::ylab(if (length(x$trait_names) == 1L) x$trait_names else "trait value") +
    ggplot2::theme_classic() +
    ggplot2::theme(text = ggplot2::element_text(size = 16))
  if (identical(x$trait_scale, "log")) p <- p + ggplot2::scale_y_log10()
  p
}

##' @rdname plot.canonical_equation
##' @export
autoplot.canonical_equation <- function(object, ...) plot.canonical_equation(object, ...)
