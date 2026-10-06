# The canonical equation of adaptive dynamics: trait values of every resident
# move up their selection gradient at a speed set by mutation supply, the
# community being re-solved to demographic equilibrium at every step. When
# selection stops, each resident is tested for evolutionary stability; one that
# sits at a fitness minimum branches -- after the waiting time it takes for a
# mutation to arise that can both invade and coexist, which depends on the
# mutational kernel and on how disruptive the fitness minimum is -- and the
# polymorphic equation continues. New phenotypes may also arrive from outside
# at a given rate, so local mutation and immigration can be run together or
# apart. Branching, the build-up of a coalition and the extinction of residents
# it drives out all fall out of one integration, and the record it keeps lets
# the community and its fitness landscape be rebuilt at any time along the way.
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
##'     Newton on the residents' gradients instead of integrating it: near a
##'     stable coalition the dynamics are stiff, and an explicit stepper held at
##'     its stability limit jitters about the point without reaching
##'     \code{gradient_tol}.}
##'   \item{\code{branch}}{what happens to a stationary resident at a fitness
##'     minimum. \code{"expected"} (default): a mutant that can invade and
##'     coexist arises at the rate
##'     \eqn{\lambda = \mathrm{rate}\, \tilde n \int M(\delta)\, p_{est}(\delta)\,
##'     \mathbf{1}[\mathrm{coexist}(\delta)]\, d\delta} along the branching
##'     direction, with \eqn{M} the mutational kernel of width
##'     \code{mutation_sd}, \eqn{p_{est} = \min(1, 2s)} the establishment
##'     probability of a rare mutant with invasion fitness \eqn{s}, and
##'     coexistence meaning the resident can invade the mutant back; time
##'     advances by \eqn{1/\lambda} and the mutant appears at the expected
##'     distance of a successful mutation. \code{"stochastic"}: the waiting time
##'     is exponential with that rate and the mutant's distance is drawn.
##'     \code{"immediate"}: split at once, \code{branch_distance} either side.
##'     \code{"none"}: stop at the first stationary coalition.}
##'   \item{\code{mutation_sd}}{standard deviation of the mutational kernel on
##'     the trait scale, one per trait; \code{NULL} takes it from \code{vcv}
##'     where given, else 2\% of the trait range.}
##'   \item{\code{branch_nodes}}{Gauss--Hermite nodes for the branching-rate
##'     integral; each costs one equilibrium solve.}
##'   \item{\code{branch_distance}}{for \code{"immediate"}: how far, as a
##'     fraction of the trait range on the trait scale, the two daughters sit
##'     from the parent.}
##'   \item{\code{immigration}}{\code{NULL}, or a list with \code{rate}
##'     (arrivals per unit evolutionary time, exponential waiting times),
##'     \code{pool} (where phenotypes come from: \code{NULL} for uniform on the
##'     trait scale within the bounds, a matrix of phenotypes to sample rows
##'     from, or a function of \code{n} returning that many rows) and
##'     \code{establish} (\code{"deterministic"}: an arrival with positive
##'     invasion fitness establishes; \code{"stochastic"}: with probability
##'     \eqn{\min(1, 2s)}).}
##'   \item{\code{seed}}{random seed for the stochastic parts; the caller's RNG
##'     state is restored afterwards.}
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
    branch = "expected", mutation_sd = NULL, branch_nodes = 12L, branch_distance = 0.02,
    max_residents = 8L, immigration = NULL, seed = NULL,
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
  if (is.logical(ret$branch) && length(ret$branch) == 1L && !is.na(ret$branch)) {
    ret$branch <- if (ret$branch) "immediate" else "none"
  }
  ret$branch <- match.arg(ret$branch, c("expected", "stochastic", "immediate", "none"))
  if (!is.null(ret$mutation_sd) && (!is.numeric(ret$mutation_sd) || any(ret$mutation_sd <= 0))) stop("mutation_sd must be positive")
  if (!is.numeric(ret$branch_nodes) || ret$branch_nodes < 2 || ret$branch_nodes != round(ret$branch_nodes)) stop("branch_nodes must be a whole number of at least 2")
  ret$branch_nodes <- as.integer(ret$branch_nodes)
  if (!is.null(ret$immigration)) {
    im <- ret$immigration
    if (!is.list(im) || is.null(im$rate) || !is.numeric(im$rate) || im$rate < 0) stop("immigration must be a list with a non-negative rate")
    im$establish <- match.arg(if (is.null(im$establish)) "deterministic" else im$establish, c("deterministic", "stochastic"))
    if (!is.null(im$pool) && !is.matrix(im$pool) && !is.function(im$pool)) stop("immigration$pool must be NULL, a matrix or a function")
    ret$immigration <- im
  }
  if (!is.null(ret$seed) && (!is.numeric(ret$seed) || length(ret$seed) != 1L)) stop("seed must be a single number")
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
    speed <- 0.5 * control$rate * weight * (g_z %*% V)
    list(f = as.numeric(speed), g_z = g_z, n = n, community = comm)
  }
  list(f = f, last = function() last, evaluations = function() evaluations,
       reset = function() { birth_rate <<- NULL; state <<- NULL })
}

## Gauss-Hermite nodes and weights for int exp(-x^2) f(x) dx, from the Jacobi
## matrix (Golub & Welsch), so no dependency is needed for a dozen nodes.
gauss_hermite <- function(n) {
  J <- matrix(0, n, n)
  off <- sqrt(seq_len(n - 1L) / 2)
  J[cbind(seq_len(n - 1L), seq_len(n - 1L) + 1L)] <- off
  J[cbind(seq_len(n - 1L) + 1L, seq_len(n - 1L))] <- off
  e <- eigen(J, symmetric = TRUE)
  list(x = e$values, w = sqrt(pi) * e$vectors[1, ]^2)
}

## Run expr with a given seed, leaving the caller's random stream as it was.
with_seed <- function(seed, expr) {
  if (is.null(seed)) return(expr)
  had <- exists(".Random.seed", envir = globalenv(), inherits = FALSE)
  old <- if (had) get(".Random.seed", envir = globalenv()) else NULL
  on.exit(if (had) assign(".Random.seed", old, envir = globalenv()) else
            if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) rm(".Random.seed", envir = globalenv()))
  set.seed(seed)
  expr
}

## The rate at which resident i, stationary at a fitness minimum, acquires a
## mutant that can invade and coexist with it, and where such a mutant lands.
## Integrates along the branching direction v (unit, trait scale) over the
## mutational kernel of width sd_v: at each node the mutant's invasion fitness
## gives its establishment probability, and one equilibrium solve with the
## mutant in the resident's place says whether the resident can invade it back
## (a protected dimorphism). Returns the rate, and the nodes with their
## success weights so a distance can be taken (expected or drawn).
canonical_branch_rate <- function(base, comm, tf, control, zm, i, v, sd_v, weight_i) {
  k <- ncol(zm)
  gh <- gauss_hermite(control$branch_nodes)
  delta <- sqrt(2) * sd_v * gh$x
  wq <- gh$w / sqrt(pi)
  success <- numeric(length(delta))
  f <- community_fitness_function(comm)
  for (j in seq_along(delta)) {
    if (abs(delta[j]) < 1e-12) next
    z_mut <- zm[i, ] + delta[j] * v
    x_mut <- tf$inv(matrix(z_mut, 1, k))
    s_mut <- as.numeric(f(x_mut))
    if (!is.finite(s_mut) || s_mut <= 0) next
    p_est <- min(1, 2 * s_mut)
    ## can the resident invade the mutant's community? (mutant in its place)
    z_other <- zm; z_other[i, ] <- z_mut
    other <- tryCatch(
      base |> community_add(trait_matrix(tf$inv(z_other), base$trait_names)) |> community_demography(),
      error = function(e) NULL)
    if (is.null(other)) next
    s_back <- as.numeric(community_fitness_function(other)(tf$inv(matrix(zm[i, ], 1, k))))
    if (!is.finite(s_back) || s_back <= 0) next
    success[j] <- p_est
  }
  mass <- wq * success
  list(rate = control$rate * weight_i * sum(mass), delta = delta, mass = mass)
}

## A phenotype from the immigration pool, on the trait scale.
immigrant_phenotype <- function(control, tf, bounds_z, k) {
  pool <- control$immigration$pool
  if (is.null(pool)) {
    return(bounds_z[, 1] + stats::runif(k) * (bounds_z[, 2] - bounds_z[, 1]))
  }
  x <- if (is.function(pool)) pool(1L) else pool[sample.int(nrow(pool), 1L), , drop = FALSE]
  as.numeric(tf$fwd(matrix(as.numeric(x), 1, k)))
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
##' \deqn{\frac{dz_i}{dt} = \tfrac{1}{2}\,\mathrm{rate}\; \tilde n_i\; V\, \nabla_z s(z_i),}
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
##' \code{"branch"}, \code{"unbranchable"}, \code{"immigrant"},
##' \code{"extinct"}, \code{"polish"}, \code{"stable"}, \code{"t_max"},
##' \code{"max_steps"}, \code{"max_residents"} --- and \code{lineage}),
##' \code{community} (the final solved community), \code{evaluations} (how
##' many equilibrium solves it cost), \code{immigration_attempts},
##' \code{outcome} and the control. \code{\link{canonical_community}} rebuilds
##' the community at any recorded time.
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
  mutation_sd <- control$mutation_sd
  if (is.null(mutation_sd)) {
    mutation_sd <- if (!is.null(control$vcv)) sqrt(diag(control$vcv)) else 0.02 * range_z
  }
  mutation_sd <- rep_len(mutation_sd, k)
  bounds_z <- tf$fwd(community$bounds)
  with_seed(control$seed, canonical_integrate(community, x0, control, tf, range_z, mutation_sd, bounds_z))
}

canonical_integrate <- function(community, x0, control, tf, range_z, mutation_sd, bounds_z) {
  trait_names <- community$trait_names
  k <- length(trait_names)

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
  immigration_attempts <- 0L
  next_arrival <- if (!is.null(control$immigration) && control$immigration$rate > 0) stats::rexp(1, control$immigration$rate) else Inf
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

  restart_step <- function(k) if (is.null(control$dt0)) 0.01 * max(range_z) / max(max(abs(k$f)), 1e-12) else control$dt0

  while (is.null(outcome)) {
    if (t >= control$t_max) { outcome <- "t_max"; break }
    if (steps >= control$max_steps) { outcome <- "max_steps"; break }
    h <- min(h, control$t_max - t, next_arrival - t)
    attempt <- step(rhs$f, z, h, k1, control)
    if (!attempt$accept) { h <- attempt$h_next; next }
    steps <- steps + 1L
    t <- t + h
    z <- attempt$y
    k1 <- attempt$k_end
    h <- attempt$h_next

    ## an immigrant arrives: it establishes if it can invade (and, in the
    ## stochastic mode, survives the lottery of being rare); the equilibrium
    ## solve and the extinction rule below decide whom it displaces
    if (t >= next_arrival) {
      immigration_attempts <- immigration_attempts + 1L
      z_new <- immigrant_phenotype(control, tf, bounds_z, k)
      s_new <- as.numeric(community_fitness_function(k1$community)(tf$inv(matrix(z_new, 1, k))))
      establishes <- is.finite(s_new) && s_new > 0 &&
        (control$immigration$establish == "deterministic" || stats::runif(1) < min(1, 2 * s_new))
      if (establishes && length(lineage) < control$max_residents) {
        z <- as.numeric(rbind(matrix(z, length(lineage), k), z_new))
        lineage <- c(lineage, next_lineage)
        event(t, "immigrant", next_lineage)
        next_lineage <- next_lineage + 1L
        rhs$reset()
        k1 <- rhs$f(z)
        h <- restart_step(k1)
      }
      next_arrival <- t + stats::rexp(1, control$immigration$rate)
    }

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

    ## near stationarity, finish with Newton on the coalition's gradients. Near a
    ## stable coalition the dynamics are stiff: the explicit stepper's step is
    ## capped by the fastest mode's stability limit, and at that step it jitters
    ## about the fixed point without ever meeting gradient_tol (DD99's pair:
    ## |g| stalls at ~6e-4 for hundreds of time units). An L-stable stepper
    ## converges, but still has to follow the approach; Newton does not.
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
      if (control$branch == "none") { outcome <- "stable"; break }
      comm <- k1$community
      mm <- length(lineage)
      zm <- matrix(z, mm, k)
      xm <- tf$inv(zm)
      V <- if (is.null(control$vcv)) diag(1, k) else control$vcv
      weight <- if (control$density) k1$n / n_total0 else rep(1, mm)
      candidates <- list()
      for (i in seq_len(mm)) {
        H <- community_fitness_hessian(comm, xm[i, , drop = FALSE])
        if (identical(tf$scale, "log")) H <- H * outer(xm[i, ], xm[i, ])
        e <- eigen((H + t(H)) / 2, symmetric = TRUE)
        if (max(e$values) <= control$classify_tol) next
        v <- e$vectors[, which.max(e$values)]
        v <- v / sqrt(sum(v^2))
        if (control$branch == "immediate") {
          candidates[[length(candidates) + 1L]] <- list(i = i, v = v, wait = 0, distance = control$branch_distance * max(range_z), split = TRUE)
        } else {
          sd_v <- sqrt(sum((v * mutation_sd)^2))
          br <- canonical_branch_rate(base, comm, tf, control, zm, i, v, sd_v, weight[i])
          if (br$rate <= 0) { event(t, "unbranchable", lineage[i]); next }
          if (control$branch == "expected") {
            wait <- 1 / br$rate
            distance <- sum(abs(br$delta) * br$mass) / sum(br$mass)
            side <- if (sum(br$mass[br$delta > 0]) >= sum(br$mass[br$delta < 0])) 1 else -1
          } else {
            wait <- stats::rexp(1, br$rate)
            j <- sample.int(length(br$delta), 1L, prob = br$mass)
            distance <- abs(br$delta[j]); side <- sign(br$delta[j])
          }
          candidates[[length(candidates) + 1L]] <- list(i = i, v = v * side, wait = wait, distance = distance, split = FALSE)
        }
      }
      if (length(candidates) == 0L) { outcome <- "stable"; break }
      if (mm + 1L > control$max_residents && control$branch != "immediate") { outcome <- "max_residents"; break }
      if (control$branch == "immediate" && mm + length(candidates) > control$max_residents) { outcome <- "max_residents"; break }
      ## the resident whose branching mutation arrives first branches; the
      ## community sits at the stationary point meanwhile
      first <- candidates[[which.min(vapply(candidates, `[[`, numeric(1), "wait"))]]
      todo <- if (control$branch == "immediate") candidates else list(first)
      if (first$wait > 0) {
        if (t + first$wait > control$t_max) { t <- control$t_max; record(t, z, k1); outcome <- "t_max"; break }
        t <- t + first$wait
        record(t, z, k1)
      }
      new_z <- zm
      new_lineage <- lineage
      for (cand in todo) {
        d <- cand$v * cand$distance
        if (cand$split) {
          new_z[cand$i, ] <- zm[cand$i, ] - d
          new_z <- rbind(new_z, zm[cand$i, ] + d)
        } else {
          new_z <- rbind(new_z, zm[cand$i, ] + d)       # the parent stays; the mutant appears
        }
        new_lineage <- c(new_lineage, next_lineage)
        event(t, "branch", lineage[cand$i])
        next_lineage <- next_lineage + 1L
      }
      z <- as.numeric(new_z)
      lineage <- new_lineage
      rhs$reset()
      k1 <- rhs$f(z)
      record(t, z, k1)
      h <- restart_step(k1)
    }
  }
  event(t, outcome)

  structure(list(
    trajectory = tibble::as_tibble(do.call(rbind, rows)),
    events = tibble::as_tibble(do.call(rbind, events)),
    community = rhs$last(),
    base = base,
    evaluations = rhs$evaluations(),
    steps = steps,
    immigration_attempts = immigration_attempts,
    outcome = outcome,
    trait_names = trait_names,
    trait_scale = tf$scale,
    mutation_sd = mutation_sd,
    control = control
  ), class = "canonical_equation")
}

##' @export
print.canonical_equation <- function(x, ...) {
  n_branch <- sum(x$events$event == "branch")
  final <- x$trajectory[x$trajectory$time == max(x$trajectory$time), ]
  n_imm <- sum(x$events$event == "immigrant")
  cat(sprintf("<canonical_equation: %s after %d steps (%d equilibrium solves), t = %.3g; %d branching event%s%s, %d resident%s at the end>\n",
              x$outcome, x$steps, x$evaluations, max(x$trajectory$time),
              n_branch, if (n_branch == 1L) "" else "s",
              if (x$immigration_attempts > 0L) sprintf(", %d of %d immigrants established", n_imm, x$immigration_attempts) else "",
              nrow(final), if (nrow(final) == 1L) "" else "s"))
  invisible(x)
}

##' The community at a recorded time of a canonical-equation run.
##'
##' Rebuilds the community from the record nearest \code{time} --- the
##' residents at their traits then --- and solves it to equilibrium, so the
##' fitness landscape, invasibility surface or classification at that moment
##' can be computed with the usual functions.
##'
##' @title Community along a canonical-equation trajectory
##' @param x A \code{canonical_equation}.
##' @param time Evolutionary time; the nearest recorded time is used.
##' @return A solved \code{community}, with \code{attr(., "time")} the record
##' used.
##' @author Daniel Falster
##' @export
canonical_community <- function(x, time) {
  tr <- x$trajectory
  times <- unique(tr$time)
  at <- times[which.min(abs(times - time))]
  rows <- tr[tr$time == at, , drop = FALSE]
  traits <- as.matrix(rows[, x$trait_names, drop = FALSE])
  comm <- x$base |>
    community_add(trait_matrix(traits, x$trait_names), birth_rate = rows$density) |>
    community_demography()
  attr(comm, "time") <- at
  comm
}

## Times worth looking at by default: the record just before and after each
## branching, or failing any, the first and last.
canonical_landscape_times <- function(x) {
  br <- x$events$time[x$events$event == "branch"]
  times <- unique(x$trajectory$time)
  if (length(br) == 0L) return(range(times))
  out <- numeric(0)
  for (b in unique(br)) {
    before <- times[times < b]
    after <- times[times > b]
    out <- c(out, if (length(before)) max(before) else b, b, if (length(after)) min(after) else b)
  }
  unique(out)
}

##' Plot trait trajectories from the canonical equation.
##'
##' \code{type = "trajectory"}: each lineage is a line of trait value against
##' evolutionary time; a branching event is where one line becomes two, an
##' immigrant a line that starts in mid-air. With several traits, one panel per
##' trait. \code{type = "landscapes"}: the invasion-fitness landscape of the
##' community at chosen \code{times} (by default the record just before, at
##' and after each branching), residents marked, one panel per time --- the
##' fitness minimum forming and the daughters moving apart. Single-trait
##' communities only.
##'
##' @title Plot a canonical-equation run
##' @param x A \code{canonical_equation}.
##' @param type \code{"trajectory"} or \code{"landscapes"}.
##' @param times For \code{"landscapes"}: evolutionary times to show.
##' @param n_evals For \code{"landscapes"}: mutant values per landscape.
##' @param object A \code{canonical_equation}, for \code{autoplot}.
##' @param ... Ignored.
##' @return A \code{ggplot} object.
##' @author Daniel Falster
##' @export
plot.canonical_equation <- function(x, type = c("trajectory", "landscapes"), times = NULL,
                                    n_evals = 101L, ...) {
  type <- match.arg(type)
  if (type == "landscapes") return(canonical_plot_landscapes(x, times, n_evals))
  tr <- x$trajectory
  long <- do.call(rbind, lapply(x$trait_names, function(nm) {
    data.frame(time = tr$time, lineage = tr$lineage, trait = nm, value = tr[[nm]])
  }))
  marks <- x$events[x$events$event %in% c("branch", "immigrant"), ]
  p <- ggplot2::ggplot(long, ggplot2::aes(x = .data[["time"]], y = .data[["value"]],
                                          group = .data[["lineage"]])) +
    ggplot2::geom_line(linewidth = 0.6)
  if (nrow(marks) > 0L) {
    p <- p + ggplot2::geom_vline(xintercept = marks$time, linetype = "dotted", colour = "grey50")
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

canonical_plot_landscapes <- function(x, times, n_evals) {
  if (length(x$trait_names) != 1L) stop("Landscape panels need a single-trait run")
  if (is.null(times)) times <- canonical_landscape_times(x)
  tf <- pip_transform(x$trait_scale)
  b <- as.numeric(x$base$bounds[1, ])
  grid <- tf$inv(seq(tf$fwd(b[1]), tf$fwd(b[2]), length.out = n_evals))
  panels <- lapply(times, function(tt) {
    comm <- canonical_community(x, tt)
    at <- attr(comm, "time")
    s <- as.numeric(comm$fitness_function(grid))
    res <- as.numeric(comm$traits[, 1])
    list(curve = data.frame(time = at, x = grid, fitness = s),
         residents = data.frame(time = at, x = res, fitness = as.numeric(comm$fitness_function(res))))
  })
  curve <- do.call(rbind, lapply(panels, `[[`, "curve"))
  residents <- do.call(rbind, lapply(panels, `[[`, "residents"))
  shown <- sort(unique(curve$time))
  labels <- make.unique(sprintf("t = %.4g", shown), sep = " #")
  curve$panel <- factor(labels[match(curve$time, shown)], levels = labels)
  residents$panel <- factor(labels[match(residents$time, shown)], levels = labels)
  p <- ggplot2::ggplot(curve, ggplot2::aes(x = .data[["x"]], y = .data[["fitness"]])) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", colour = "grey60") +
    ggplot2::geom_line(linewidth = 0.6) +
    ggplot2::geom_point(data = residents, colour = "firebrick", size = 2.5) +
    ggplot2::facet_wrap(~ panel) +
    ggplot2::xlab(x$trait_names) + ggplot2::ylab("invasion fitness") +
    ggplot2::theme_classic() +
    ggplot2::theme(text = ggplot2::element_text(size = 14))
  if (identical(x$trait_scale, "log")) p <- p + ggplot2::scale_x_log10()
  p
}

##' @rdname plot.canonical_equation
##' @export
autoplot.canonical_equation <- function(object, ...) plot.canonical_equation(object, ...)
