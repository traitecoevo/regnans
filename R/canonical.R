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
# chosen for how few evaluations it needs, not for accuracy in the trajectory.
# The stepping is odelia's: its OdeSolver takes the right-hand side as an R
# closure, steps once at a time up to a bound (so branching and immigration
# clocks are landed on exactly), re-seeds the state at a new length when the
# community changes, and offers explicit embedded pairs or RODAS for the stiff
# case of many traits or species (#43), whose Jacobian of the right-hand side
# is by finite differences over the residents today and exact once the model
# supplies equilibrium sensitivities. RODAS is the default because trait
# models are stiff as a rule -- selection on different traits and species runs
# on widely separated time scales, and near a coalition an explicit pair is
# held at its stability limit -- while on smooth stretches it costs 1.5-2x the
# explicit pairs in equilibrium solves (scripts/canonical-stepper-benchmark.R),
# a price an exact Jacobian removes.

##' Control the canonical equation.
##'
##' @title Canonical-equation settings
##' @param control A list of values to modify from the defaults:
##' \describe{
##'   \item{\code{rate}}{the mutation rate \eqn{\mu}: mutations per unit
##'     evolutionary time in the initial community (the resident's share scales
##'     with its relative density when \code{density} is on). It sets the time
##'     scale of the speed and of the waiting times alike.}
##'   \item{\code{vcv}}{mutational covariance on the trait scale, a k by k
##'     matrix. \code{NULL} uses \code{diag(mutation_sd^2)}, so the speed and
##'     the branching waiting time are driven by the same kernel.}
##'   \item{\code{density}}{multiply each resident's speed by its equilibrium
##'     density relative to the initial total (mutation supply scales with
##'     population size).}
##'   \item{\code{stepper}}{an \code{odelia} stepper: \code{"rodas"} (linearly
##'     implicit RODAS4(3), the default: trait models are stiff as a rule, and
##'     an L-stable step is not held at a stability limit; its Jacobian is
##'     formed by finite differences across the residents, one equilibrium
##'     solve per unknown), \code{"rkck"} (explicit Cash--Karp 4(5), the
##'     cheapest on smooth stretches) or \code{"dopri"} (explicit
##'     Dormand--Prince 5(4)). See \code{\link[odelia]{OdeSolver}}.}
##'   \item{\code{t_max}, \code{max_steps}}{limits on evolutionary time
##'     (\code{Inf} by default: a run ends when the coalition is stable or a
##'     limit is reached) and accepted steps.}
##'   \item{\code{rtol}, \code{atol}}{tolerances of the step-size control, on
##'     the trait scale: a step is accepted when its error estimate is within
##'     \code{atol + rtol * abs(z)}.}
##'   \item{\code{dt0}, \code{dt_max}}{initial and largest step; \code{NULL}
##'     chooses \code{dt0} from the initial speed.}
##'   \item{\code{gradient_tol}}{selection has stopped when every resident's
##'     gradient on the trait scale is below this.}
##'   \item{\code{polish}, \code{polish_tol}}{once every gradient is below
##'     \code{polish_tol}, finish the approach to the stationary coalition by
##'     solving for it (\code{\link{community_solve_singularity}}) instead of
##'     integrating it: near a stable coalition the dynamics are
##'     stiff, and an explicit stepper held at its stability limit jitters
##'     about the point without reaching \code{gradient_tol}. A solution more
##'     than a tenth of the trait range from the current residents is not a
##'     local finish and is refused; the integration carries on.}
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
##'   \item{\code{establishment_factor}}{a rare mutant or immigrant with
##'     invasion fitness \eqn{s > 0} establishes with probability
##'     \eqn{\min(1, c\,s)}. The default \eqn{c = 2} is Haldane's
##'     approximation for a discrete-generation Poisson branching process;
##'     fitness is a per-capita rate for DD99, a log ratio for the other
##'     reference models and a log net-reproduction ratio for plant, so the
##'     right constant is model-dependent and the waiting times are exact only up
##'     to it.}
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
    stepper = "rodas",
    t_max = Inf, max_steps = 2000L,
    rtol = 1e-3, atol = 1e-5, dt0 = NULL, dt_max = Inf,
    gradient_tol = 1e-4,
    polish = TRUE, polish_tol = 1e-2,
    branch = "expected", mutation_sd = NULL, branch_nodes = 12L, branch_distance = 0.02,
    establishment_factor = 2,
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
  ret$stepper <- match.arg(ret$stepper, c("rodas", "rkck", "dopri"))
  for (nm in c("rate", "rtol", "atol", "gradient_tol", "polish_tol", "branch_distance",
               "extinct_fraction", "classify_tol", "establishment_factor")) {
    v <- ret[[nm]]
    if (!is.numeric(v) || length(v) != 1L || !is.finite(v) || v <= 0) {
      stop(nm, " must be a single positive number")
    }
  }
  if (!is.numeric(ret$t_max) || length(ret$t_max) != 1L || is.na(ret$t_max) || ret$t_max <= 0) stop("t_max must be positive")
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
## starts and the whole last evaluation (speeds, gradients, densities,
## community) kept for the caller: odelia guarantees that after an accepted
## step or a re-seed the last evaluation was at the solver's state. `ode` is
## the same map in the (t, y) form the solver calls; an evaluation whose
## equilibrium did not converge or whose speed is not finite is refused with
## odelia::domain_error(), which rejects the step and retries it smaller.
canonical_rhs <- function(base, tf, control, k, n_total0, V) {
  birth_rate <- NULL
  state <- NULL
  last <- NULL
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
    g <- community_fitness_gradient(comm)                       # m x k, raw units
    g_z <- if (identical(tf$scale, "log")) g * x else g           # w.r.t. the trait scale
    weight <- if (control$density) n / n_total0 else rep(1, m)
    speed <- 0.5 * control$rate * weight * (g_z %*% V)
    last <<- list(f = as.numeric(speed), g_z = g_z, n = n, community = comm)
  }
  ode <- function(t, y) {
    out <- f(y)
    if (!isTRUE(attr(out$community, "converged")) || !all(is.finite(out$f))) {
      odelia::domain_error("the demographic equilibrium did not converge")
    }
    out$f
  }
  list(f = f, ode = ode, last = function() last, evaluations = function() evaluations,
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
## mutational kernel of width sd_v: one vectorised fitness call gives every
## node's invasion fitness and so its establishment probability, and for the
## nodes that can invade, one equilibrium solve with the mutant in the
## resident's place (warm-started from the residents' densities) says whether
## the resident can invade it back -- a protected dimorphism. A solve that fails
## is counted and warned about, not silently taken as "no coexistence".
## Returns the rate and the nodes with their success weights, so a distance can
## be taken (expected or drawn).
canonical_branch_rate <- function(base, comm, tf, control, zm, i, v, sd_v, weight_i) {
  k <- ncol(zm)
  gh <- gauss_hermite(control$branch_nodes)
  delta <- sqrt(2) * sd_v * gh$x
  wq <- gh$w / sqrt(pi)
  z_mut <- sweep(outer(delta, v), 2, zm[i, ], "+")                 # nodes x k
  s_mut <- as.numeric(community_fitness_function(comm)(tf$inv(z_mut)))
  p_est <- ifelse(is.finite(s_mut) & s_mut > 0 & abs(delta) > 1e-12,
                  pmin(1, control$establishment_factor * s_mut), 0)
  success <- numeric(length(delta))
  failed <- 0L
  n_res <- as.numeric(comm$birth_rate)
  for (j in which(p_est > 0)) {
    z_other <- zm; z_other[i, ] <- z_mut[j, ]
    other <- tryCatch(
      base |>
        community_add(trait_matrix(tf$inv(z_other), base$trait_names), birth_rate = n_res) |>
        community_demography(),
      error = function(e) NULL)
    if (is.null(other) || !isTRUE(attr(other, "converged"))) { failed <- failed + 1L; next }
    s_back <- as.numeric(community_fitness_function(other)(tf$inv(matrix(zm[i, ], 1, k))))
    if (is.finite(s_back) && s_back > 0) success[j] <- p_est[j]
  }
  if (failed > 0L) {
    warning(sprintf("%d of %d coexistence solves failed while computing a branching rate; the rate is a lower bound",
                    failed, sum(p_est > 0)))
  }
  mass <- wq * success
  list(rate = control$rate * weight_i * sum(mass), delta = delta, mass = mass, failed = failed)
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

## Finish the approach to a stationary coalition by solving for it: the
## residents, at densities n, handed to community_solve_singularity(), which
## roots every resident's gradient with the equilibrium re-solved inside and
## the Jacobian across the residents from the derivative layer. Returns the
## new z and the equilibrium solves it cost, with z NULL if the solve did not
## converge or lost the coalition, in which case the integration simply
## carries on. Polishing is a local finish: a small gradient does not mean the
## root is near (a flat landscape, residents close to merging), and the
## solver's full steps can carry it to a different singular coalition, so a
## root more than `reach` of the trait range from z on any trait is refused.
canonical_polish <- function(base, z, n, tf, k, range_z, control, reach = 0.1) {
  sol <- tryCatch(
    suppressWarnings(community_solve_singularity(
      base, x0 = tf$inv(matrix(z, ncol = k)), tol = control$gradient_tol / 10, maxit = 50, birth_rate = n)),
    error = function(e) e)
  if (inherits(sol, "error")) {
    return(list(z = NULL, evaluations = if (is.null(sol$evaluations)) 0L else sol$evaluations))
  }
  evaluations <- attr(sol, "evaluations")
  z_new <- as.numeric(tf$fwd(sol$traits))
  m <- length(z) / k
  if (!isTRUE(attr(sol, "converged")) ||
      any(abs(z_new - z) > rep(reach * range_z, each = m))) {
    return(list(z = NULL, evaluations = evaluations))
  }
  list(z = z_new, evaluations = evaluations)
}

## odelia's step-size control in the canonical equation's terms: the error
## level is atol + rtol * |z| (odelia's a_y = 1, a_dydt = 0 defaults) and the
## largest step dt_max; the step the solver tries first is set by the caller
## from the initial speed.
canonical_ode_control <- function(control) {
  ctl <- odelia::OdeControl$new()
  ctl$set_tol_rel(control$rtol)
  ctl$set_tol_abs(control$atol)
  ctl$set_step_size_max(control$dt_max)
  ctl
}

##' Integrate the canonical equation of adaptive dynamics.
##'
##' Every resident's trait moves up its selection gradient,
##' \deqn{\frac{dz_i}{dt} = \tfrac{1}{2}\,\mu\; \tilde n_i\; V\, \nabla_z s(z_i)}
##' (Dieckmann & Law 1996), on the community's trait scale \eqn{z}, with the
##' community re-solved to demographic equilibrium at every evaluation;
##' \eqn{\mu} is \code{rate}, \eqn{\tilde n_i} the resident's equilibrium
##' density relative to the initial total and \eqn{V} the mutational
##' covariance (see \code{\link{canonical_control}}). A resident whose own
##' gradient has vanished is tested: a positive eigenvalue of its fitness
##' Hessian means it sits at a fitness minimum. By default it then branches
##' after the waiting time for a mutation that can invade and coexist, the
##' parent staying put and the mutant appearing at the kernel's successful
##' distance along that eigenvector; branching and immigration are competing
##' clocks, and whichever comes first happens. When every resident is at a
##' fitness maximum the coalition is evolutionarily stable and the run stops
##' (unless immigration is on, in which case it runs to \code{t_max}). A
##' resident whose density falls below \code{extinct_fraction} of the total is
##' removed.
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
##' \code{"branch"}, \code{"immigrant"}, \code{"immigrant_refused"} (it
##' could have established but the community is at \code{max_residents}),
##' \code{"extinct"}, \code{"polish"}, \code{"stable"}, \code{"t_max"},
##' \code{"max_steps"}, \code{"max_residents"} --- and \code{lineage}),
##' \code{community} (the final solved community), \code{evaluations} (how
##' many equilibrium solves it cost), \code{steps} and \code{rejections}
##' (accepted steps and attempts the step-size control rejected),
##' \code{immigration_attempts},
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
  x0 <- if (is.matrix(x0)) x0 else matrix(x0, ncol = k)
  if (ncol(x0) != k) stop("x0 must have ", k, " column(s)")
  m <- nrow(x0)
  base <- community_clear_residents(community)
  V <- if (is.null(control$vcv)) diag(mutation_sd^2, k) else control$vcv
  immigrating <- !is.null(control$immigration) && control$immigration$rate > 0

  ## initial solve sets the density scale and the first right-hand side
  first <- base |> community_add(trait_matrix(x0, trait_names)) |> community_demography()
  n_total0 <- sum(as.numeric(first$birth_rate))
  if (!is.finite(n_total0) || n_total0 <= 0) stop("The starting community has no positive equilibrium density")
  rhs <- canonical_rhs(base, tf, control, k, n_total0, V)

  z <- as.numeric(tf$fwd(x0))
  lineage <- seq_len(m)
  next_lineage <- m + 1L
  t <- 0
  solver <- odelia::OdeSolver$new(rhs$ode, z, t0 = t, method = control$stepper, autonomous = TRUE,
                                  control = canonical_ode_control(control))
  k1 <- rhs$last()
  restart_step <- function(k) if (is.null(control$dt0)) min(control$dt_max, 0.01 * max(range_z) / max(max(abs(k$f)), 1e-12)) else control$dt0
  solver$set_step_size(restart_step(k1))

  rows <- list()
  events <- list()
  immigration_attempts <- 0L
  next_arrival <- if (immigrating) stats::rexp(1, control$immigration$rate) else Inf
  pending <- NULL        # a scheduled branching: list(time, i, v, distance)

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
  ## the solver is re-seeded at the current z (and t), which evaluates the
  ## right-hand side there once
  reseed <- function() {
    solver$set_state(z, t)
    k1 <<- rhs$last()
  }
  ## the community changed: new densities and Jacobian, new first right-hand
  ## side, a fresh step, and any scheduled branching is void
  changed <- function() {
    rhs$reset()
    reseed()
    solver$set_step_size(restart_step(k1))
    pending <<- NULL
  }
  resident_stationary <- function(k) apply(abs(k$g_z), 1, max) < control$gradient_tol

  record(t, z, k1)
  outcome <- NULL
  steps <- 0L
  polish_evaluations <- 0L

  while (is.null(outcome)) {
    if (t >= control$t_max) { outcome <- "t_max"; break }
    if (steps >= control$max_steps) { outcome <- "max_steps"; break }

    ## due events happen before any step is taken
    if (!is.null(pending) && t >= pending$time) {
      mm <- length(lineage)
      zm <- matrix(z, mm, k)
      z <- as.numeric(rbind(zm, zm[pending$i, ] + pending$v * pending$distance))  # the parent stays; the mutant appears
      event(t, "branch", lineage[pending$i])
      lineage <- c(lineage, next_lineage)
      next_lineage <- next_lineage + 1L
      changed()
      record(t, z, k1)
      next
    }
    if (t >= next_arrival) {
      immigration_attempts <- immigration_attempts + 1L
      z_new <- immigrant_phenotype(control, tf, bounds_z, k)
      s_new <- as.numeric(community_fitness_function(k1$community)(tf$inv(matrix(z_new, 1, k))))
      can_invade <- is.finite(s_new) && s_new > 0 &&
        (control$immigration$establish == "deterministic" ||
           stats::runif(1) < min(1, control$establishment_factor * s_new))
      if (can_invade && length(lineage) >= control$max_residents) {
        event(t, "immigrant_refused", NA_integer_)
      } else if (can_invade) {
        z <- as.numeric(rbind(matrix(z, length(lineage), k), z_new))
        event(t, "immigrant", next_lineage)
        lineage <- c(lineage, next_lineage)
        next_lineage <- next_lineage + 1L
        changed()
        record(t, z, k1)
      }
      next_arrival <- t + stats::rexp(1, control$immigration$rate)
      next
    }

    ## one error-controlled step, stopping at whichever comes first: a
    ## scheduled branching, an arrival, or the end
    bound <- min(control$t_max, next_arrival, if (is.null(pending)) Inf else pending$time)
    if (bound <= t) stop("internal error: the step bound is not ahead")   # the clocks above guarantee this
    solver$step(bound)
    steps <- steps + 1L
    t <- solver$time()
    z <- solver$state()
    k1 <- rhs$last()

    ## extinction: a resident the others have driven out stops counting
    mm <- length(lineage)
    frac <- k1$n / sum(k1$n)
    gone <- which(frac < control$extinct_fraction)
    if (length(gone) > 0L && length(gone) < mm) {
      for (i in gone) event(t, "extinct", lineage[i])
      keep <- setdiff(seq_len(mm), gone)
      z <- as.numeric(matrix(z, mm, k)[keep, , drop = FALSE])
      lineage <- lineage[keep]
      changed()
    }
    record(t, z, k1)

    ## near stationarity of the whole coalition, finish with Newton on its
    ## gradients: the dynamics there are stiff -- the explicit stepper's step is
    ## capped by the fastest mode's stability limit and jitters about the point
    ## without ever meeting gradient_tol -- and even an L-stable stepper must
    ## follow the approach, which Newton skips
    if (control$polish && max(abs(k1$g_z)) < control$polish_tol &&
        max(abs(k1$g_z)) >= control$gradient_tol) {
      polished <- canonical_polish(base, z, k1$n, tf, k, range_z, control)
      polish_evaluations <- polish_evaluations + polished$evaluations
      if (!is.null(polished$z)) {
        z <- polished$z
        reseed()
        event(t, "polish")
        record(t, z, k1)
      }
    }

    ## a resident whose own gradient has vanished is tested for a fitness
    ## minimum; branching is scheduled as a clock competing with immigration
    stationary <- resident_stationary(k1)
    if (control$branch == "none") {
      if (all(stationary) && !immigrating) { outcome <- "stable"; break }
      next
    }
    if (is.null(pending) && any(stationary)) {
      comm <- k1$community
      mm <- length(lineage)
      zm <- matrix(z, mm, k)
      xm <- tf$inv(zm)
      weight <- if (control$density) k1$n / n_total0 else rep(1, mm)
      candidates <- list()
      for (i in which(stationary)) {
        H <- community_fitness_hessian(comm, xm[i, , drop = FALSE])
        if (identical(tf$scale, "log")) H <- H * outer(xm[i, ], xm[i, ])
        e <- eigen((H + t(H)) / 2, symmetric = TRUE)
        if (max(e$values) <= control$classify_tol) next
        v <- e$vectors[, which.max(e$values)]
        v <- v / sqrt(sum(v^2))
        if (control$branch == "immediate") {
          candidates[[length(candidates) + 1L]] <- list(i = i, v = v, rate = Inf,
                                                        distance = control$branch_distance * max(range_z))
        } else {
          sd_v <- sqrt(sum((v * mutation_sd)^2))
          br <- canonical_branch_rate(base, comm, tf, control, zm, i, v, sd_v, weight[i])
          if (br$rate <= 0) next
          candidates[[length(candidates) + 1L]] <- list(i = i, v = v, rate = br$rate,
                                                        delta = br$delta, mass = br$mass)
        }
      }
      if (length(candidates) == 0L) {
        if (all(stationary) && !immigrating) { outcome <- "stable"; break }
      } else if (mm >= control$max_residents) {
        if (!immigrating) { outcome <- "max_residents"; break }
      } else if (control$branch == "immediate") {
        ## every candidate splits now, daughters either side of the parent
        new_z <- zm
        new_lineage <- lineage
        for (cand in candidates) {
          d <- cand$v * cand$distance
          new_z[cand$i, ] <- zm[cand$i, ] - d
          new_z <- rbind(new_z, zm[cand$i, ] + d)
          new_lineage <- c(new_lineage, next_lineage)
          event(t, "branch", lineage[cand$i])
          next_lineage <- next_lineage + 1L
        }
        z <- as.numeric(new_z)
        lineage <- new_lineage
        changed()
        record(t, z, k1)
      } else {
        ## independent clocks: the first fires after 1/sum(rates) in
        ## expectation, and belongs to resident i with probability rate_i/sum
        rates <- vapply(candidates, `[[`, numeric(1), "rate")
        total <- sum(rates)
        if (control$branch == "expected") {
          wait <- 1 / total
          cand <- candidates[[which.max(rates)]]
          distance <- sum(abs(cand$delta) * cand$mass) / sum(cand$mass)
          side <- if (sum(cand$mass[cand$delta > 0]) >= sum(cand$mass[cand$delta < 0])) 1 else -1
        } else {
          wait <- stats::rexp(1, total)
          cand <- candidates[[sample.int(length(candidates), 1L, prob = rates)]]
          j <- sample.int(length(cand$delta), 1L, prob = cand$mass)
          distance <- abs(cand$delta[j]); side <- sign(cand$delta[j])
        }
        pending <- list(time = t + wait, i = cand$i, v = cand$v * side, distance = distance)
      }
    }
  }
  if (!is.null(outcome) && outcome %in% c("t_max", "max_steps") && t != rows[[length(rows)]]$time[1]) record(t, z, k1)
  event(t, outcome)

  structure(list(
    trajectory = tibble::as_tibble(do.call(rbind, rows)),
    events = tibble::as_tibble(do.call(rbind, events)),
    community = rhs$last(),
    base = base,
    evaluations = rhs$evaluations() + polish_evaluations,
    steps = steps,
    rejections = unname(solver$counts()[["n_rejections"]]),
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
  n_ref <- sum(x$events$event == "immigrant_refused")
  cat(sprintf("<canonical_equation: %s after %d steps (%d equilibrium solves), t = %.3g; %d branching event%s%s, %d resident%s at the end>\n",
              x$outcome, x$steps, x$evaluations, max(x$trajectory$time),
              n_branch, if (n_branch == 1L) "" else "s",
              if (x$immigration_attempts > 0L) sprintf(", %d of %d immigrants established%s", n_imm, x$immigration_attempts,
                                                        if (n_ref > 0L) sprintf(" (%d refused at max_residents)", n_ref) else "") else "",
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
