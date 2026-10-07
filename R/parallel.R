# Independent model evaluations -- the points of a finite-difference stencil,
# residents on a grid, points in a scouting design -- go through one map so
# that a user's future::plan() parallelises them without any function here
# knowing how. With no plan it is lapply().
#
# What crosses to a worker is a community without a solve (a base community,
# traits, seed densities, the equilibrium solver's state); the worker solves it
# and returns numbers. A solved community does not cross: on plant its fitness
# function closes over the SCM, an external pointer that only a forked worker
# can use. So the functions run on workers are defined at top level and take
# their data as arguments (a closure would carry its enclosing frame, and any
# solved community in it, along unseen), and under a plan whose workers are
# not forks the map refuses work that holds a pointer.
#
# The work is deterministic, so workers are given no random seeds and a plan
# never moves the caller's random-number stream.

## The number of workers the active future plan provides (1 without a plan,
## or when parallel is FALSE).
regnans_workers <- function(parallel = TRUE) {
  if (!parallel || !requireNamespace("future", quietly = TRUE)) {
    return(1L)
  }
  n <- as.integer(future::nbrOfWorkers())
  if (n > 1L && !requireNamespace("future.apply", quietly = TRUE)) {
    stop("A future plan with ", n, " workers is set, but regnans maps over it with ",
         "future.apply, which is not installed; install it, or set future::plan(future::sequential)")
  }
  n
}

regnans_map <- function(X, FUN, ..., parallel = TRUE) {
  if (regnans_workers(parallel) == 1L) {
    return(lapply(X, FUN, ...))
  }
  if (!regnans_plan_forks() && regnans_pointers(list(X, FUN, ...)) > 0L) {
    stop("This work holds an external pointer (on plant, the SCM inside a solved ",
         "community's fitness function), which the workers of the current future plan, ",
         "separate R processes, cannot use. Use future::plan(future::multicore), whose ",
         "workers are forks of this session, or future::plan(future::sequential)")
  }
  future.apply::future_lapply(X, FUN, ..., future.seed = NULL)
}

## Are the active plan's workers forks of this process, sharing its memory?
regnans_plan_forks <- function() {
  inherits(future::plan("next"), "multicore")
}

## The number of external pointers reachable from x, which no serialisation
## can carry to another process.
regnans_pointers <- function(x) {
  n <- 0L
  serialize(x, NULL, refhook = function(o) {
    if (typeof(o) == "externalptr") n <<- n + 1L
    NULL
  })
  n
}

## Split n jobs into one contiguous chunk per worker, so that work which
## benefits from its neighbour's result (warm-started equilibria along a trait
## grid) keeps that benefit inside each chunk.
regnans_chunks <- function(n, parallel = TRUE) {
  if (n == 0L) return(list())
  k <- min(regnans_workers(parallel), n)
  unname(split(seq_len(n), ceiling(seq_len(n) * k / n)))
}
