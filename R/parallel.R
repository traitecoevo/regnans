# Independent model evaluations -- residents on a grid, points in a scouting
# design -- go through one map so that a user's future::plan() parallelises
# them without any function here knowing how. With no plan (or the packages
# absent) it is lapply().

## Whether a future plan with more than one worker is in effect.
regnans_workers <- function(parallel = TRUE) {
  if (!parallel || !requireNamespace("future", quietly = TRUE) ||
      !requireNamespace("future.apply", quietly = TRUE)) {
    return(1L)
  }
  as.integer(future::nbrOfWorkers())
}

regnans_map <- function(X, FUN, ..., parallel = TRUE) {
  if (regnans_workers(parallel) > 1L) {
    future.apply::future_lapply(X, FUN, ..., future.seed = TRUE)
  } else {
    lapply(X, FUN, ...)
  }
}

## Split n jobs into one contiguous chunk per worker, so that work which
## benefits from its neighbour's result (warm-started equilibria along a trait
## grid) keeps that benefit inside each chunk.
regnans_chunks <- function(n, parallel = TRUE) {
  if (n == 0L) return(list())
  k <- min(regnans_workers(parallel), n)
  unname(split(seq_len(n), ceiling(seq_len(n) * k / n)))
}
