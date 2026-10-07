# Two-worker future plans for the tests of the parallel map.
#
# "multicore" forks this session, so its workers run whatever code is loaded
# here. "multisession" starts fresh R processes, which would load regnans from
# the library -- under devtools::load_all() not the code under test -- so each
# worker loads this source tree instead (a cluster of two such processes is
# what plan(multisession) builds). Tests assert the workers are separate
# processes running this source, rather than trusting the plan.
local_two_workers <- function(type = c("multicore", "multisession"), env = parent.frame()) {
  type <- match.arg(type)
  testthat::skip_if_not_installed("future")
  testthat::skip_if_not_installed("future.apply")
  testthat::skip_if_not_installed("pkgload")
  if (type == "multicore") {
    testthat::skip_on_os("windows")
    testthat::skip_if_not(future::supportsMulticore())
    old <- future::plan(future::multicore, workers = 2)
  } else {
    cl <- parallel::makePSOCKcluster(2)
    withr::defer(parallel::stopCluster(cl), envir = env)
    if (pkgload::is_dev_package("regnans")) {
      parallel::clusterCall(cl, function(path) {
        suppressMessages(pkgload::load_all(path, compile = FALSE, helpers = FALSE, quiet = TRUE))
        NULL
      }, getNamespaceInfo("regnans", "path"))
    }
    old <- future::plan(future::cluster, workers = cl)
  }
  withr::defer(future::plan(old), envir = env)
  invisible(type)
}

# Evaluate `code` under a two-worker plan of the given type, asserting first
# that the map sees both workers (a test that fell back to lapply() would
# otherwise pass wherever the answer does not depend on the plan).
with_two_workers <- function(type, code) {
  local_two_workers(type)
  testthat::expect_identical(regnans_workers(), 2L)
  code
}

# Which process ran each of two jobs, and from which source.
worker_identity <- function() {
  regnans_map(1:2, function(i) {
    Sys.sleep(0.2)
    list(pid = Sys.getpid(), path = getNamespaceInfo("regnans", "path"))
  })
}
