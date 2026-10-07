# Closed-form oracles for singular coalitions, computed independently of the
# package's machinery (no harness, no derivative layer).
#
# DD99 is Lotka-Volterra with Gaussian kernels: for residents x the equilibrium
# solves A n = K with A_ij = C(x_i - x_j), and invasion fitness is
#   s(y) = r (1 - sum_j C(y - x_j) n_j / K(y)),
# so with sum_j C(x_i - x_j) n_j = K(x_i) at equilibrium, resident i's
# selection gradient is
#   g_i = r [sum_j (x_i - x_j) C(x_i - x_j) n_j / (sigma_C^2 K(x_i)) - (x_i - x0)/sigma_K^2].
# For a mirror-image pair +-a the condition g = 0 reduces to
#   C(2a) / (1 + C(2a)) = sigma_C^2 / (2 sigma_K^2),
# so a* = sigma_C sqrt(-log(c / (1 - c)) / 2), c = sigma_C^2 / (2 sigma_K^2),
# which exists iff sigma_C < sigma_K. A symmetric triple (-b, 0, b) has a
# middle resident with no gradient by symmetry and b* a scalar root.
#
# GK98 (Levene soft selection, patches -d, 0, d, equal capacities): a
# mirror-image pair +-a holds equal shares, and resident a's gradient is
# proportional to sum_j (mu_j - a) f_j(a) / D_j with D_j = f_j(a) + f_j(-a).
#
# GM99 (seed size, safe-site lottery) has no symmetry, so its dimorphism is the
# asymmetric oracle: residents at different distances from anything, at
# densities ~3.4 : 1. Coded here from the model's definition,
#   W(y) = (R/y) s(y) E[c(y) / (c(y) + Z)],  Z = sum_j k_j c(x_j),  k_j ~ Poisson(N_j),
#   s(y) = 1 - 2 exp(-beta y),  c(y) = exp(alpha y),
# with the densities solving W_i = 1 for every resident and the gradient
#   d log W / dy = -1/y + s'(y)/s(y) + alpha E[c Z / (c + Z)^2] / E[c / (c + Z)].
#
# Lotka-Volterra (below the GM99 functions) is the oracle for protected
# coexistence: invasion fitness when rare is closed form, and an infeasible
# equilibrium (negative densities) is easy to make.

dd99_kernels <- function(sigma_C, sigma_K, x0, K0) {
  list(C = function(d) exp(-d^2 / (2 * sigma_C^2)),
       K = function(x) K0 * exp(-(x - x0)^2 / (2 * sigma_K^2)))
}

dd99_coalition_density <- function(x, sigma_C = 0.4, sigma_K = 1, x0 = 0, K0 = 500) {
  kn <- dd99_kernels(sigma_C, sigma_K, x0, K0)
  solve(kn$C(outer(x, x, "-")), kn$K(x))
}

dd99_coalition_gradient <- function(x, sigma_C = 0.4, sigma_K = 1, x0 = 0, r = 1, K0 = 500) {
  kn <- dd99_kernels(sigma_C, sigma_K, x0, K0)
  n <- dd99_coalition_density(x, sigma_C, sigma_K, x0, K0)
  D <- outer(x, x, "-")
  r * (as.numeric((D * kn$C(D)) %*% n) / (sigma_C^2 * kn$K(x)) - (x - x0) / sigma_K^2)
}

## d2s/dy2 at each resident: s = r (1 - A/K) with A(y) = sum_j C(y - x_j) n_j.
dd99_coalition_curvature <- function(x, sigma_C = 0.4, sigma_K = 1, x0 = 0, r = 1, K0 = 500) {
  kn <- dd99_kernels(sigma_C, sigma_K, x0, K0)
  n <- dd99_coalition_density(x, sigma_C, sigma_K, x0, K0)
  vapply(x, function(y) {
    d <- y - x
    A <- sum(kn$C(d) * n)
    A1 <- sum(-d / sigma_C^2 * kn$C(d) * n)
    A2 <- sum((d^2 / sigma_C^4 - 1 / sigma_C^2) * kn$C(d) * n)
    K <- kn$K(y)
    K1 <- -(y - x0) / sigma_K^2 * K
    K2 <- ((y - x0)^2 / sigma_K^4 - 1 / sigma_K^2) * K
    -r * (A2 / K - 2 * A1 * K1 / K^2 - A * K2 / K^2 + 2 * A * K1^2 / K^3)
  }, numeric(1))
}

dd99_pair_root <- function(sigma_C = 0.4, sigma_K = 1) {
  c0 <- sigma_C^2 / (2 * sigma_K^2)
  sigma_C * sqrt(-log(c0 / (1 - c0)) / 2)
}

dd99_triple_root <- function(sigma_C = 0.3, sigma_K = 1) {
  g_outer <- function(b) dd99_coalition_gradient(c(-b, 0, b), sigma_C, sigma_K)[3]
  stats::uniroot(g_outer, c(0.3, 1.2), tol = 1e-13)$root
}

gk98_pair_root <- function(d = 1.5, sigma = 1) {
  mu <- c(-d, 0, d)
  f <- function(x) exp(-(x - mu)^2 / (2 * sigma^2))
  stats::uniroot(function(a) sum((mu - a) * f(a) / (f(a) + f(-a))), c(0.1, 3), tol = 1e-13)$root
}

## d2S/dy2 at resident a of the GK98 pair: S = log sum_j f_j(y) / D_j, S' = 0.
gk98_pair_curvature <- function(a, d = 1.5, sigma = 1) {
  mu <- c(-d, 0, d)
  f <- function(x) exp(-(x - mu)^2 / (2 * sigma^2))
  w <- 1 / (f(a) + f(-a))
  sum(w * ((a - mu)^2 / sigma^4 - 1 / sigma^2) * f(a)) / sum(w * f(a))
}

gm99_estab <- function(y, x, N, alpha) {
  k <- lapply(N, function(n) 0:ceiling(n + 12 * sqrt(n + 1) + 30))
  grid <- as.matrix(expand.grid(k))
  w <- Reduce(`*`, lapply(seq_along(N), function(j) stats::dpois(grid[, j], N[j])))
  Z <- as.numeric(grid %*% exp(alpha * x))
  c_y <- exp(alpha * y)
  c(g0 = sum(w * c_y / (c_y + Z)), g1 = sum(w * c_y * Z / (c_y + Z)^2))
}
gm99_log_fitness <- function(y, x, N, alpha, beta, R = 1) {
  log(R / y) + log(1 - 2 * exp(-beta * y)) + log(gm99_estab(y, x, N, alpha)[["g0"]])
}
## residents' equilibrium densities: every resident's W = 1. Newton in log N
## by hand, since the root below runs nleqslv and nleqslv cannot nest.
gm99_coalition_density <- function(x, alpha, beta, R = 1, N0 = rep(1, length(x))) {
  f <- function(lN) vapply(x, gm99_log_fitness, 0, x = x, N = exp(lN), alpha = alpha, beta = beta, R = R)
  lN <- log(N0)
  for (it in 1:50) {
    f0 <- f(lN)
    if (max(abs(f0)) < 1e-13) break
    J <- vapply(seq_along(lN), function(j) (f(replace(lN, j, lN[j] + 1e-7)) - f0) / 1e-7, f0)
    step <- solve(J, f0)
    lN <- lN - step * min(1, 1 / max(abs(step)))
  }
  stopifnot(max(abs(f(lN))) < 1e-12)
  exp(lN)
}
## a mutant's gradient d log W / dy, residents x held at densities N
gm99_mutant_gradient <- function(y, x, N, alpha, beta) {
  g <- gm99_estab(y, x, N, alpha)
  -1 / y + 2 * beta * exp(-beta * y) / (1 - 2 * exp(-beta * y)) + alpha * g[["g1"]] / g[["g0"]]
}
## each resident's selection gradient: the mutant gradient at y = x_i
gm99_coalition_gradient <- function(x, alpha, beta, R = 1) {
  N <- gm99_coalition_density(x, alpha, beta, R)
  vapply(x, gm99_mutant_gradient, 0, x = x, N = N, alpha = alpha, beta = beta)
}
## each resident's curvature d2 log W / dy2 at y = x_i, the residents held
gm99_coalition_curvature <- function(x, alpha, beta, R = 1, h = 1e-5) {
  N <- gm99_coalition_density(x, alpha, beta, R)
  vapply(x, function(y) {
    (gm99_mutant_gradient(y + h, x, N, alpha, beta) -
       gm99_mutant_gradient(y - h, x, N, alpha, beta)) / (2 * h)
  }, 0)
}
gm99_pair_root <- function(alpha = 7, beta = 15, R = 1, x0 = c(0.2, 0.7)) {
  sol <- nleqslv::nleqslv(x0, gm99_coalition_gradient, alpha = alpha, beta = beta, R = R,
                          control = list(ftol = 1e-11, xtol = 1e-13))
  sol$x
}

## Lotka-Volterra for protected coexistence.
##
## lv_gauss_harness: a continuous trait, K(y) = exp(-y^2/2), competition
## A(y, x) = 1 + gamma (1 - exp(-(y - x)^2 / 2s^2)), stronger between residents
## than within them for gamma > 0 (priority effects), weaker for gamma < 0.
lv_gauss_harness <- function(gamma, s = 0.3) {
  A <- function(y, x) 1 + gamma * (1 - exp(-outer(y, x, "-")^2 / (2 * s^2)))
  K <- function(y) exp(-y^2 / 2)
  harness_explicit(
    fitness = function(x_mut, x_res, n_res, pars) 1 - as.numeric(A(x_mut, x_res) %*% n_res) / K(x_mut),
    equilibrium = function(x_res, pars) as.numeric(solve(A(x_res, x_res), K(x_res))),
    pars = list(), trait_names = "x", label = "lv")
}
## lv_matrix_harness: residents are the species 1, 2, ... of a fixed
## competition matrix A, the trait being the species' index (so only invasion
## at a resident's own trait is defined) and K = 1, so that a resident rare in
## a monoculture of j (n_j = 1) has fitness 1 - A_ij; the "model" equilibrium
## is solve(A, 1), negative densities and all.
lv_matrix_harness <- function(A) {
  idx <- function(x) as.integer(round(x))
  harness_explicit(
    fitness = function(x_mut, x_res, n_res, pars)
      1 - as.numeric(A[idx(x_mut), idx(x_res), drop = FALSE] %*% n_res),
    equilibrium = function(x_res, pars)
      as.numeric(solve(A[idx(x_res), idx(x_res), drop = FALSE], rep(1, length(x_res)))),
    pars = list(), trait_names = "x", label = "lv_matrix")
}
lv_matrix_community <- function(A) {
  m <- nrow(A)
  community_start(bounds(x = c(0, m + 1)), trait_scale = "linear", harness = lv_matrix_harness(A)) |>
    community_add(trait_matrix(seq_len(m), "x"), birth_rate = rep(1, m)) |>
    community_demography()
}

## Central-difference Jacobian of a vector function, for oracles built on
## closed forms (fine steps: the function is exact).
oracle_jacobian <- function(f, x, h = 1e-6) {
  vapply(seq_along(x), function(j) {
    e <- replace(numeric(length(x)), j, h)
    (f(x + e) - f(x - e)) / (2 * h)
  }, numeric(length(f(x))))
}
