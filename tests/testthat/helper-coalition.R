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

## Central-difference Jacobian of a vector function, for oracles built on
## closed forms (fine steps: the function is exact).
oracle_jacobian <- function(f, x, h = 1e-6) {
  vapply(seq_along(x), function(j) {
    e <- replace(numeric(length(x)), j, h)
    (f(x + e) - f(x - e)) / (2 * h)
  }, numeric(length(f(x))))
}
