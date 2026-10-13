# Closed-form oracles for coexistence boundaries in parameter space, computed
# independently of the package.
#
# GK98 (Levene soft selection) with three patches mu = (-d, 0, d), capacities
# K = (1, kappa, 1) and sigma = 1, in the (d, kappa) plane: the number of
# species in the evolutionarily stable community is 1, 2 or 3.
#
#   1 -> 2: the generalist at 0 is a branching point iff S''(0) > 0 with
#           S''(0) = sum_j c_j mu_j^2 - 1 = 2 d^2 / (2 + kappa) - 1, so the
#           boundary is d = sqrt((2 + kappa) / 2), a pitchfork (the pair +-a*
#           emerges from 0 continuously).
#   2 -> 3: the mirror-image pair +-a* holds equal shares (N = K_tot / 2
#           each) and each resident's gradient vanishes,
#               sum_j K_j (mu_j - a) f_j(a) / (f_j(a) + f_j(-a)) = 0,
#           so a* is a scalar root; the boundary is where a mutant at the
#           middle patch, y* = 0 by symmetry, is exactly marginal against the
#           pair, S(0; +-a*) = 0, a scalar root in d. There y* = 0 is a strict
#           maximum of S distinct from the residents, which are each at an ESS,
#           and the third species' density grows from zero beyond it: the
#           transversal case. Three is the most (one per patch).
#
# DD99 (Gaussian-Gaussian Lotka-Volterra) has one boundary, sigma_C =
# sigma_K: below it the pair's curvature is positive and S(0) > 0 for every
# ratio, and fitness tends to r far from the residents, so no finite
# coalition resists invasion and the box's edge sets the count.

gk98_patch <- function(x, mu, sigma = 1) exp(-(x - mu)^2 / (2 * sigma^2))

## the mirror-image pair's half-separation a* at spacing d
gk98_kappa_pair_root <- function(d, kappa = 1) {
  mu <- c(-d, 0, d); K <- c(1, kappa, 1)
  g <- function(a) sum(K * (mu - a) * gk98_patch(a, mu) / (gk98_patch(a, mu) + gk98_patch(-a, mu)))
  stats::uniroot(g, c(1e-8, 3 * d), tol = 1e-14)$root
}

## invasion fitness of mutants y against the pair +-a, each at K_tot / 2
gk98_kappa_pair_fitness <- function(y, a, d, kappa = 1) {
  mu <- c(-d, 0, d); K <- c(1, kappa, 1)
  D <- (gk98_patch(a, mu) + gk98_patch(-a, mu)) * sum(K) / 2
  vapply(y, function(yy) log(sum(K * gk98_patch(yy, mu) / D)), 0)
}

## the two boundaries along d at a given kappa
gk98_kappa_boundaries <- function(kappa = 1) {
  c1 <- sqrt((2 + kappa) / 2)
  c2 <- stats::uniroot(function(d) gk98_kappa_pair_fitness(0, gk98_kappa_pair_root(d, kappa), d, kappa),
                       c(c1 + 1e-3, 6), tol = 1e-13)$root
  c(c1 = c1, c2 = c2)
}

## the symmetric triple (-b, 0, b): densities by iterating the soft-selection
## recursion (Eq. 19), b the root of the outer resident's gradient
gk98_kappa_triple <- function(d, kappa = 1) {
  mu <- c(-d, 0, d); K <- c(1, kappa, 1)
  density <- function(x) {
    F <- outer(mu, x, function(m, xx) gk98_patch(xx, m))
    n <- rep(sum(K) / 3, 3)
    for (i in seq_len(20000)) {
      nn <- as.numeric(t(F) %*% (K / as.numeric(F %*% n))) * n
      if (max(abs(nn - n)) < 1e-14) break
      n <- nn
    }
    n
  }
  S <- function(y, x, n) {
    D <- as.numeric(outer(mu, x, function(m, xx) gk98_patch(xx, m)) %*% n)
    log(sum(K * gk98_patch(y, mu) / D))
  }
  outer_gradient <- function(b) {
    x <- c(-b, 0, b); n <- density(x); h <- 1e-7
    (S(b + h, x, n) - S(b - h, x, n)) / (2 * h)
  }
  b <- stats::uniroot(outer_gradient, c(0.3, 2 * d), tol = 1e-12)$root
  list(traits = c(-b, 0, b), density = density(c(-b, 0, b)))
}

## the GK98 community the label starts from: seeded, since the empty
## soft-selection community has fitness <= 0 everywhere
gk98_kappa_start <- function(d, kappa = 1, extinct = 1e-8, ...) {
  community_start(bounds(x = c(-4, 4)), trait_scale = "linear",
                  harness = harness_gk98(sigma = 1, mu = c(-d, 0, d), K = c(1, kappa, 1)),
                  demography_control = demographic_step_control(list(equilibrium_extinct_birth_rate = extinct)),
                  ...) |>
    community_add(trait_matrix(0.5, "x"), birth_rate = 1)
}

## DD99 seeded at x = 0.5: its empty community's fitness is flat
dd99_start <- function(sigma_C = 0.4, sigma_K = 1, ...) {
  community_start(bounds(x = c(-2, 2)), trait_scale = "linear",
                  harness = harness_dd99(sigma_C = sigma_C, sigma_K = sigma_K), ...) |>
    community_add(trait_matrix(0.5, "x"), birth_rate = 100)
}

## a singular strategy the dynamics flees: s(y; x) = c (y - x) x - (y - x)^2
## has x* = 0 with H = -2 (an ESS) and J = c > 0 (a Garden of Eden)
garden_of_eden_start <- function(c = 1) {
  h <- harness_explicit(
    fitness = function(x_mut, x_res, n_res, pars) pars$c * (x_mut - x_res[1]) * x_res[1] - (x_mut - x_res[1])^2,
    equilibrium = function(x_res, pars) rep(1, length(x_res)),
    pars = list(c = c), trait_names = "x", label = "goe")
  community_start(bounds(x = c(-1, 1)), trait_scale = "linear", harness = h) |>
    community_add(trait_matrix(0.3, "x"), birth_rate = 1)
}

## Two resource peaks at -1 and +1 of heights 1 - q and 1 + q, and
## competition between species stronger than within them (gamma > 0), so a
## monoculture on one peak resists a newcomer on the other: alternative
## stable communities. A monoculture at x (a CSS, at its peak) has
## N = K(x), and a newcomer at y has fitness 1 - A(y, x) K(x) / K(y), with
## A(y, x) ~ 1 + gamma between the peaks. So the monoculture on the -1 peak
## resists invasion iff (1 - q) / (1 + q) > 1 / (1 + gamma), that on +1 iff
## (1 + q) / (1 - q) > 1 / (1 + gamma): both, with gamma = 1, for
## |q| < 1/3. Equilibria by iterating the one-generation map, so a newcomer
## that wins excludes the resident as the dynamics would.
bistable_start <- function(q = 0, gamma = 1, s = 0.3, w = 0.3, x = 0.5) {
  h <- harness_explicit(
    fitness = function(x_mut, x_res, n_res, pars) {
      K <- function(y) (1 + pars$q) * exp(-(y - 1)^2 / (2 * pars$s^2)) +
        (1 - pars$q) * exp(-(y + 1)^2 / (2 * pars$s^2))
      A <- 1 + pars$gamma * (1 - exp(-outer(x_mut, x_res, "-")^2 / (2 * pars$w^2)))
      1 - as.numeric(A %*% n_res) / K(x_mut)
    },
    equilibrium = function(x_res, pars) stop("iterate the map"),
    pars = list(q = q, gamma = gamma, s = s, w = w), trait_names = "x", label = "bistable")
  community_start(bounds(x = c(-2, 2)), trait_scale = "linear", harness = h,
                  demography_control = demographic_step_control(list(
                    equilibrium_solver_name = "equilibrium_iteration", equilibrium_nsteps = 500,
                    equilibrium_eps = 1e-10))) |>
    community_add(trait_matrix(x, "x"), birth_rate = 0.5)
}

## A carrying capacity with a dip of depth h at y = 1.5, and a narrow
## competition kernel, so the landscape of a resident at 0 has a minimum
## near the dip that crosses zero as the dip fills: an invasion-type
## boundary point at which y* is a minimum, not a peak.
dip_start <- function(h = 0.95) {
  K <- function(y, h) exp(-y^2 / 2) * (1 - h * exp(-(y - 1.5)^2 / (2 * 0.1^2)))
  hr <- harness_explicit(
    fitness = function(x_mut, x_res, n_res, pars) {
      C <- exp(-outer(x_mut, x_res, "-")^2 / (2 * 0.5^2))
      1 - as.numeric(C %*% n_res) / K(x_mut, pars$h)
    },
    equilibrium = function(x_res, pars) K(x_res, pars$h),
    pars = list(h = h), trait_names = "x", label = "dip")
  community_start(bounds(x = c(-3, 3)), trait_scale = "linear", harness = hr) |>
    community_add(trait_matrix(0, "x"), birth_rate = 1) |>
    community_demography()
}

## GK98 through (d, kappa) as a parameter map. It runs on the workers of a
## plan, which load the package and not these helpers, so it is written with
## package functions only.
gk98_by_d_kappa <- function(p) {
  community_start(bounds(x = c(-4, 4)), trait_scale = "linear",
                  harness = harness_gk98(sigma = 1, mu = c(-p[1], 0, p[1]), K = c(1, p[2], 1)),
                  demography_control = demographic_step_control(list(equilibrium_extinct_birth_rate = 1e-8))) |>
    community_add(trait_matrix(0.5, "x"), birth_rate = 1)
}

## Where the GK98 boundaries end. The single resident at 0 has
##   S(y) = -y^2 / 2 + log(sum_j c_j exp(y mu_j)),
## minus y^2/2 plus the cumulant generating function of the optima under the
## capacities c. On the 1 -> 2 locus the second cumulant is 1, so the quartic
## term is the fourth cumulant, (2 + kappa) / 2 - 3: the pitchfork is
## supercritical (a pair emerges) for kappa < 4 and subcritical beyond, where
## the locus stops being a coexistence boundary and the 2 -> 3 boundary ends
## on it, at (sqrt(3), 4).
gk98_locus_fitness <- function(y, kappa) {
  d <- sqrt((2 + kappa) / 2)
  -y^2 / 2 + log((2 * cosh(y * d) + kappa) / (2 + kappa))
}
## the 2 -> 3 boundary's pair reaches half-separation `a` (where a
## continuation's separation test fires) at this kappa
gk98_merge_kappa <- function(a) {
  stats::uniroot(function(k) gk98_kappa_pair_root(gk98_kappa_boundaries(k)[["c2"]], k) - a,
                 c(3.5, 3.98), tol = 1e-13)$root
}
## on the 1 -> 2 locus, the highest fitness farther than `reach` from the
## resident reaches `fitness_tol` at this kappa
gk98_subcritical_kappa <- function(reach, fitness_tol) {
  high <- function(k) {
    y <- seq(reach, 4, length.out = 4001)
    max(gk98_locus_fitness(y, k))
  }
  stats::uniroot(function(k) high(k) - fitness_tol, c(4 + 1e-6, 4.5), tol = 1e-13)$root
}
