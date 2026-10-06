// GM99 -- Geritz, van der Meijden & Metz 1999 seed-size safe-site model
// "Evolutionary Dynamics of Seed Size and Seedling Competitive Ability."
//
// Plants compete for safe sites; seeds arrive by a Poisson process; seedlings
// undergo size-asymmetric lottery competition within each site.
//   trait x = seed size, 0 < x <= R (R = resources per site)
//   pre-competitive survival   s(x) = max(0, 1 - 2 exp(-beta x))
//   competitive ability        c(x) = exp(alpha x)
//   establishment of mutant x' against residents {x_j} at densities {N_j}:
//     g = E[ c(x') / (c(x') + sum_j k_j c(x_j)) ],  k_j ~ Poisson(N_j) indep.
//   invasion fitness (lifetime reproductive ratio, = 1 at equilibrium):
//     W = (R/x') s(x') g          (gm99_fitness returns log W; resident ~0)
//   resident equilibrium: W_i = 1 for every resident; single resident solves
//     (R/x) s(x) (1 - exp(-N))/N = 1   (exists iff (R/x) s(x) > 1).
//
// Only the products alpha*R and beta*R matter. No closed-form singular strategy
// (numerical); branching diagnosed by the sign of the invasion-fitness second
// derivative. This is the model implemented in Daniel's MATLAB.
// Distinct from the GK98 1998 soft-selection model.

#include <Rcpp.h>
#include <cmath>
#include <vector>
using namespace Rcpp;

namespace gm99 {

static inline double s_of(double x, double beta) {
  double v = 1.0 - 2.0 * std::exp(-beta * x);
  return v > 0.0 ? v : 0.0;
}
static inline double c_of(double x, double alpha) {
  return std::exp(alpha * x);
}

// Expectations over independent k_j ~ Poisson(N_j), with Z = sum_j k_j cres_j:
//   g0 = E[ c / (c + Z) ]           establishment probability of the mutant
//   g1 = E[ c Z / (c + Z)^2 ]       (d/dy) g = alpha g1,  since dc/dy = alpha c
//   g2 = E[ c Z (Z - c) / (c + Z)^3 ]   (d2/dy2) g = alpha^2 g2
//
// Cartesian odometer over per-resident truncated Poisson pmfs. Truncation at
// N + 6*sqrt(N+1) + 6 captures > 1-1e-10 of the Poisson mass (tighter than
// the previous N + 10*sqrt(N+1) + 20, giving ~2x fewer terms per dimension).
static void estab_moments(double c_mut, const std::vector<double>& cres,
                          const std::vector<double>& N,
                          double& g0, double& g1, double& g2) {
  int n = cres.size();
  g0 = 1.0; g1 = 0.0; g2 = 0.0;
  if (n == 0) return;
  std::vector<std::vector<double> > P(n);
  std::vector<int> Kmax(n);
  for (int j = 0; j < n; j++) {
    int kmax = (int)std::ceil(N[j] + 6.0 * std::sqrt(N[j] + 1.0) + 6.0);
    P[j].resize(kmax + 1);
    double pk = std::exp(-N[j]);
    for (int k = 0; k <= kmax; k++) { P[j][k] = pk; pk *= N[j] / (k + 1); }
    Kmax[j] = kmax;
  }
  std::vector<int> k(n, 0);
  double sum = 0.0, sum1 = 0.0, sum2 = 0.0;
  while (true) {
    double w = 1.0, denom = c_mut;
    for (int j = 0; j < n; j++) { w *= P[j][k[j]]; denom += k[j] * cres[j]; }
    double Z = denom - c_mut;
    sum  += w * c_mut / denom;
    sum1 += w * c_mut * Z / (denom * denom);
    sum2 += w * c_mut * Z * (Z - c_mut) / (denom * denom * denom);
    int j = 0;
    while (j < n) { if (++k[j] <= Kmax[j]) break; k[j] = 0; j++; }
    if (j == n) break;
  }
  g0 = sum; g1 = sum1; g2 = sum2;
}

static double estab(double c_mut, const std::vector<double>& cres,
                    const std::vector<double>& N) {
  double g0, g1, g2;
  estab_moments(c_mut, cres, N, g0, g1, g2);
  return g0;
}

// log W = log R - log y + log s(y) + log g(y); first and second derivatives in
// the mutant seed size y. NaN where the mutant is not viable (s(y) <= 0).
static void derivs(double y, const std::vector<double>& cres,
                   const std::vector<double>& N, double R, double alpha,
                   double beta, double& d1, double& d2) {
  double s = s_of(y, beta);
  if (s <= 0.0 || y <= 0.0) { d1 = NA_REAL; d2 = NA_REAL; return; }
  double e = std::exp(-beta * y);
  double s1 = 2.0 * beta * e, s2 = -2.0 * beta * beta * e;
  double g0, g1, g2;
  estab_moments(c_of(y, alpha), cres, N, g0, g1, g2);
  double dlg = alpha * g1 / g0;
  d1 = -1.0 / y + s1 / s + dlg;
  d2 = 1.0 / (y * y) + (s2 * s - s1 * s1) / (s * s)
     + alpha * alpha * g2 / g0 - dlg * dlg;
}

// Evaluate log fitness for each resident at densities Nv.
static void eval_logW(const std::vector<double>& seeds_i,
                      const std::vector<double>& cres,
                      const std::vector<double>& Nv,
                      std::vector<double>& logW) {
  int nr = cres.size();
  for (int i = 0; i < nr; i++) {
    double W = seeds_i[i] * estab(cres[i], cres, Nv);
    logW[i] = std::log(W > 1e-300 ? W : 1e-300);
  }
}

// Solve J * delta = rhs by Gaussian elimination with partial pivoting (in-place).
// J is nr x nr, rhs is length nr; result overwrites rhs.
static void gauss_solve(std::vector<std::vector<double> >& A,
                        std::vector<double>& rhs, int nr) {
  for (int col = 0; col < nr; col++) {
    int pivot = col;
    for (int row = col + 1; row < nr; row++)
      if (std::abs(A[row][col]) > std::abs(A[pivot][col])) pivot = row;
    std::swap(A[col], A[pivot]);
    std::swap(rhs[col], rhs[pivot]);
    double d = A[col][col];
    if (std::abs(d) < 1e-14) continue;
    for (int row = col + 1; row < nr; row++) {
      double f = A[row][col] / d;
      for (int k = col; k < nr; k++) A[row][k] -= f * A[col][k];
      rhs[row] -= f * rhs[col];
    }
  }
  for (int i = nr - 1; i >= 0; i--) {
    if (std::abs(A[i][i]) < 1e-14) continue;
    for (int j = i + 1; j < nr; j++) rhs[i] -= A[i][j] * rhs[j];
    rhs[i] /= A[i][i];
  }
}

} // namespace gm99

//' GM99 seed-size model: log invasion fitness of mutants
//'
//' @param x_mut numeric vector of mutant seed sizes
//' @param x_res numeric vector of resident seed sizes
//' @param n_res numeric vector of resident equilibrium densities (seeds/site)
//' @param pars list with R, alpha, beta
//' @return numeric vector of log invasion fitness (=0 for resident at equilibrium;
//'   -Inf for non-viable mutants where s(x) = 0)
//' @keywords internal
// [[Rcpp::export]]
NumericVector gm99_fitness(NumericVector x_mut, NumericVector x_res,
                           NumericVector n_res, List pars) {
  double R = pars["R"], alpha = pars["alpha"], beta = pars["beta"];
  int nm = x_mut.size();
  int nr = x_res.size();
  NumericVector out(nm);
  std::vector<double> cres(nr), N(nr);
  for (int j = 0; j < nr; j++) { cres[j] = gm99::c_of(x_res[j], alpha); N[j] = n_res[j]; }
  for (int i = 0; i < nm; i++) {
    double mp = x_mut[i];
    double s = gm99::s_of(mp, beta);
    if (s <= 0.0 || mp <= 0.0) { out[i] = R_NegInf; continue; }
    double seeds = R / mp * s;
    double g = (nr == 0) ? 1.0 : gm99::estab(gm99::c_of(mp, alpha), cres, N);
    out[i] = std::log(seeds * g);
  }
  return out;
}

//' GM99 seed-size model: resident equilibrium densities (seeds/site)
//'
//' Single resident: solve (R/x) s(x) (1 - exp(-N))/N = 1 by bisection (N* = 0 if
//' the strategy is non-viable). Many residents: Newton's method on log W_i = 0
//' in log-N space, with a numerical Jacobian. Newton converges quadratically so
//' only ~10-30 steps are needed regardless of resident count, vs. thousands for
//' the previous fixed-point iteration.
//'
//' @param x_res numeric vector of resident seed sizes
//' @param pars list with R, alpha, beta
//' @param max_iter maximum Newton iterations (multi-resident case)
//' @param eps convergence tolerance on max |log W_i|
//' @return numeric vector of equilibrium densities
//' @keywords internal
// [[Rcpp::export]]
NumericVector gm99_equilibrium(NumericVector x_res, List pars,
                               int max_iter = 200, double eps = 1e-10) {
  double R = pars["R"], alpha = pars["alpha"], beta = pars["beta"];
  int nr = x_res.size();
  NumericVector n(nr);
  if (nr == 0) return n;
  if (nr == 1) {
    double m = x_res[0], s = gm99::s_of(m, beta);
    double seeds = (s > 0.0 && m > 0.0) ? R / m * s : 0.0;
    if (seeds <= 1.0) { n[0] = 0.0; return n; }
    // (R/x) s(x) (1 - e^-N)/N = 1 is monotone decreasing in N; bisect [0, seeds]
    double lo = 1e-12, hi = seeds;
    for (int it = 0; it < 100; it++) {
      double mid = 0.5 * (lo + hi);
      double f = seeds * (1.0 - std::exp(-mid)) / mid - 1.0;
      if (f > 0.0) lo = mid; else hi = mid;
    }
    n[0] = 0.5 * (lo + hi);
    return n;
  }

  // Multi-resident: Newton's method on F(log N) = log W(N) = 0.
  // Each step costs na*(na+1) estab evaluations (na = active residents) and
  // converges quadratically, so only ~10-30 steps are typically needed.
  // Residents with log W_i < 0 when N_i is negligible are flagged as
  // excluded and held at 0 to avoid Jacobian singularities.
  std::vector<double> cres(nr), seeds_i(nr);
  for (int j = 0; j < nr; j++) {
    cres[j] = gm99::c_of(x_res[j], alpha);
    NumericVector xj(1); xj[0] = x_res[j];
    n[j] = gm99_equilibrium(xj, pars)[0];
    if (n[j] <= 0.0) n[j] = 0.01;
    double s = gm99::s_of(x_res[j], beta);
    seeds_i[j] = (s > 0.0 && x_res[j] > 0.0) ? R / x_res[j] * s : 0.0;
  }

  const double h = 1e-5;          // log-N perturbation for numerical Jacobian
  const double N_extinct = 1e-8;  // treat residents below this as excluded
  std::vector<bool> active(nr, true);
  std::vector<double> Nv(nr), logW(nr), logW_plus(nr);
  std::vector<std::vector<double> > J(nr, std::vector<double>(nr));

  for (int it = 0; it < max_iter; it++) {
    // Build active index map and pack N into Nv (excluded held at 0)
    for (int j = 0; j < nr; j++) Nv[j] = active[j] ? n[j] : 0.0;

    gm99::eval_logW(seeds_i, cres, Nv, logW);

    // Flag any active resident with near-zero N and negative log W as excluded
    bool changed = false;
    for (int i = 0; i < nr; i++) {
      if (active[i] && n[i] < N_extinct && logW[i] < 0.0) {
        active[i] = false; n[i] = 0.0; Nv[i] = 0.0; changed = true;
      }
    }
    if (changed) gm99::eval_logW(seeds_i, cres, Nv, logW);

    // Convergence check on active residents only
    double maxlogW = 0.0;
    for (int i = 0; i < nr; i++)
      if (active[i]) maxlogW = std::max(maxlogW, std::abs(logW[i]));
    if (maxlogW < eps) break;

    // Build active index list
    std::vector<int> idx;
    for (int j = 0; j < nr; j++) if (active[j]) idx.push_back(j);
    int na = idx.size();

    // Numerical Jacobian restricted to active residents: J_ij = d(logW_i)/d(logN_j)
    std::vector<std::vector<double> > Ja(na, std::vector<double>(na));
    std::vector<double> rhs(na);
    for (int jj = 0; jj < na; jj++) {
      int j = idx[jj];
      Nv[j] = n[j] * (1.0 + h);
      gm99::eval_logW(seeds_i, cres, Nv, logW_plus);
      for (int ii = 0; ii < na; ii++) Ja[ii][jj] = (logW_plus[idx[ii]] - logW[idx[ii]]) / h;
      Nv[j] = n[j];
    }
    for (int ii = 0; ii < na; ii++) rhs[ii] = -logW[idx[ii]];

    // Solve Ja * delta = rhs; apply log-N step with |delta| <= 2 guard
    gm99::gauss_solve(Ja, rhs, na);
    double scale = 1.0;
    for (int ii = 0; ii < na; ii++) if (std::abs(rhs[ii]) * scale > 2.0) scale = 2.0 / std::abs(rhs[ii]);

    for (int ii = 0; ii < na; ii++) {
      int i = idx[ii];
      double newN = n[i] * std::exp(scale * rhs[ii]);
      n[i] = newN < 1e-30 ? 1e-30 : newN;
    }
  }
  return n;
}

//' GM99 seed-size model: gradient of log invasion fitness with respect to the mutant seed size
//'
//' @inheritParams gm99_fitness
//' @return numeric matrix, one row per mutant and one column (NaN for non-viable mutants)
//' @keywords internal
// [[Rcpp::export]]
NumericMatrix gm99_fitness_gradient(NumericVector x_mut, NumericVector x_res,
                                    NumericVector n_res, List pars) {
  double R = pars["R"], alpha = pars["alpha"], beta = pars["beta"];
  int nr = x_res.size();
  std::vector<double> cres(nr), N(nr);
  for (int j = 0; j < nr; j++) { cres[j] = gm99::c_of(x_res[j], alpha); N[j] = n_res[j]; }
  NumericMatrix out(x_mut.size(), 1);
  for (int i = 0; i < x_mut.size(); i++) {
    double d1, d2;
    gm99::derivs(x_mut[i], cres, N, R, alpha, beta, d1, d2);
    out(i, 0) = d1;
  }
  return out;
}

//' GM99 seed-size model: second derivative of log invasion fitness with respect to the mutant seed size
//'
//' @param x_mut a single mutant seed size
//' @inheritParams gm99_fitness
//' @return a 1 x 1 numeric matrix (NaN for a non-viable mutant)
//' @keywords internal
// [[Rcpp::export]]
NumericMatrix gm99_fitness_hessian(NumericVector x_mut, NumericVector x_res,
                                   NumericVector n_res, List pars) {
  if (x_mut.size() != 1) stop("gm99_fitness_hessian takes a single mutant");
  double R = pars["R"], alpha = pars["alpha"], beta = pars["beta"];
  int nr = x_res.size();
  std::vector<double> cres(nr), N(nr);
  for (int j = 0; j < nr; j++) { cres[j] = gm99::c_of(x_res[j], alpha); N[j] = n_res[j]; }
  double d1, d2;
  gm99::derivs(x_mut[0], cres, N, R, alpha, beta, d1, d2);
  NumericMatrix out(1, 1);
  out(0, 0) = d2;
  return out;
}
