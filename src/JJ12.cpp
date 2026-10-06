// JJ12 -- migratory-bird arrival-time model
// Johansson & Jonzen 2012, Ecol. Lett. 15:881-888; simplified analytic form of
// Brannstrom, Johansson & von Festenberg 2013, Games 4:304-328, section 4.
//
// Trait x = arrival time. Birds compete for K territories; early arrival raises
// competitive ability C(x) = exp(-a x), reproduction R(x) is Gaussian about the
// seasonal optimum x_opt, survival p in (0,1).
//
//   R(x) = R0 exp(-(x - x_opt)^2 / (2 sigma^2))
//   C(x) = exp(-a x)
//   resident demography (discrete):  n_{t+1} = K R(x) + p n_t
//   single-resident equilibrium:     n* = K R(x) / (1 - p)
//   invasion fitness of mutant x':    w = K R(x') C(x') / (sum_j n_j C(x_j)) + p
//   (jj12_fitness returns log w; for one resident at n* this reduces to
//    w = (1-p) R(x')C(x') / (R(x)C(x)) + p, and w(x,x) = 1 exactly.)
//
// Analytic oracle: singular strategy x* = x_opt - a*sigma^2 (a CSS, no branching).

#include <Rcpp.h>
#include <cmath>
#include <vector>
using namespace Rcpp;

namespace jj12 {

static inline double R_of(double x, double R0, double x_opt, double sigma) {
  double d = x - x_opt;
  return R0 * std::exp(-(d * d) / (2.0 * sigma * sigma));
}

static inline double C_of(double x, double a) {
  return std::exp(-a * x);
}

} // namespace jj12

//' JJ12 bird model: log invasion fitness of mutants
//'
//' @param x_mut numeric vector of mutant trait values (arrival times)
//' @param x_res numeric vector of resident trait values
//' @param n_res numeric vector of resident equilibrium densities
//' @param pars list with a, x_opt, sigma, R0, K, p
//' @return numeric vector of log invasion fitness, one per mutant
//' @keywords internal
// [[Rcpp::export]]
NumericVector jj12_fitness(NumericVector x_mut, NumericVector x_res,
                           NumericVector n_res, List pars) {
  double a = pars["a"], x_opt = pars["x_opt"], sigma = pars["sigma"],
         R0 = pars["R0"], K = pars["K"], p = pars["p"];
  int nm = x_mut.size();
  int nr = x_res.size();
  NumericVector out(nm);

  if (nr == 0) {
    // Fundamental fitness of a lone strategy: log of its equilibrium
    // reproduction (the bird model is viable everywhere R > 0); a sensible
    // default only, not used on the main path.
    for (int i = 0; i < nm; i++) {
      out[i] = std::log(K * jj12::R_of(x_mut[i], R0, x_opt, sigma) / (1.0 - p));
    }
    return out;
  }

  double denom = 0.0;
  for (int j = 0; j < nr; j++) denom += n_res[j] * jj12::C_of(x_res[j], a);

  for (int i = 0; i < nm; i++) {
    double w = K * jj12::R_of(x_mut[i], R0, x_opt, sigma) *
                   jj12::C_of(x_mut[i], a) / denom + p;
    out[i] = std::log(w);
  }
  return out;
}

//' JJ12 bird model: resident demographic equilibrium densities
//'
//' Single resident: closed form n* = K R(x)/(1-p). Multiple residents: iterate
//' the territory-competition recursion to its fixed point (a single limiting
//' resource, so this generically resolves to competitive exclusion of all but
//' the strategy maximising R(x)C(x)).
//'
//' @param x_res numeric vector of resident trait values
//' @param pars list with a, x_opt, sigma, R0, K, p
//' @param max_iter maximum fixed-point iterations (multi-resident case)
//' @param eps convergence tolerance (multi-resident case)
//' @return numeric vector of equilibrium densities, one per resident
//' @keywords internal
// [[Rcpp::export]]
NumericVector jj12_equilibrium(NumericVector x_res, List pars,
                               int max_iter = 5000, double eps = 1e-12) {
  double a = pars["a"], x_opt = pars["x_opt"], sigma = pars["sigma"],
         R0 = pars["R0"], K = pars["K"], p = pars["p"];
  int nr = x_res.size();
  NumericVector n(nr);
  if (nr == 0) return n;
  if (nr == 1) {
    n[0] = K * jj12::R_of(x_res[0], R0, x_opt, sigma) / (1.0 - p);
    return n;
  }

  std::vector<double> Rv(nr), Cv(nr);
  for (int j = 0; j < nr; j++) {
    Rv[j] = jj12::R_of(x_res[j], R0, x_opt, sigma);
    Cv[j] = jj12::C_of(x_res[j], a);
    n[j] = K * Rv[j] / (1.0 - p) / nr; // initial guess
  }

  for (int it = 0; it < max_iter; it++) {
    double denom = 0.0;
    for (int j = 0; j < nr; j++) denom += n[j] * Cv[j];
    double maxchange = 0.0;
    NumericVector nn(nr);
    for (int j = 0; j < nr; j++) {
      nn[j] = K * (n[j] * Cv[j] / denom) * Rv[j] + p * n[j];
      maxchange = std::max(maxchange, std::abs(nn[j] - n[j]));
    }
    n = nn;
    if (maxchange < eps) break;
  }
  return n;
}

// --- derivatives in the mutant direction -------------------------------------
// Against residents, w(y) = q(y) + p with q = K R(y) C(y) / denom, and
//   q' = q g,  g = -(y - x_opt)/sigma^2 - a,   q'' = q (g^2 - 1/sigma^2)
// so for S = log w:  S' = q g / w,  S'' = q (g^2 - 1/sigma^2) / w - S'^2.
// For a lone strategy S = log(K R(y)/(1-p)):  S' = -(y - x_opt)/sigma^2,
// S'' = -1/sigma^2.

namespace jj12 {

static void derivs(double y, double denom, bool lone, double a, double x_opt,
                   double sigma, double R0, double K, double p,
                   double& d1, double& d2) {
  double s2 = sigma * sigma;
  if (lone) {
    d1 = -(y - x_opt) / s2;
    d2 = -1.0 / s2;
    return;
  }
  double q = K * R_of(y, R0, x_opt, sigma) * C_of(y, a) / denom;
  double w = q + p;
  double g = -(y - x_opt) / s2 - a;
  d1 = q * g / w;
  d2 = q * (g * g - 1.0 / s2) / w - d1 * d1;
}

static double denominator(const NumericVector& x_res, const NumericVector& n_res,
                          double a) {
  double denom = 0.0;
  for (int j = 0; j < x_res.size(); j++) denom += n_res[j] * C_of(x_res[j], a);
  return denom;
}

} // namespace jj12

//' JJ12 bird model: gradient of log invasion fitness with respect to the mutant trait
//'
//' @inheritParams jj12_fitness
//' @return numeric matrix, one row per mutant and one column
//' @keywords internal
// [[Rcpp::export]]
NumericMatrix jj12_fitness_gradient(NumericVector x_mut, NumericVector x_res,
                                    NumericVector n_res, List pars) {
  double a = pars["a"], x_opt = pars["x_opt"], sigma = pars["sigma"],
         R0 = pars["R0"], K = pars["K"], p = pars["p"];
  bool lone = x_res.size() == 0;
  double denom = lone ? 0.0 : jj12::denominator(x_res, n_res, a);
  NumericMatrix out(x_mut.size(), 1);
  for (int i = 0; i < x_mut.size(); i++) {
    double d1, d2;
    jj12::derivs(x_mut[i], denom, lone, a, x_opt, sigma, R0, K, p, d1, d2);
    out(i, 0) = d1;
  }
  return out;
}

//' JJ12 bird model: second derivative of log invasion fitness with respect to the mutant trait
//'
//' @param x_mut a single mutant trait value
//' @inheritParams jj12_fitness
//' @return a 1 x 1 numeric matrix
//' @keywords internal
// [[Rcpp::export]]
NumericMatrix jj12_fitness_hessian(NumericVector x_mut, NumericVector x_res,
                                   NumericVector n_res, List pars) {
  if (x_mut.size() != 1) stop("jj12_fitness_hessian takes a single mutant");
  double a = pars["a"], x_opt = pars["x_opt"], sigma = pars["sigma"],
         R0 = pars["R0"], K = pars["K"], p = pars["p"];
  bool lone = x_res.size() == 0;
  double denom = lone ? 0.0 : jj12::denominator(x_res, n_res, a);
  double d1, d2;
  jj12::derivs(x_mut[0], denom, lone, a, x_opt, sigma, R0, K, p, d1, d2);
  NumericMatrix out(1, 1);
  out(0, 0) = d2;
  return out;
}
