// GK98 -- Geritz, Kisdi, Meszena & Metz 1998 soft-selection model
// "Evolutionarily singular strategies and the adaptive growth and branching of
// the evolutionary tree", Evol. Ecol. 12:35-57 (the worked Levene example).
//
// m patches, each with capacity K_j and Gaussian within-patch survival peaked at
// optimum mu_j:
//   f_j(x) = exp(-(x - mu_j)^2 / (2 sigma^2))        (amplitude cancels)
//   resident equilibrium (Eq 19):  N'_i = sum_j K_j f_j(x_i) N_i / (sum_h f_j(x_h) N_h)
//   invasion fitness (Eq 21):  S(y) = log sum_j K_j f_j(y) / (sum_h f_j(x_h) N_h)
//   (single resident reduces to Eq B1:
//        S(y) = log sum_j c_j f_j(y)/f_j(x),  c_j = K_j / sum K)
//
// Analytic oracles (symmetric m=3, mu=(-d,0,d), equal K): x* = sum c_j mu_j = 0;
// convergence stable for all d/sigma; branching point (lacks ESS) iff
// d/sigma > sqrt(3/2) ~= 1.2247, a CSS below that.

#include <Rcpp.h>
#include <cmath>
#include <vector>
using namespace Rcpp;

namespace gk98 {

static inline double f_of(double x, double mu, double sigma) {
  double d = x - mu;
  return std::exp(-(d * d) / (2.0 * sigma * sigma));
}

} // namespace gk98

//' GK98 soft-selection model: log invasion fitness of mutants
//'
//' @param x_mut numeric vector of mutant trait values
//' @param x_res numeric vector of resident trait values
//' @param n_res numeric vector of resident equilibrium abundances
//' @param pars list with sigma, mu (patch optima), K (patch capacities)
//' @return numeric vector of log invasion fitness (=0 for resident at equilibrium)
//' @keywords internal
// [[Rcpp::export]]
NumericVector gk98_fitness(NumericVector x_mut, NumericVector x_res,
                           NumericVector n_res, List pars) {
  double sigma = pars["sigma"];
  NumericVector mu = pars["mu"], Kp = pars["K"];
  int m = mu.size();
  int nm = x_mut.size();
  int nr = x_res.size();
  NumericVector out(nm);

  if (nr == 0) {
    // fundamental fitness of a lone strategy: log sum_j c_j f_j(y)
    double Ktot = 0.0;
    for (int j = 0; j < m; j++) Ktot += Kp[j];
    for (int i = 0; i < nm; i++) {
      double s = 0.0;
      for (int j = 0; j < m; j++)
        s += (Kp[j] / Ktot) * gk98::f_of(x_mut[i], mu[j], sigma);
      out[i] = std::log(s);
    }
    return out;
  }

  // patch denominators D_j = sum_h f_j(x_h) N_h
  std::vector<double> D(m, 0.0);
  for (int j = 0; j < m; j++)
    for (int h = 0; h < nr; h++)
      D[j] += gk98::f_of(x_res[h], mu[j], sigma) * n_res[h];

  for (int i = 0; i < nm; i++) {
    double s = 0.0;
    for (int j = 0; j < m; j++)
      s += Kp[j] * gk98::f_of(x_mut[i], mu[j], sigma) / D[j];
    out[i] = std::log(s);
  }
  return out;
}

//' GK98 soft-selection model: resident equilibrium abundances
//'
//' Iterates the soft-selection recursion (Eq 19) to its fixed point; total
//' abundance is conserved at sum(K). Single resident: N* = sum(K).
//'
//' @param x_res numeric vector of resident trait values
//' @param pars list with sigma, mu (patch optima), K (patch capacities)
//' @param max_iter maximum fixed-point iterations
//' @param eps convergence tolerance
//' @return numeric vector of equilibrium abundances
//' @keywords internal
// [[Rcpp::export]]
NumericVector gk98_equilibrium(NumericVector x_res, List pars,
                               int max_iter = 5000, double eps = 1e-12) {
  double sigma = pars["sigma"];
  NumericVector mu = pars["mu"], Kp = pars["K"];
  int m = mu.size();
  int nr = x_res.size();
  NumericVector n(nr);
  if (nr == 0) return n;
  double Ktot = 0.0;
  for (int j = 0; j < m; j++) Ktot += Kp[j];
  if (nr == 1) {
    n[0] = Ktot;
    return n;
  }
  for (int h = 0; h < nr; h++) n[h] = Ktot / nr; // initial guess

  std::vector<std::vector<double> > F(m, std::vector<double>(nr));
  for (int j = 0; j < m; j++)
    for (int h = 0; h < nr; h++) F[j][h] = gk98::f_of(x_res[h], mu[j], sigma);

  for (int it = 0; it < max_iter; it++) {
    std::vector<double> D(m, 0.0);
    for (int j = 0; j < m; j++)
      for (int h = 0; h < nr; h++) D[j] += F[j][h] * n[h];
    double maxchange = 0.0;
    NumericVector nn(nr);
    for (int h = 0; h < nr; h++) {
      double v = 0.0;
      for (int j = 0; j < m; j++) v += Kp[j] * F[j][h] * n[h] / D[j];
      nn[h] = v;
      maxchange = std::max(maxchange, std::abs(nn[h] - n[h]));
    }
    n = nn;
    if (maxchange < eps) break;
  }
  return n;
}

// --- derivatives in the mutant direction -------------------------------------
// With S(y) = log A(y), A = sum_j w_j f_j(y) and f_j'(y) = f_j b_j,
// b_j = -(y - mu_j)/sigma^2, f_j'' = f_j (b_j^2 - 1/sigma^2):
//   S'  = sum_j w_j f_j b_j / A
//   S'' = sum_j w_j f_j (b_j^2 - 1/sigma^2) / A - S'^2
// where w_j = K_j / D_j against residents, or c_j for a lone strategy.

namespace gk98 {

static void weights(const NumericVector& x_res, const NumericVector& n_res,
                    const NumericVector& mu, const NumericVector& Kp,
                    double sigma, std::vector<double>& w) {
  int m = mu.size(), nr = x_res.size();
  w.assign(m, 0.0);
  if (nr == 0) {
    double Ktot = 0.0;
    for (int j = 0; j < m; j++) Ktot += Kp[j];
    for (int j = 0; j < m; j++) w[j] = Kp[j] / Ktot;
    return;
  }
  for (int j = 0; j < m; j++) {
    double D = 0.0;
    for (int h = 0; h < nr; h++) D += f_of(x_res[h], mu[j], sigma) * n_res[h];
    w[j] = Kp[j] / D;
  }
}

static void derivs(double y, const std::vector<double>& w, const NumericVector& mu,
                   double sigma, double& d1, double& d2) {
  int m = mu.size();
  double A = 0.0, A1 = 0.0, A2 = 0.0;
  for (int j = 0; j < m; j++) {
    double f = f_of(y, mu[j], sigma);
    double b = -(y - mu[j]) / (sigma * sigma);
    A  += w[j] * f;
    A1 += w[j] * f * b;
    A2 += w[j] * f * (b * b - 1.0 / (sigma * sigma));
  }
  d1 = A1 / A;
  d2 = A2 / A - d1 * d1;
}

} // namespace gk98

//' GK98 soft-selection model: gradient of log invasion fitness with respect to the mutant trait
//'
//' @inheritParams gk98_fitness
//' @return numeric matrix, one row per mutant and one column
//' @keywords internal
// [[Rcpp::export]]
NumericMatrix gk98_fitness_gradient(NumericVector x_mut, NumericVector x_res,
                                    NumericVector n_res, List pars) {
  double sigma = pars["sigma"];
  NumericVector mu = pars["mu"], Kp = pars["K"];
  std::vector<double> w;
  gk98::weights(x_res, n_res, mu, Kp, sigma, w);
  NumericMatrix out(x_mut.size(), 1);
  for (int i = 0; i < x_mut.size(); i++) {
    double d1, d2;
    gk98::derivs(x_mut[i], w, mu, sigma, d1, d2);
    out(i, 0) = d1;
  }
  return out;
}

//' GK98 soft-selection model: second derivative of log invasion fitness with respect to the mutant trait
//'
//' @param x_mut a single mutant trait value
//' @inheritParams gk98_fitness
//' @return a 1 x 1 numeric matrix
//' @keywords internal
// [[Rcpp::export]]
NumericMatrix gk98_fitness_hessian(NumericVector x_mut, NumericVector x_res,
                                   NumericVector n_res, List pars) {
  if (x_mut.size() != 1) stop("gk98_fitness_hessian takes a single mutant");
  double sigma = pars["sigma"];
  NumericVector mu = pars["mu"], Kp = pars["K"];
  std::vector<double> w;
  gk98::weights(x_res, n_res, mu, Kp, sigma, w);
  double d1, d2;
  gk98::derivs(x_mut[0], w, mu, sigma, d1, d2);
  NumericMatrix out(1, 1);
  out(0, 0) = d2;
  return out;
}
