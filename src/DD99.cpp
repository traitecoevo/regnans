// DD99 -- Dieckmann & Doebeli 1999 competition model
// "On the origin of species by sympatric speciation", Nature 400:354-357.
//
// Continuous-time logistic competition for a Gaussian resource:
//   K(x)   = K0 exp(-(x - x0)^2 / (2 sigma_K^2))    carrying capacity
//   C(d)   = exp(-d^2 / (2 sigma_C^2))              competition kernel (C(0)=1)
//   single-resident equilibrium:  N* = K(x)
//   many-resident equilibrium:    solve A N = K, A_ij = C(x_i - x_j)
//   invasion fitness (a per-capita growth RATE, not a ratio; ~0 at the resident):
//       s(y) = r (1 - sum_i N_i C(y - x_i) / K(y))
//
// Analytic oracles: singular strategy x* = x0; a branching point (fitness
// minimum) iff sigma_C < sigma_K, an ESS (maximum) iff sigma_C > sigma_K.

#include <Rcpp.h>
#include <cmath>
#include <vector>
using namespace Rcpp;

namespace dd99 {

static inline double K_of(double x, double K0, double x0, double sK) {
  double d = x - x0;
  return K0 * std::exp(-(d * d) / (2.0 * sK * sK));
}
static inline double C_of(double d, double sC) {
  return std::exp(-(d * d) / (2.0 * sC * sC));
}

} // namespace dd99

//' DD99 model: invasion fitness of mutants (a per-capita growth rate)
//'
//' @param x_mut numeric vector of mutant trait values
//' @param x_res numeric vector of resident trait values
//' @param n_res numeric vector of resident equilibrium densities
//' @param pars list with r, K0, x0, sigma_K, sigma_C
//' @return numeric vector of invasion fitness (=0 for a resident at equilibrium)
//' @keywords internal
// [[Rcpp::export]]
NumericVector dd99_fitness(NumericVector x_mut, NumericVector x_res,
                           NumericVector n_res, List pars) {
  double r = pars["r"], K0 = pars["K0"], x0 = pars["x0"],
         sK = pars["sigma_K"], sC = pars["sigma_C"];
  int nm = x_mut.size();
  int nr = x_res.size();
  NumericVector out(nm);
  for (int i = 0; i < nm; i++) {
    double comp = 0.0;
    for (int j = 0; j < nr; j++)
      comp += n_res[j] * dd99::C_of(x_mut[i] - x_res[j], sC);
    double Ky = dd99::K_of(x_mut[i], K0, x0, sK);
    out[i] = r * (1.0 - comp / Ky);
  }
  return out;
}

//' DD99 model: resident demographic equilibrium densities
//'
//' Single resident: N* = K(x). Many residents: solve the linear system A N = K
//' with A_ij = C(x_i - x_j) (negative densities, i.e. strategies that cannot
//' coexist, are clamped to zero).
//'
//' @param x_res numeric vector of resident trait values
//' @param pars list with r, K0, x0, sigma_K, sigma_C
//' @return numeric vector of equilibrium densities
//' @keywords internal
// [[Rcpp::export]]
NumericVector dd99_equilibrium(NumericVector x_res, List pars) {
  double K0 = pars["K0"], x0 = pars["x0"],
         sK = pars["sigma_K"], sC = pars["sigma_C"];
  int nr = x_res.size();
  NumericVector n(nr);
  if (nr == 0) return n;
  if (nr == 1) {
    n[0] = dd99::K_of(x_res[0], K0, x0, sK);
    return n;
  }
  // dense linear solve via Gaussian elimination with partial pivoting
  std::vector<std::vector<double> > A(nr, std::vector<double>(nr + 1));
  for (int i = 0; i < nr; i++) {
    for (int j = 0; j < nr; j++) A[i][j] = dd99::C_of(x_res[i] - x_res[j], sC);
    A[i][nr] = dd99::K_of(x_res[i], K0, x0, sK);
  }
  for (int col = 0; col < nr; col++) {
    int piv = col;
    for (int i = col + 1; i < nr; i++)
      if (std::abs(A[i][col]) > std::abs(A[piv][col])) piv = i;
    std::swap(A[col], A[piv]);
    double d = A[col][col];
    for (int j = col; j <= nr; j++) A[col][j] /= d;
    for (int i = 0; i < nr; i++) {
      if (i == col) continue;
      double f = A[i][col];
      for (int j = col; j <= nr; j++) A[i][j] -= f * A[col][j];
    }
  }
  for (int i = 0; i < nr; i++) n[i] = std::max(0.0, A[i][nr]);
  return n;
}

// --- multi-trait (nD) DD99 (cf. Ito & Dieckmann 2007) ----------------------
// Traits and optima are k-dimensional; the resource and competition kernels are
// products of per-dimension Gaussians:
//   K(x) = K0 exp(-sum_d (x_d - x0_d)^2 / (2 sigma_K_d^2))
//   C(x',x) = exp(-sum_d (x'_d - x_d)^2 / (2 sigma_C_d^2))
// x_mut, x_res are matrices (rows = individuals, cols = traits); x0, sigma_K,
// sigma_C are length-k vectors. Same invasion fitness and equilibrium as the 1D
// case. Singular strategy x* = x0; branching in dimension d iff sigma_C_d <
// sigma_K_d.

namespace dd99 {

static inline double K_nd(const NumericMatrix& X, int i, double K0,
                          const NumericVector& x0, const NumericVector& sK) {
  int k = X.ncol();
  double e = 0.0;
  for (int d = 0; d < k; d++) {
    double z = X(i, d) - x0[d];
    e += z * z / (2.0 * sK[d] * sK[d]);
  }
  return K0 * std::exp(-e);
}

static inline double C_nd(const NumericMatrix& A, int i,
                          const NumericMatrix& B, int j,
                          const NumericVector& sC) {
  int k = A.ncol();
  double e = 0.0;
  for (int d = 0; d < k; d++) {
    double z = A(i, d) - B(j, d);
    e += z * z / (2.0 * sC[d] * sC[d]);
  }
  return std::exp(-e);
}

} // namespace dd99

//' DD99 model (nD): invasion fitness of mutants (a per-capita growth rate)
//'
//' @param x_mut numeric matrix of mutant trait values (rows = mutants, cols = traits)
//' @param x_res numeric matrix of resident trait values
//' @param n_res numeric vector of resident equilibrium densities
//' @param pars list with r, K0, x0, sigma_K, sigma_C (the last three length-k)
//' @return numeric vector of invasion fitness (=0 for a resident at equilibrium)
//' @keywords internal
// [[Rcpp::export]]
NumericVector dd99_nd_fitness(NumericMatrix x_mut, NumericMatrix x_res,
                              NumericVector n_res, List pars) {
  double r = pars["r"], K0 = pars["K0"];
  NumericVector x0 = pars["x0"], sK = pars["sigma_K"], sC = pars["sigma_C"];
  int nm = x_mut.nrow();
  int nr = x_res.nrow();
  NumericVector out(nm);
  for (int i = 0; i < nm; i++) {
    double Ky = dd99::K_nd(x_mut, i, K0, x0, sK);
    double comp = 0.0;
    for (int j = 0; j < nr; j++)
      comp += n_res[j] * dd99::C_nd(x_mut, i, x_res, j, sC);
    out[i] = r * (1.0 - comp / Ky);
  }
  return out;
}

//' DD99 model (nD): resident demographic equilibrium densities
//'
//' @param x_res numeric matrix of resident trait values
//' @param pars list with r, K0, x0, sigma_K, sigma_C (the last three length-k)
//' @return numeric vector of equilibrium densities
//' @keywords internal
// [[Rcpp::export]]
NumericVector dd99_nd_equilibrium(NumericMatrix x_res, List pars) {
  double K0 = pars["K0"];
  NumericVector x0 = pars["x0"], sK = pars["sigma_K"], sC = pars["sigma_C"];
  int nr = x_res.nrow();
  NumericVector n(nr);
  if (nr == 0) return n;
  if (nr == 1) { n[0] = dd99::K_nd(x_res, 0, K0, x0, sK); return n; }
  std::vector<std::vector<double> > A(nr, std::vector<double>(nr + 1));
  for (int i = 0; i < nr; i++) {
    for (int j = 0; j < nr; j++) A[i][j] = dd99::C_nd(x_res, i, x_res, j, sC);
    A[i][nr] = dd99::K_nd(x_res, i, K0, x0, sK);
  }
  for (int col = 0; col < nr; col++) {
    int piv = col;
    for (int i = col + 1; i < nr; i++)
      if (std::abs(A[i][col]) > std::abs(A[piv][col])) piv = i;
    std::swap(A[col], A[piv]);
    double d = A[col][col];
    for (int j = col; j <= nr; j++) A[col][j] /= d;
    for (int i = 0; i < nr; i++) {
      if (i == col) continue;
      double f = A[i][col];
      for (int j = col; j <= nr; j++) A[i][j] -= f * A[col][j];
    }
  }
  for (int i = 0; i < nr; i++) n[i] = std::max(0.0, A[i][nr]);
  return n;
}

// --- derivatives in the mutant direction -------------------------------------
// With E_j(y) = n_j C(y - x_j) / K(y) the fitness is s(y) = r (1 - sum_j E_j),
// and for product Gaussian kernels
//   dE_j/dy_d       = E_j a_jd,   a_jd = (y_d - x0_d)/sK_d^2 - (y_d - x_jd)/sC_d^2
//   d2E_j/dy_d dy_e = E_j (a_jd a_je + delta_de (1/sK_d^2 - 1/sC_d^2))
// so grad s = -r sum_j E_j a_j and hess s = -r sum_j E_j (a_j a_j^T + diag(1/sK^2 - 1/sC^2)).
// The 1D model is the k = 1 case and shares the code.

namespace dd99 {

static NumericMatrix nd_gradient(const NumericMatrix& x_mut, const NumericMatrix& x_res,
                                 const NumericVector& n_res, double r, double K0,
                                 const NumericVector& x0, const NumericVector& sK,
                                 const NumericVector& sC) {
  int nm = x_mut.nrow(), nr = x_res.nrow(), k = x_mut.ncol();
  NumericMatrix out(nm, k);
  for (int i = 0; i < nm; i++) {
    double Ky = K_nd(x_mut, i, K0, x0, sK);
    for (int j = 0; j < nr; j++) {
      double E = n_res[j] * C_nd(x_mut, i, x_res, j, sC) / Ky;
      for (int d = 0; d < k; d++) {
        double a = (x_mut(i, d) - x0[d]) / (sK[d] * sK[d])
                 - (x_mut(i, d) - x_res(j, d)) / (sC[d] * sC[d]);
        out(i, d) -= r * E * a;
      }
    }
  }
  return out;
}

static NumericMatrix nd_hessian(const NumericMatrix& x_mut, const NumericMatrix& x_res,
                                const NumericVector& n_res, double r, double K0,
                                const NumericVector& x0, const NumericVector& sK,
                                const NumericVector& sC) {
  int nr = x_res.nrow(), k = x_mut.ncol();
  NumericMatrix H(k, k);
  double Ky = K_nd(x_mut, 0, K0, x0, sK);
  std::vector<double> a(k);
  for (int j = 0; j < nr; j++) {
    double E = n_res[j] * C_nd(x_mut, 0, x_res, j, sC) / Ky;
    for (int d = 0; d < k; d++)
      a[d] = (x_mut(0, d) - x0[d]) / (sK[d] * sK[d])
           - (x_mut(0, d) - x_res(j, d)) / (sC[d] * sC[d]);
    for (int d = 0; d < k; d++)
      for (int e = 0; e < k; e++) {
        double curv = (d == e) ? 1.0 / (sK[d] * sK[d]) - 1.0 / (sC[d] * sC[d]) : 0.0;
        H(d, e) -= r * E * (a[d] * a[e] + curv);
      }
  }
  return H;
}

static NumericMatrix as_column(const NumericVector& v) {
  NumericMatrix m(v.size(), 1);
  for (int i = 0; i < v.size(); i++) m(i, 0) = v[i];
  return m;
}

} // namespace dd99

//' DD99 model: gradient of invasion fitness with respect to the mutant trait
//'
//' @inheritParams dd99_fitness
//' @return numeric matrix, one row per mutant and one column (the trait)
//' @keywords internal
// [[Rcpp::export]]
NumericMatrix dd99_fitness_gradient(NumericVector x_mut, NumericVector x_res,
                                    NumericVector n_res, List pars) {
  double r = pars["r"], K0 = pars["K0"];
  NumericVector x0 = NumericVector::create(pars["x0"]),
                sK = NumericVector::create(pars["sigma_K"]),
                sC = NumericVector::create(pars["sigma_C"]);
  return dd99::nd_gradient(dd99::as_column(x_mut), dd99::as_column(x_res),
                           n_res, r, K0, x0, sK, sC);
}

//' DD99 model: Hessian of invasion fitness with respect to the mutant trait
//'
//' @param x_mut a single mutant trait value
//' @inheritParams dd99_fitness
//' @return a 1 x 1 numeric matrix
//' @keywords internal
// [[Rcpp::export]]
NumericMatrix dd99_fitness_hessian(NumericVector x_mut, NumericVector x_res,
                                   NumericVector n_res, List pars) {
  if (x_mut.size() != 1) stop("dd99_fitness_hessian takes a single mutant");
  double r = pars["r"], K0 = pars["K0"];
  NumericVector x0 = NumericVector::create(pars["x0"]),
                sK = NumericVector::create(pars["sigma_K"]),
                sC = NumericVector::create(pars["sigma_C"]);
  return dd99::nd_hessian(dd99::as_column(x_mut), dd99::as_column(x_res),
                          n_res, r, K0, x0, sK, sC);
}

//' DD99 model (nD): gradient of invasion fitness with respect to the mutant traits
//'
//' @inheritParams dd99_nd_fitness
//' @return numeric matrix, one row per mutant and one column per trait
//' @keywords internal
// [[Rcpp::export]]
NumericMatrix dd99_nd_fitness_gradient(NumericMatrix x_mut, NumericMatrix x_res,
                                       NumericVector n_res, List pars) {
  double r = pars["r"], K0 = pars["K0"];
  NumericVector x0 = pars["x0"], sK = pars["sigma_K"], sC = pars["sigma_C"];
  return dd99::nd_gradient(x_mut, x_res, n_res, r, K0, x0, sK, sC);
}

//' DD99 model (nD): Hessian of invasion fitness with respect to the mutant traits
//'
//' @param x_mut a one-row numeric matrix: the mutant trait values
//' @inheritParams dd99_nd_fitness
//' @return a k x k numeric matrix
//' @keywords internal
// [[Rcpp::export]]
NumericMatrix dd99_nd_fitness_hessian(NumericMatrix x_mut, NumericMatrix x_res,
                                      NumericVector n_res, List pars) {
  if (x_mut.nrow() != 1) stop("dd99_nd_fitness_hessian takes a single mutant");
  double r = pars["r"], K0 = pars["K0"];
  NumericVector x0 = pars["x0"], sK = pars["sigma_K"], sC = pars["sigma_C"];
  return dd99::nd_hessian(x_mut, x_res, n_res, r, K0, x0, sK, sC);
}
