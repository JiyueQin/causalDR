#include <Rcpp.h>
#include <algorithm>
#include <cmath>
#include <vector>
using namespace Rcpp;

// Augmentation quantities for Lambda_3 under left truncation (AIPW-trunc).
// For each subject i and each target time index k in k_target (0-based grid indices):
//   R_diagonal(i, k)   = R(t_k, t_k | a, Z_i)   = P(T1 < t_k <= T2 | a, Z_i)
//   rho_diagonal(i, k) = rho(t_k, t_k | a, Z_i) = P(T1 < t_k, T2 = t_k | a, Z_i)
//   m3(i, k), m4(i, k): truncation augmentation terms (only computed when need_m[i]=True)
// where, under the conditional model of T2 given T1 = u,
//   R(v, t) = sum_{u <= v, u < t} b(u, t),
//   b(u, t) = S(u- | a, Z) dLambda1(u | a, Z) exp{-Lambda3(t- | u, a, Z)},
//   Lambda3(t- | u, a, Z) = {B3(t-) - B3(u)} exp{lp3(u, a, Z)},
//   dLambda3(t | u, a, Z) = dB3(t) exp{lp3(u, a, Z)},
// and rho(v, t) is defined in the same way with the extra factor dLambda3(t | u).
// With g(v) = dG(v)/{S(v) G(v)^2} and P(x) = sum_{v < x} g(v) (C(x) in
// dLambda_trunc(); S is for arm a, which equals S at the observed treatment for the
// subjects with need_m),
//   m3 = R(Q, t)/{S(Q) G(Q)} - sum_{v in [Q, min(X1, t))} R(v, t) g(v)
//      = sum_{u <= Q} b(u, t)/{S(Q) G(Q)}
//        - sum_{u < min(X1, t)} b(u, t) [P{min(X1, t)} - P{max(u, Q)}],
// so that only the jump points of Lambda1 need to be visited.
//
// Called from dLambda_trunc() in R/trunc_estimation.R. It is written in C++ because
// R(t, t) and rho(t, t) need a sum over the jump points u of Lambda1 for every
// subject and every target time, a triple loop that does not vectorize well in R.
//
// Grid indices are 0-based, as in C++: index j refers to grid[j + 1] in R, so the
// caller passes R's 1-based indices minus 1 (e.g., k3 - 1L). The returned matrices
// need no conversion, since their columns follow the order of k_target.
//
// Arguments (n x m matrices are subjects x grid points):
//   S         NumericMatrix n x m, S(v | a, Z) for the arm a of interest
//   dLambda1  NumericMatrix n x m, the jumps of Lambda1(v | a, Z) for the arm a
//   lp3       NumericMatrix n x m, lp3(i, j) = lp3(u_j, a, Z_i), the linear
//             predictor of Lambda3 for subject i given T1 = u_j, the j-th grid point,
//             for the arm a
//   B3        NumericMatrix n x m, baseline cumulative hazard of Lambda3 on the grid
//   dB3       NumericMatrix n x m, the jumps of B3
//   G         NumericMatrix n x m, G(v | A, Z) at the observed treatment
//   dG        NumericMatrix n x m, the jumps of G
//   index_Q   IntegerVector of length n, 0-based grid indices of Q
//   index_X1  IntegerVector of length n, 0-based grid indices of X1
//   need_m    LogicalVector of length n, whether m3 and m4 are needed, i.e.,
//             1(A = a, Delta = 1); they enter the estimating equation only through
//             V3 = 1(A = a) Delta3/{pi_a S_D}, which is 0 otherwise, so they are
//             skipped (left at 0) for the other subjects to save time
//   k_target  IntegerVector of length K, 0-based grid indices at which the
//             quantities are computed: the target times t_k at which the jump of
//             Lambda3 is estimated, i.e., observed Lambda3 event times or jump times
//             of B3, up to the largest time of interest; elsewhere the jump is 0, as
//             dN3 = 0 and dB3 = 0 there, so that rho and m4 are 0
//   S_min     double in (0, 1), lower bound of S in denominators
//   G_min     double in (0, 1), lower bound of G in denominators
//
// Returns a List of four n x K NumericMatrix: R_diagonal, rho_diagonal, m3 and m4,
// with columns in the order of k_target; m3 and m4 are 0 where need_m is FALSE or Q >= t_k.
// [[Rcpp::export]]
List trunc_aug_cpp(NumericMatrix S, NumericMatrix dLambda1, NumericMatrix lp3,
                   NumericMatrix B3, NumericMatrix dB3,
                   NumericMatrix G, NumericMatrix dG,
                   IntegerVector index_Q, IntegerVector index_X1, LogicalVector need_m,
                   IntegerVector k_target, double S_min, double G_min) {
  int n = S.nrow();
  int m = S.ncol();
  int K = k_target.size();
  NumericMatrix R_diagonal(n, K), rho_diagonal(n, K), m3(n, K), m4(n, K);
  // P[j] = sum_{v < j} g(v), cumulative sum of g on the grid
  std::vector<double> P(m + 1);
  // jump points u of Lambda1, with S(u-) dLambda1(u) and exp{lp3(u)}
  std::vector<int> jump_idx;
  std::vector<double> jump_mass, exp_lp3;

  for (int i = 0; i < n; ++i) {
    // Step 1: b(u, t) is nonzero only at the jump points u of Lambda1, so collect
    // them once per subject, with the parts of b(u, t) that do not depend on t:
    // P(T1 = u | a, Z) = S(u-) dLambda1(u) and exp{lp3(u, a, Z)}. The jump points
    // are stored in increasing order, which Step 3 uses to stop early.
    jump_idx.clear(); jump_mass.clear(); exp_lp3.clear();
    for (int j = 0; j < m; ++j) {
      if (dLambda1(i, j) > 0) {
        double S_minus = (j > 0) ? S(i, j - 1) : 1.0;
        jump_idx.push_back(j);
        jump_mass.push_back(S_minus * dLambda1(i, j));
        exp_lp3.push_back(std::exp(lp3(i, j)));
      }
    }
    int n_jump = jump_idx.size();

    // Step 2 (only if m3 and m4 are needed): the cumulative sum P of
    // g(v) = dG(v)/{S(v) G(v)^2}, with P[j] = sum over grid indices v < j, so that
    // any partial sum of g over [v1, v2) is P[v2] - P[v1]; and S(Q) G(Q). S and G
    // are bounded below by S_min and G_min, as in dLambda_trunc().
    double SG_Q = 1.0;
    if (need_m[i]) {
      P[0] = 0.0;
      for (int j = 0; j < m; ++j) {
        double S_j = std::max(S(i, j), S_min);
        double G_j = std::max(G(i, j), G_min);
        P[j + 1] = P[j] + ((dG(i, j) > 0) ? dG(i, j) / (S_j * G_j * G_j) : 0.0);
      }
      SG_Q = std::max(S(i, index_Q[i]), S_min) * std::max(G(i, index_Q[i]), G_min);
    }

    // Step 3: for each target time t_k, accumulate R(t_k, t_k), rho(t_k, t_k), m3 and
    // m4 in one pass over the jump points u < t_k.
    for (int kk = 0; kk < K; ++kk) {
      int k = k_target[kk];
      double B3_minus = (k > 0) ? B3(i, k - 1) : 0.0;   // B3(t_k-)
      double dB3_k = dB3(i, k);
      // the truncation terms are 0 unless Q < t_k
      bool add_m = need_m[i] && (index_Q[i] < k);
      int j_end = std::min(index_X1[i], k);                  // grid index of min(X1, t_k)
      double R = 0.0, rho = 0.0, m3_ik = 0.0, m4_ik = 0.0;

      for (int l = 0; l < n_jump; ++l) {
        int u = jump_idx[l];
        // R(t, t) and rho(t, t) sum over u < t; the jump points are increasing
        if (u >= k) break;
        // b(u, t_k) = P(T1 = u, T2 >= t_k | a, Z), using
        // P(T2 >= t_k | T1 = u) = exp[-{B3(t_k-) - B3(u)} exp{lp3(u)}]
        double b = jump_mass[l] * std::exp(-(B3_minus - B3(i, u)) * exp_lp3[l]);
        // times dLambda3(t_k | u) = dB3(t_k) exp{lp3(u)}: P(T1 = u, T2 = t_k | a, Z)
        double b_rho = b * exp_lp3[l] * dB3_k;
        R += b;
        rho += b_rho;
        if (add_m) {
          // first term of m3, m4: R(Q, t)/{S(Q) G(Q)} collects u <= Q
          if (u <= index_Q[i]) {
            m3_ik += b / SG_Q;
            m4_ik += b_rho / SG_Q;
          }
          // second term: exchanging the sums over v and u, each u < min(X1, t) is
          // weighted by the sum of g(v) over v in [max(u, Q), min(X1, t))
          if (u < j_end) {
            double P_diff = P[j_end] - P[std::max(u, index_Q[i])];
            m3_ik -= b * P_diff;
            m4_ik -= b_rho * P_diff;
          }
        }
      }
      R_diagonal(i, kk) = R;
      rho_diagonal(i, kk) = rho;
      m3(i, kk) = m3_ik;
      m4(i, kk) = m4_ik;
    }
  }
  return List::create(_["R_diagonal"] = R_diagonal, _["rho_diagonal"] = rho_diagonal,
                      _["m3"] = m3, _["m4"] = m4);
}
