#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix get_integral_cpp (NumericVector t,
                                NumericMatrix X,
                                NumericMatrix Y) {
  int n = X.nrow();
  int m = X.ncol();
  NumericMatrix Z(n, m);
  
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < m; ++j) {
      double x_ij = X(i, j);
      double sum = 0.0;
      for (int k = 0; k < m; ++k) {
        if (t[k] <= x_ij) {
          sum += Y(i, k);
        }
      }
      Z(i, j) = sum;
    }
  }
  return Z;
}
