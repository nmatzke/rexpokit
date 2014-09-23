#include <Rcpp.h>

// Sparsity of densely stored matrix
int sparse_count_zeros(int m, int n, double *x)
{
  int count = 0;
  int i, j;
  
  for (j=0; j<n; j++)
  {
    for (i=0; i<m; i++)
    {
      if (x[i + m*j] == 0.)
        count++;
    }
  }
  
  return count;
}

// Reimplentation of SparseM's dense-to-coo converter
extern "C" {
SEXP rexpokit_as_coo(SEXP x_)
{
  Rcpp::NumericMatrix x(x_);
  
  int ct = 0, i, j;
  const int n = x.nrow();
  const int sparsity = sparse_count_zeros(n, n, x.begin());
  const int len = n*n - sparsity;
  
  Rcpp::NumericMatrix ret(len, 3);
  
  for (i=0; i<n; i++)
  {
    for (j=0; j<n; j++)
    {
      if (x(i, j) != 0.)
      {
        ret(0, ct) = i + 1;       // ia
        ret(1, ct) = j + 1;       // ja
        ret(2, ct) = x(i, j);     // a
        
        ct++;
      }
    }
  }
  
  
  return ret;
}
}
