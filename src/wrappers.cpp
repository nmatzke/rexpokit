#include <Rcpp.h>
#include "expokit.h"

extern "C"
{

// Computes matrix exponential for dense input H
SEXP R_dgpadm(SEXP ideg, SEXP m_, SEXP t, SEXP H, SEXP ldh)
{
  int m = INTEGER(m_)[0];
  int ns = 0, iflag = 0;
  int lwsp = 4*m*m + INTEGER(ideg)[0] + 1;
  
  Rcpp::NumericVector wsp(lwsp);
  Rcpp::IntegerVector ipiv(m);
  
  Rcpp::IntegerVector iexph(1);
  
  Rcpp::List ret; 
  
  dgpadm_(INTEGER(ideg), &m, REAL(t), REAL(H), INTEGER(ldh), REAL(wsp), &lwsp, INTEGER(ipiv), INTEGER(iexph), &ns, &iflag);
  
  ret["wsp"] = wsp; 
  ret["ind"] = iexph;
  
  return ret;
}

}
