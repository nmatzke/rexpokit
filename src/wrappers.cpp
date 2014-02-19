#include <Rcpp.h>

extern "C"
{

void dgpadm_(int *ideg, int *m, double *t, double *H, int *ldh, double *wsp, int *lwsp, double *ipiv, int *iexph, int *ns, int *iflag);

SEXP R_dgpadm(SEXP ideg, SEXP m_, SEXP t, SEXP H, SEXP ldh)
{
  int m = INTEGER(m_)[0];
  int ns = 0, iflag = 0;
  int lwsp = 4*m*m + INTEGER(ideg)[0] + 1;
  
  Rcpp::NumericVector wsp(lwsp);
  Rcpp::NumericVector ipiv(m);
  
  Rcpp::IntegerVector iexph(1);
  
  Rcpp::List ret; 
  
  dgpadm_(INTEGER(ideg), &m, REAL(t), REAL(H), INTEGER(ldh), REAL(wsp), &lwsp, REAL(ipiv), INTEGER(iexph), &ns, &iflag);
  
  ret["wsp"] = wsp; 
  ret["ind"] = iexph;
  
  return ret;
}

}
