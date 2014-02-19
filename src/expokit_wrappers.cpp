extern"C" {
// These functions are defined in my_matexp.f, which calls the following functions/

// DMEXPV contains an additional check ensuring sums to 1, for Markov-chain applications

// This exponentiates a large, sparse matrix (sparse = lots of zeros)
// Before input, the matrix should be transposed and
// put into coordinate list (COO) format: 
// http://en.wikipedia.org/wiki/Sparse_matrix#Coordinate_list_.28COO.29
//
// i.e.:
// ia = row number
// ja = col number
// a  = corresponding value of that cell (zeros are excluded)
void wrapalldmexpv_(int * n,int* m,double * t,double* v,double * w,double* tol,
	double* anorm,double* wsp,int * lwsp,int* iwsp,int *liwsp, int * itrace,int *iflag,
		int *ia, int *ja, double *a, int *nz, double * res);

// This returns just one row (?) of the transition matrix, useful in Lagrange;
// same inputs as wrapalldmexpv_.
void wrapsingledmexpv_(int * n,int* m,double * t,double* v,double * w,double* tol,
	double* anorm,double* wsp,int * lwsp,int* iwsp,int *liwsp, int * itrace,int *iflag,
			int *ia, int *ja, double *a, int *nz, double * res);

// DGEXPV contains an additional check ensuring sums to 1, for Markov-chain applications
void wrapalldgexpv_(int * n,int* m,double * t,double* v,double * w,double* tol,
	double* anorm,double* wsp,int * lwsp,int* iwsp,int *liwsp, int * itrace,int *iflag,
		int *ia, int *ja, double *a, int *nz, double * res);

// This returns just one row (?) of the transition matrix, useful in Lagrange;
// same inputs as wrapalldmexpv_.
void wrapsingledgexpv_(int * n,int* m,double * t,double* v,double * w,double* tol,
	double* anorm,double* wsp,int * lwsp,int* iwsp,int *liwsp, int * itrace,int *iflag,
			int *ia, int *ja, double *a, int *nz, double * res);


// The myDMEXPV etc. functions provide direct access to EXPOKIT functions;
// This should be faster, especially for sparse matrices
// Here, you input v (starting probabilities) and it fills in w, which are the
// output probabilities (in output list item #5)
void myDMEXPV_(int* n, int* m, double* t, double* v, double* w, double* tol,
	double* anorm, double* wsp, int* lwsp, int* iwsp, int* liwsp, int* itrace, int* iflag,
		int* ia, int* ja, double* a, int* nz );

void myDGEXPV_(int* n, int* m, double* t, double* v, double* w, double* tol,
	double* anorm, double* wsp, int* lwsp, int* iwsp, int* liwsp, int* itrace, int* iflag,
		int* ia, int* ja, double* a, int* nz );


// This exponentiates a small, dense (no zeros) matrix with the padm
// approximation.  Here, the input matrix is just the 
// matrix, tranposed and then put into a list of numbers, e.g.:
// 
/*
# R code
# transpose the matrix in R
tmat = t(mat)
# convert into list of numbers (as numeric or double)
tmat_nums = as.double(tmat)
*/
//void wrapdgpadm_(int * ideg,int * m,double * t,double * H,int * ldh,
//	double * wsp,int * lwsp,int * ipiv,int * iexph,int *ns,int *iflag );

} // end extern "C"

