#' EXPOKIT dmexpv wrapper function
#'
#' This function wraps the .C call to EXPOKIT for the dmexpv function.  Only the output probability
#' matrix is returned.
#'
#' @param n number of rows in Q matrix
#' @param m n-1
#' @param t the value to exponentiate the rate matrix by (often e.g. a time value)
#' @param v variable to store some results in; should have n elements (and perhaps start with 1)
#' @param w same length as v
#' @param tol tolerance for approximations; usually set to 0.01
#' @param anorm the norm of the Q matrix
#' @param lwsp length of workspace (wsp); for dmexpv, lwsp=n*(m+2)+5*(m+2)^2+ideg+1
#' @param wsp workspace to store some results in; should be a double with lwsp elements
#' @param liwsp length of integer workspace; for dmexpv, liwsp=m+2
#' @param iwsp integer workspace
#' @param itrace option, set to 0
#' @param iflag option, set to 0
#' @param ia i indices of Qmat nonzero values
#' @param ja j indices of Qmat nonzero values
#' @param a nonzero values of Qmat (ia, ja, a are columns of a COO-formatted Q matrix)
#' @param nz number of non-zeros in Qmat
#' @param res space for output probability matrix (n x n)
#'
#' EXPOKIT needs the input matrix to be transposed compared to normal.
#' COO format is required for EXPOKIT.
#' @return \code{tmpoutmat} the output matrix for the (first) input t-value
#' @seealso \code{\link{expokit_dmexpv_Qmat}}
#' @seealso \code{\link{expokit_wrapalldmexpv_tvals}}
#' @export
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu}
#' @examples
#' # Make a square instantaneous rate matrix (Q matrix)
#' # This matrix is taken from Peter Foster's (2001) "The Idiot's Guide
#' # to the Zen of Likelihood in a Nutshell in Seven Days for Dummies,
#' # Unleashed" at:
#' # \url{http://www.bioinf.org/molsys/data/idiots.pdf}
#' #
#' # The Q matrix includes the stationary base freqencies, which Pmat 
#' # converges to as t becomes large.
#' Qmat = matrix(c(-1.218, 0.504, 0.336, 0.378, 0.126, -0.882, 0.252, 0.504, 0.168, 
#' 0.504, -1.05, 0.378, 0.126, 0.672, 0.252, -1.05), nrow=4, byrow=TRUE)
#' 
#' # Make a series of t values
#' tvals = c(0.001, 0.005, 0.01, 0.05, 0.1, 0.5, 1, 2, 5, 14)
#' 
#' # DMEXPV and DGEXPV are designed for large, sparse Q matrices (sparse = lots of zeros).
#' # DMEXPV is specifically designed for Markov chains and so may be slower, but more accurate.
#' 
#' # DMEXPV, single t-value
#' expokit_wrapalldmexpv_tvals(Qmat=Qmat, tvals=tvals[1], transpose_needed=TRUE)
#' expokit_wrapalldmexpv_tvals(Qmat=Qmat, tvals=2)
#' 
#' # This function runs a for-loop itself (sadly, we could not get mapply() to work
#' # on a function that calls dmexpv/dgexpv), returning a list of probability matrices.
#' 
#' # DMEXPV functions
#' list_of_P_matrices_dmexpv = expokit_wrapalldmexpv_tvals(Qmat=Qmat, 
#' tvals=tvals, transpose_needed=TRUE)
#' list_of_P_matrices_dmexpv
#' 
expokit_dmexpv_wrapper <- function(n, m, t, v, tol, anorm, wsp, lwsp, iwsp, liwsp, itrace, iflag, ia, ja, a, nz)
{
	ret <- .Call("R_dmexpv", 
		           as.integer(n), as.integer(m), as.double(t), 
		           as.double(v), as.double(tol), 
		           as.double(anorm), as.double(wsp), as.integer(lwsp), 
		           as.integer(iwsp), as.integer(liwsp), 
		           as.integer(ia), as.integer(ja), 
		           as.double(a), as.integer(nz),
		           PACKAGE="rexpokit")
	output_Pmat = matrix(ret$res, nrow=n, byrow=TRUE)
	
	return(output_Pmat)
}






#' EXPOKIT dmexpv wrapper function, return just output probs
#'
#' This function wraps the .C call to EXPOKIT for the dmexpv function.  Only the output probabilities
#' not the Pmat probability matrix, are returned.
#'
#' @param n number of rows in Q matrix
#' @param m n-1
#' @param t the value to exponentiate the rate matrix by (often e.g. a time value)
#' @param v variable to store some results in; should have n elements (and perhaps start with 1)
#' @param w same length as v
#' @param tol tolerance for approximations; usually set to 0.01
#' @param anorm the norm of the Q matrix
#' @param lwsp length of workspace (wsp); for dmexpv, lwsp=n*(m+2)+5*(m+2)^2+ideg+1
#' @param wsp workspace to store some results in; should be a double with lwsp elements
#' @param liwsp length of integer workspace; for dmexpv, liwsp=m+2
#' @param iwsp integer workspace
#' @param itrace option, set to 0
#' @param iflag option, set to 0
#' @param ia i indices of Qmat nonzero values
#' @param ja j indices of Qmat nonzero values
#' @param a nonzero values of Qmat (ia, ja, a are columns of a COO-formatted Q matrix)
#' @param nz number of non-zeros in Qmat
#'
#' EXPOKIT needs the input matrix to be transposed compared to normal.
#' COO format is required for EXPOKIT.
#' @return \code{w_output_probs} the output probabilities (= \code{myDMEXPV} variable \code{w}, or the fifth output
#' in the output from .Call("mydmexpv_", ...), given the (first) input t-value.
#' @seealso \code{\link{expokit_dmexpv_Qmat}}
#' @seealso \code{\link{expokit_wrapalldmexpv_tvals}}
#' @export
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu}
#' @examples
#' # Make a square instantaneous rate matrix (Q matrix)
#' # This matrix is taken from Peter Foster's (2001) "The Idiot's Guide
#' # to the Zen of Likelihood in a Nutshell in Seven Days for Dummies,
#' # Unleashed" at:
#' # \url{http://www.bioinf.org/molsys/data/idiots.pdf}
#' #
#' # The Q matrix includes the stationary base freqencies, which Pmat 
#' # converges to as t becomes large.
#' Qmat = matrix(c(-1.218, 0.504, 0.336, 0.378, 0.126, -0.882, 0.252, 0.504, 0.168, 
#' 0.504, -1.05, 0.378, 0.126, 0.672, 0.252, -1.05), nrow=4, byrow=TRUE)
#' 
#' # Make a series of t values
#' tvals = c(0.001, 0.005, 0.01, 0.05, 0.1, 0.5, 1, 2, 5, 14)
#' 
#' # DMEXPV and DGEXPV are designed for large, sparse Q matrices (sparse = lots of zeros).
#' # DMEXPV is specifically designed for Markov chains and so may be slower, but more accurate.
#' 
#' # DMEXPV, single t-value
#' expokit_wrapalldmexpv_tvals(Qmat=Qmat, tvals=tvals[1], transpose_needed=TRUE)
#' expokit_wrapalldmexpv_tvals(Qmat=Qmat, tvals=2)
#' 
#' # This function runs a for-loop itself (sadly, we could not get mapply() to work
#' # on a function that calls dmexpv/dgexpv), returning a list of probability matrices.
#' 
#' # DMEXPV functions
#' list_of_P_matrices_dmexpv = expokit_wrapalldmexpv_tvals(Qmat=Qmat, 
#' tvals=tvals, transpose_needed=TRUE)
#' list_of_P_matrices_dmexpv
#' 
expokit_mydmexpv_wrapper <- function(n, m, t, v, tol, anorm, wsp, lwsp, iwsp, liwsp, itrace, iflag, ia, ja, a, nz)
{
	ret <- .Call("R_dmexpv", 
		           as.integer(n), as.integer(m), as.double(t), 
		           as.double(v), as.double(tol), 
		           as.double(anorm), as.double(wsp), as.integer(lwsp), 
		           as.integer(iwsp), as.integer(liwsp), 
		           as.integer(ia), as.integer(ja), 
		           as.double(a), as.integer(nz),
		           PACKAGE="rexpokit")
	
	w_output_probs = matrix(ret$w, ncol=n, byrow=TRUE)
	
	return(w_output_probs)
}





#' EXPOKIT dgexpv wrapper function, return just output probs
#'
#' This function wraps the .C call to EXPOKIT for the dgexpv function.  Only the output probabilities
#' not the Pmat probability matrix, are returned.
#'
#' @param n number of rows in Q matrix
#' @param m n-1
#' @param t the value to exponentiate the rate matrix by (often e.g. a time value)
#' @param v variable to store some results in; should have n elements (and perhaps start with 1)
#' @param w same length as v
#' @param tol tolerance for approximations; usually set to 0.01
#' @param anorm the norm of the Q matrix
#' @param lwsp length of workspace (wsp); for dgexpv, lwsp=n*(m+2)+5*(m+2)^2+ideg+1
#' @param wsp workspace to store some results in; should be a double with lwsp elements
#' @param liwsp length of integer workspace; for dgexpv, liwsp=m+2
#' @param iwsp integer workspace
#' @param itrace option, set to 0
#' @param iflag option, set to 0
#' @param ia i indices of Qmat nonzero values
#' @param ja j indices of Qmat nonzero values
#' @param a nonzero values of Qmat (ia, ja, a are columns of a COO-formatted Q matrix)
#' @param nz number of non-zeros in Qmat
#'
#' EXPOKIT needs the input matrix to be transposed compared to normal.
#' COO format is required for EXPOKIT.
#' @return \code{w_output_probs} the output probabilities (= \code{myDGEXPV} variable \code{w}, or the fifth output
#' in the output from .Call("mydgexpv_", ...), given the (first) input t-value.
#' @seealso \code{\link{expokit_dgexpv_Qmat}}
#' @seealso \code{\link{expokit_wrapalldgexpv_tvals}}
#' @export
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu}
#' @examples
#' # Make a square instantaneous rate matrix (Q matrix)
#' # This matrix is taken from Peter Foster's (2001) "The Idiot's Guide
#' # to the Zen of Likelihood in a Nutshell in Seven Days for Dummies,
#' # Unleashed" at:
#' # \url{http://www.bioinf.org/molsys/data/idiots.pdf}
#' #
#' # The Q matrix includes the stationary base freqencies, which Pmat 
#' # converges to as t becomes large.
#' Qmat = matrix(c(-1.218, 0.504, 0.336, 0.378, 0.126, -0.882, 0.252, 0.504, 0.168, 
#' 0.504, -1.05, 0.378, 0.126, 0.672, 0.252, -1.05), nrow=4, byrow=TRUE)
#' 
#' # Make a series of t values
#' tvals = c(0.001, 0.005, 0.01, 0.05, 0.1, 0.5, 1, 2, 5, 14)
#' 
#' # dgexpv and DGEXPV are designed for large, sparse Q matrices (sparse = lots of zeros).
#' # dgexpv is specifically designed for Markov chains and so may be slower, but more accurate.
#' 
#' # dgexpv, single t-value
#' expokit_wrapalldgexpv_tvals(Qmat=Qmat, tvals=tvals[1], transpose_needed=TRUE)
#' expokit_wrapalldgexpv_tvals(Qmat=Qmat, tvals=2)
#' 
#' # This function runs a for-loop itself (sadly, we could not get mapply() to work
#' # on a function that calls dgexpv/dgexpv), returning a list of probability matrices.
#' 
#' # dgexpv functions
#' list_of_P_matrices_dgexpv = expokit_wrapalldgexpv_tvals(Qmat=Qmat, 
#' tvals=tvals, transpose_needed=TRUE)
#' list_of_P_matrices_dgexpv
#' 
expokit_mydgexpv_wrapper <- function(n, m, t, v, w, tol, anorm, wsp, lwsp, iwsp, liwsp, itrace, iflag, ia, ja, a, nz)
{
	res2 = NULL
	
	# This must be mydgexpv_, not mydgexpv_ !!!!
	
	res2 <- .C("mydgexpv_", as.integer(n), as.integer(m), as.double(t), as.double(v), as.double(w), as.double(tol), as.double(anorm), as.double(wsp), as.integer(lwsp), as.integer(iwsp), as.integer(liwsp), as.integer(itrace), as.integer(iflag), as.integer(ia), as.integer(ja), as.double(a), as.integer(nz))
	
	w_output_probs = matrix(res2[[5]], ncol=n, byrow=TRUE)
	
	return(w_output_probs)
}


