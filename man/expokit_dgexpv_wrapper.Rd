\name{expokit_dgexpv_wrapper}
\alias{expokit_dgexpv_wrapper}
\title{EXPOKIT dgexpv wrapper function}
\usage{
  expokit_dgexpv_wrapper(n, m, timeval, v, w, tol, anorm,
    wsp, lwsp, iwsp, liwsp, itrace, iflag, ia, ja, a, nz,
    res)
}
\arguments{
  \item{n}{number of rows in Q matrix}

  \item{m}{n-1}

  \item{timeval}{the value to exponentiate the rate matrix
  by (often e.g. a time value)}

  \item{v}{variable to store some results in; should have n
  elements (and perhaps start with 1)}

  \item{w}{same length as v}

  \item{tol}{tolerance for approximations; usually set to
  0.01}

  \item{anorm}{the norm of the Q matrix}

  \item{lwsp}{length of workspace (wsp); for dgexpv,
  lwsp=n*(m+2)+5*(m+2)^2+ideg+1}

  \item{wsp}{workspace to store some results in; should be
  a double with lwsp elements}

  \item{liwsp}{length of integer workspace; for dgexpv,
  liwsp=m+2}

  \item{iwsp}{integer workspace}

  \item{itrace}{option, set to 0}

  \item{iflag}{option, set to 0}

  \item{ia}{i indices of Qmat nonzero values}

  \item{ja}{j indices of Qmat nonzero values}

  \item{a}{nonzero values of Qmat (ia, ja, a are columns of
  a COO-formatted Q matrix)}

  \item{nz}{number of non-zeros in Qmat}

  \item{res}{space for output probability matrix (n x n)

  EXPOKIT needs the input matrix to be transposed compared
  to normal) COO format is required for EXPOKIT.}
}
\value{
  \code{tmpoutmat} the output matrix for the (first) input
  t-value
}
\description{
  This function wraps the .C call to EXPOKIT for the dgexpv
  function.  Only the output probability matrix is
  returned.
}
\details{
  NOTE: DGEXPV vs. DMEXPV. According to the EXPOKIT
  documentation, DGEXPV should be faster than DMEXPV,
  however DMEXPV runs an accuracy check appropriate for
  Markov chains, which is not done in DGEXPV.
}
\examples{
# Example building the inputs from scratch:

# Make a square instantaneous rate matrix (Q matrix)
# This matrix is taken from Peter Foster's (2001) "The Idiot's Guide
# to the Zen of Likelihood in a Nutshell in Seven Days for Dummies,
# Unleashed" at:
# \\url{http://www.bioinf.org/molsys/data/idiots.pdf}
#
# The Q matrix includes the stationary base freqencies, which Pmat
# converges to as t becomes large.
Qmat = matrix(c(-1.218, 0.504, 0.336, 0.378, 0.126, -0.882, 0.252, 0.504, 0.168,
0.504, -1.05, 0.378, 0.126, 0.672, 0.252, -1.05), nrow=4, byrow=TRUE)

# Make a series of t values
tvals = c(0.001, 0.005, 0.01, 0.05, 0.1, 0.5, 1, 2, 5, 14)
timeval = tvals[2]

	ideg = as.integer(6)
	n=nrow(Qmat)
	m=n-1
	# t=as.numeric(2.1)

	# v should have as many elements as n; first element = 1 (?)
	v=double(n)
	v[1] = 1

	# w is the same length
	w = double(length=n)
	tol=as.numeric(0.01)

	# length of wsp
	#lwsp = as.integer(n*(m+1)+n+(m+2)^2 + 4*(m+2)^2+ideg+1)
	#lwsp = as.integer(n*(m+1)+n+(m+2)^2 + 5*(m+2)^2+ideg+1)
	lwsp = as.integer(n*(m+2)+5*(m+2)^2+ideg+1)

	#lwsp = 100
	wsp = double(length=lwsp)

	# length of iwsp
	liwsp = max(m+2, 7)
	iwsp = integer(length=liwsp)

	res = double(length=n*n)

	#matvec = matrix(data=Q, nrow=n, byrow=TRUE)
	matvec = Qmat
	tmatvec = t(matvec)
	rowSums(tmatvec)
	colSums(tmatvec)

	# type="O" is being used here, this is supposed to be the
	# default for norm(), although it throws an error if not
	# specified
	#
	# From the help:
	# type - character string, specifying the type of matrix norm to be
	# computed. A character indicating the type of norm desired.
	# 	"O", "o" or "1"
	# 		specifies the one norm, (maximum absolute column sum);
	anorm = as.numeric(norm(matvec, type="O"))
	#anorm = 1


	itrace = 0
	iflag = 0


	#a = as.numeric(tmatvec)
	#a = as.numeric(matvec)
	tmpmat = tmatvec
	tmpmat_in_REXPOKIT_coo_fmt = mat2coo(tmpmat)
	ia = tmpmat_in_REXPOKIT_coo_fmt[,"ia"]
	ja = tmpmat_in_REXPOKIT_coo_fmt[,"ja"]
	a = tmpmat_in_REXPOKIT_coo_fmt[,"a"]

 # Number of non-zeros
 nz = nrow(Qmat) * ncol(Qmat)

# Run the wrapper function

tmpoutmat = expokit_dgexpv_wrapper(n, m, timeval, v, w, tol, anorm, wsp,
lwsp, iwsp, liwsp, itrace, iflag, ia, ja, a, nz, res)

print(tmpoutmat)
}
\author{
  Nicholas J. Matzke \email{matzke@berkeley.edu}
}
\seealso{
  \code{\link{expokit_dgexpv_Qmat}}

  \code{\link{expokit_wrapalldgexpv_tvals}}
}

