\name{expokit_mydmexpv_wrapper}
\alias{expokit_mydmexpv_wrapper}
\title{EXPOKIT dmexpv wrapper function, return just output probs}
\usage{
  expokit_mydmexpv_wrapper(n, m, t, v, w, tol, anorm, wsp,
    lwsp, iwsp, liwsp, itrace, iflag, ia, ja, a, nz)
}
\arguments{
  \item{n}{number of rows in Q matrix}

  \item{m}{n-1}

  \item{t}{the value to exponentiate the rate matrix by
  (often e.g. a time value)}

  \item{v}{variable to store some results in; should have n
  elements (and perhaps start with 1)}

  \item{w}{same length as v}

  \item{tol}{tolerance for approximations; usually set to
  0.01}

  \item{anorm}{the norm of the Q matrix}

  \item{lwsp}{length of workspace (wsp); for dmexpv,
  lwsp=n*(m+2)+5*(m+2)^2+ideg+1}

  \item{wsp}{workspace to store some results in; should be
  a double with lwsp elements}

  \item{liwsp}{length of integer workspace; for dmexpv,
  liwsp=m+2}

  \item{iwsp}{integer workspace}

  \item{itrace}{option, set to 0}

  \item{iflag}{option, set to 0}

  \item{ia}{i indices of Qmat nonzero values}

  \item{ja}{j indices of Qmat nonzero values}

  \item{a}{nonzero values of Qmat (ia, ja, a are columns of
  a COO-formatted Q matrix)}

  \item{nz}{number of non-zeros in Qmat

  EXPOKIT needs the input matrix to be transposed compared
  to normal. COO format is required for EXPOKIT.}
}
\value{
  \code{w_output_probs} the output probabilities (=
  \code{myDMEXPV} variable \code{w}, or the fifth output in
  the output from .Call("mydmexpv_", ...), given the
  (first) input t-value.
}
\description{
  This function wraps the .C call to EXPOKIT for the dmexpv
  function.  Only the output probabilities not the Pmat
  probability matrix, are returned.
}
\examples{
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

# DMEXPV and DGEXPV are designed for large, sparse Q matrices (sparse = lots of zeros).
# DMEXPV is specifically designed for Markov chains and so may be slower, but more accurate.

# DMEXPV, single t-value
expokit_wrapalldmexpv_tvals(Qmat=Qmat, tvals=tvals[1], transpose_needed=TRUE)
expokit_wrapalldmexpv_tvals(Qmat=Qmat, tvals=2)

# This function runs a for-loop itself (sadly, we could not get mapply() to work
# on a function that calls dmexpv/dgexpv), returning a list of probability matrices.

# DMEXPV functions
list_of_P_matrices_dmexpv = expokit_wrapalldmexpv_tvals(Qmat=Qmat,
tvals=tvals, transpose_needed=TRUE)
list_of_P_matrices_dmexpv
}
\author{
  Nicholas J. Matzke \email{matzke@berkeley.edu}
}
\seealso{
  \code{\link{expokit_dmexpv_Qmat}}

  \code{\link{expokit_wrapalldmexpv_tvals}}
}

