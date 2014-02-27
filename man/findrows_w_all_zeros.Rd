\name{findrows_w_all_zeros}
\alias{findrows_w_all_zeros}
\title{Check if a Q matrix has rows with all zeros}
\usage{
  findrows_w_all_zeros(matvec)
}
\arguments{
  \item{matvec}{Q transition matrix}
}
\value{
  A list of TRUE/FALSE, as long as the number of rows.
  \code{TRUE}=the is all zeros, \code{FALSE}=the row has
  nonzero values.
}
\description{
  Q matrices with all-zero rows will crash
  .Call(wrapalldmexpv_, ...) and .Call(wrapalldgexpv_,
  ...), and therefore will crash
  expokit_wrapalldmexpv_tvals() and
  expokit_wrapalldgexpv_tvals() when these are set (the
  default) to return the full P matrix.  These functions
  work fine with zero rows if \code{inputprobs_for_fast} is
  supplied, meaning that only the output probabilities of
  each state are returned.
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

# DGEXPV, single t-value
expokit_wrapalldgexpv_tvals(Qmat=Qmat, tvals=tvals[1], transpose_needed=TRUE)
expokit_wrapalldgexpv_tvals(Qmat=Qmat, tvals=2)

# This function runs the for-loop itself (sadly, we could not get mapply() to work
# on a function that calls dmexpv/dgexpv), returning a list of probability matrices.

# DGEXPV functions
list_of_P_matrices_dgexpv = expokit_wrapalldgexpv_tvals(Qmat=Qmat,
tvals=tvals, transpose_needed=TRUE)
list_of_P_matrices_dgexpv
}
\author{
  Nicholas J. Matzke \email{matzke@berkeley.edu}
}
\seealso{
  \code{\link{expokit_wrapalldmexpv_tvals}}

  \code{\link{expokit_wrapalldgexpv_tvals}}
}

