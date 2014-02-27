\name{fermat.test}
\alias{fermat.test}
\title{Test an integer for primality with Fermat's little theorem.}
\usage{
  fermat.test(n)
}
\arguments{
  \item{n}{the integer to test for primality}
}
\value{
  Whether the integer passes the Fermat test for a
  randomized \eqn{0 < a < n} no workyATcallGraphPrimitives
}
\description{
  Fermat's little theorem states that if \eqn{n} is a prime
  number and \eqn{a} is any positive integer less than
  \eqn{n}, then \eqn{a} raised to the \eqn{n}th power is
  congruent to \eqn{a\ modulo\ n}{a modulo n}.
}
\note{
  \code{fermat.test} doesn't work for integers above
  approximately fifteen because modulus loses precision.
}
\author{
  Peter Danenberg \email{pcd@roxygen.org}
}
\references{
  \url{http://en.wikipedia.org/wiki/Fermat's_little_theorem}
}

