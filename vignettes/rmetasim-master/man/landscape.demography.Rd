\name{landscape.demography}
\alias{landscape.demography}
\title{Calculate demographic parameters}

\description{
  Calculate demographic parameters from a landscape:  CURRENTLY BROKEN!}
\usage{
   landscape.demography(Rland)
}
\arguments{
  \item{Rland}{the Rmetasim landscape object}
}
\value{
  A list of length populations+1.  The first 1..populations elements
  are lists comprised of lambda, the equilibrium stage-structure, the
  actual stage structure, a \eqn{\chi^2} value for the test of difference
  between predicted and actual, and an estimate of significance for that
  test.  The last element of the main list is the same as the previous
  ones except it refers to the entire landscape
  }

\keyword{misc}
