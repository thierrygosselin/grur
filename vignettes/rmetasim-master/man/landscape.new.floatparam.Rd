\name{landscape.new.floatparam}
\alias{landscape.new.floatparam}
\title{Create a set of floating point parameters}

\description{
  Create a set of floating point parameters for a Rmetasim landscape.}
\usage{
  ## must be called AFTER landscape.new.empty()
   landscape.new.floatparam(rland,s=0)}
\arguments{
  \item{rland}{skeletion of landscape object, required}
  \item{s}{selfing (default=0), the selfing rate of the species}
}
\examples{
  ## Defaults
  exampleland <- landscape.new.empty()
  exampleland <- landscape.new.floatparam(exampleland)
  exampleland$floatparam

  ## .5 selfing rate
  exampleland <- landscape.new.empty()
  exampleland <- landscape.new.floatparam(exampleland,s=0.5)
  exampleland$floatparam

  rm(exampleland)
}
\keyword{misc}
