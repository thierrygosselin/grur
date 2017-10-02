\name{SimulationComponents}
\alias{landscape.advance}
\alias{landscape.carry}
\alias{landscape.extinct}
\alias{landscape.reproduce}
\alias{landscape.survive}
\title{Code components to simulate a landscape}

\description{
  These functions can be used to construct custom simulations of
  landscapes.  Each conducts only a single generations worth of change}
\usage{
   landscape.advance(Rland, seed = -1)
   landscape.carry(Rland, seed = -1)
   landscape.extinct(Rland, seed = -1)
   landscape.reproduce(Rland, seed = -1)
   landscape.survive(Rland, seed = -1)
}
\arguments{
  \item{Rland}{the Rmetasim landscape object}
  \item{seed}{seed for random number generator}
}
\details{
  \code{landscape.advance()} merely advances the generation counter and selects
  the new generations demographic conditions if such conditions can
  vary.  The other functions implement carrying capacity, local
  extinction, reproduction, and survival/growth, respectively.  The
  function \code{landscape.simulate()} bundles the functionality of these
  components into a single function (and executes it slightly faster all within linked C++ code).
}

\seealso{
  \code{landscape.simulate}
  }

\keyword{misc}
