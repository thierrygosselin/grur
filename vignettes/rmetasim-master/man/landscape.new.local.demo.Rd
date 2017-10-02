\name{landscape.new.local.demo}
\alias{landscape.new.local.demo}
\title{Create a Local Demography}

\description{
  Create a local demography for an Rmetasim Landscape object}
\usage{
  ## must be called AFTER integer, switch, and float params have been created
  ## S, R, and M matricies must be square matricies of size
  ##  rland$intparam$stages by rland$intparam$stages

   landscape.new.local.demo(rland,S,R,M,k=0)}
\arguments{
  \item{rland}{partially created landscape object, required}
  \item{S}{Survivablity matrix for demograpy, required}
  \item{R}{female Reproduction matrix for demography, required}
  \item{M}{Male reporduction matrix for demography, required}
  \item{k}{flag for type of matrix, 0=demgraphy at zero population
    density, 1=demography at carrying capacity}
}

\details{The local demography objects encapsulate demography within a
  particular region.  Multiple such objects can be defined to account
  for different demographies across space.  The flag, k, can indicate
  whether the matrices represent demography at zero population growth
  and at carrying capacity, if density-dependence is modeled}

\examples{
  
  exampleS <- matrix(c(0.1, 0, 0.5, 0.3), nrow = 2)
  exampleR <- matrix(c(0, 1.1, 0, 0), nrow = 2)
  exampleM <- matrix(c(0, 0, 0, 1), nrow = 2)
  
  exampleland <- landscape.new.empty()
  exampleland <- landscape.new.intparam(exampleland, s=2)
  exampleland <- landscape.new.floatparam(exampleland)
  exampleland <- landscape.new.switchparam(exampleland)
  exampleland <- landscape.new.local.demo(exampleland,exampleS,exampleR,exampleM)

  exampleland$demography$localdem

  rm(exampleS)
  rm(exampleR)
  rm(exampleM)
  rm(exampleland)
}
\keyword{misc}
