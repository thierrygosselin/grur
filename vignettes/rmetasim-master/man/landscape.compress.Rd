\name{landscape.compress}
\alias{landscape.compress}
\title{Function to resolve inconsistencies within a landscape, deprecated}

\description{
  Deprecated, 'clean.landscape()' does the same.
}
\usage{
   landscape.compress(Rland)
}
\arguments{
  \item{Rland}{the Rmetasim landscape object}
}

\examples{
  exampleland <- landscape.new.example()
  exampleland <- landscape.simulate(exampleland, 4)
  exampleland.clean <- landscape.compress(exampleland)
  rm(exampleland)
}
\keyword{misc}
