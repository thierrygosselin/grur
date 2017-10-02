\name{landscape.clean}
\alias{landscape.clean}
\title{Function to resolve inconsistencies within a landscape}

\description{Converts a landscape to internal format and back.  This can
resolve inconsistencies in a 'hand-built' landscape}
  
  \usage{
   landscape.clean(rland)
}
\arguments{
  \item{rland}{the Rmetasim landscape object}
}

\examples{
  exampleland <- landscape.new.example()
  exampleland <- landscape.simulate(exampleland, 4)
  exampleland.clean <- landscape.clean(exampleland)
  rm(exampleland)
}
\keyword{misc}
