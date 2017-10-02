\name{landscape.populations}
\alias{landscape.populations}
\title{return a vector of population IDs from a landscape}

\description{return a vector of population IDs from a landscape
}
\usage{
  landscape.populations(Rland)
}
\arguments{
  \item{Rland}{the Rmetasim landscape object}
}


\value{
  a vector
}
\details{  Returns a vector of length \code{dim(rland$individuals)[1]} where rland is a
  landscape object.  The vector classifies individuals into populations
  (or habitats)  }

\seealso{
  landscape.locus, landscape.ploidy
  }

\examples{
  exampleland <- landscape.new.example()
  exampleland <- landscape.simulate(exampleland, 4)
  plot(table(landscape.populations(exampleland)),
       main="Distribution of population size in landscape")
  rm(exampleland)
}
\keyword{misc}
