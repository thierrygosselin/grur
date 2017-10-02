\name{landscape.locusvec}
\alias{landscape.locusvec}
\title{return a vector with the locus ids for each column in the individuals component of a landscape}
\description{return a vector with the locus ids for each column in the individuals component of a landscape}

\usage{
  landscape.locusvec(Rland)
}
\arguments{
  \item{Rland}{the Rmetasim landscape object}
}


\value{
vector
}

\seealso{
  landscape.populations
  }

\examples{
  exampleland <- landscape.new.example()
  exampleland <- landscape.simulate(exampleland, 4)
  landscape.locusvec(exampleland)
  rm(exampleland)
}
\keyword{misc}
