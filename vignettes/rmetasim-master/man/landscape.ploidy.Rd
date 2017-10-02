\name{landscape.ploidy}
\alias{landscape.ploidy}
\title{return a vector with the ploidy of each locus}
\description{return a vector with the ploidy of each locus in the order they appear in the landscape}

\usage{
  landscape.ploidy(Rland)
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
  landscape.ploidy(exampleland)
  rm(exampleland)
}
\keyword{misc}
