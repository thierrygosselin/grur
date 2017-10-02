\name{landscape.locus}
\alias{landscape.locus}
\title{return a matrix containing genotypes for a particular locus}

\description{return a matrix containing genotypes for a particular locus
}
\usage{
  landscape.locus(Rland, lnum=1)
}
\arguments{
  \item{lnum}{the locus to return}
  \item{Rland}{the Rmetasim landscape object}
}


\value{
  matrix
}

\details{ Returns a matrix with rows = \code{dim(rland$individuals)[1]}.  The
first three columns correspond to the class (and two placeholder
variables) of an individual.  Here rland is a
  landscape object.  The remaining   columns (1 if haploid, 2 if
  diploid) contain the allele indices for the various loci}

\seealso{
  landscape.populations
  }

\examples{
  exampleland <- landscape.new.example()
  exampleland <- landscape.simulate(exampleland, 4)
  print("Allele frequencies at locus 1")
  table(landscape.locus(exampleland,1)[,c(-1:-(landscape.democol()))]) 
  rm(exampleland)
}
\keyword{misc}
