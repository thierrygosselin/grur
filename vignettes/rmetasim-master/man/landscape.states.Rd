\name{landscape.states}
\alias{landscape.states}

\title{return a matrix containing actual genotypes for a particular locus}

\description{return a matrix containing the states of the alleles in genotypes for a particular locus
}

\usage{
  landscape.states(Rland,lnum=1)
}

\arguments{
  \item{lnum}{the locus to return}
  \item{Rland}{the Rmetasim landscape object}
}


\value{
  matrix
}

\details{ Returns a matrix with rows = \code{dim(rland$individuals)[1]}.  The
 columns 1:landscape.democol() correspond to demographic variables for an individual.  The columns are: state, placeholder, birthyear, id, mother's id, and father's id.  Here rland is a
  landscape object.  The remaining   columns (1 if haploid, 2 if
  diploid) contain the states of the alleles for the selected loci
}

\seealso{
  landscape.locus
  }

\examples{
  exampleland <- landscape.new.example()
  exampleland <- landscape.simulate(exampleland, 4)
  print("Allele frequencies at locus 1")
  table(landscape.states(exampleland,1)[,c(-1:-landscape.democol())]) 
  rm(exampleland)
}

\keyword{misc}
