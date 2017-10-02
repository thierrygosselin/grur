\name{landscape.mismatchdist}
\alias{landscape.mismatchdist}
\title{Calculate a mismatch distribution for a locus in a landscape}

\description{Calculate mismatch distribution from a landscape based upon
  the number of seqregating sites.}

\usage{
  landscape.mismatchdist(Rland,lnum=1)
}
\arguments{
  \item{lnum}{locus number to calculate mismatch upon}
  \item{Rland}{the Rmetasim landscape object}
}
\details{Calculates a mismatch distribution for DNA-sequence-based loci.}
  
\value{
A matrix with num loci columns and num populations rows.  Each element
reflects the estimated theta for that population x locus combination  }
\examples{
  exampleland <- landscape.new.example()
  exampleland <- landscape.simulate(exampleland, 4)
  misdist <- landscape.mismatchdist(exampleland,3) #will produce       
  misdist                                          #ridiculous output
}
\keyword{misc}
