\name{landscape.allelecount}
\alias{landscape.allelecount}
\title{Calculate allele numbers (frequency in the statistical sense) at each locus in each population}

\description{
  Calculate allele counts}
\usage{
  landscape.allelecount(Rland, tbl.out=FALSE)
}
\arguments{
  \item{Rland}{the Rmetasim landscape object}
  \item{tbl.out}{Create a table rather than data.frame output}
}
\value{
  Depends on the value of tbl.out.  See above.
}
\seealso{
landscape.allelefreq, landscape.obs.het, landscape.exp.het, landscape.Fwright, landscape.Fst
}
\examples{
  exampleland <- landscape.new.example()
  exampleland <- landscape.simulate(exampleland, 4)
  landscape.allelefreq(exampleland,tbl.out=TRUE)
  landscape.allelefreq(exampleland,tbl.out=FALSE)
  rm(exampleland)
}
\keyword{misc}
