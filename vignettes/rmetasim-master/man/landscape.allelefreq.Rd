\name{landscape.allelefreq}
\alias{landscape.allelefreq}
\title{Calculate allele frequencies at each locus in each population}

\description{
  Calculate allele frequencies}
\usage{
   landscape.allelefreq(Rland,tbl.out=FALSE)
}
\arguments{
  \item{Rland}{the Rmetasim landscape object}
  \item{tbl.out}{return as a (three-dimensional) table if TRUE.  If
  FALSE, return as a dataframe with categorical variables denoting the
  locus, population and allele.}
}
\value{
  Depends on the value of tbl.out.  See above.
}
\seealso{
landscape.obs.het, landscape.exp.het, landscape.Fwright, landscape.Fst
}
\examples{
  exampleland <- landscape.new.example()
  exampleland <- landscape.simulate(exampleland, 4)
  landscape.allelefreq(exampleland,tbl.out=TRUE)
  landscape.allelefreq(exampleland,tbl.out=FALSE)
  rm(exampleland)
}
\keyword{misc}
