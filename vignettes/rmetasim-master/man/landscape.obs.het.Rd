\name{landscape.obs.het}
\alias{landscape.obs.het}
\title{Calculate observed heterozygosity}

\description{
  Calculate observed heterozygosity from a landscape}
\usage{
   landscape.obs.het(Rland)
}
\arguments{
  \item{Rland}{the Rmetasim landscape object}
}
\value{
A matrix with num loci columns and num populations rows.  Each element
reflects the observed heterozygosity for that population x locus combination  }

\seealso{
landscape.exp.het, landscape.Fst
}

\examples{
  exampleland <- landscape.new.example()
  exampleland <- landscape.simulate(exampleland, 4)
  obshet <- landscape.obs.het(exampleland)
  rm(exampleland)
}
\keyword{misc}
