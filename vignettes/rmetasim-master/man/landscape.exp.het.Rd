\name{landscape.exp.het}
\alias{landscape.exp.het}
\title{Calculate expected heterozygosity}

\description{
  Calculate expected heterozygosity from a landscape}
\usage{
   landscape.exp.het(Rland)
}
\arguments{
  \item{Rland}{the Rmetasim landscape object}
}
\value{
A matrix with num loci columns and num populations rows.  Each element
reflects the expected heterozygosity for that population x locus combination  }
\details{
  Calculates the expected heterozygosity in each population:
  \deqn{1-\Sigma_{i_k} p_i^2}{1 - sum(p^2)}
  where \eqn{p} is a vector of allele frequencies for a locus in a population.
}
\seealso{
landscape.obs.het, Fst.landscape
}
\examples{
  exampleland <- landscape.new.example()
  exampleland <- landscape.simulate(exampleland, 4)
  exphet <- landscape.exp.het(exampleland)
  rm(exampleland)
}
\keyword{misc}
