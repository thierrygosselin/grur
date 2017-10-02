\name{landscape.theta.h}
\alias{landscape.theta.h}
\title{Calculate theta using heterozygosity}

\description{Calculate theta from a landscape based upon
  heterozygosity.}

\usage{
  landscape.theta.h(rland)
}
\arguments{
  \item{rland}{the Rmetasim landscape object}
}
\details{Uses routines in the package 'ape'}
  
\value{
A matrix with num loci columns and num populations rows.  Each element
reflects the estimated theta for that population x locus combination  }
\seealso{
landscape.theta.k, landscape.theta.s
}
\examples{

  exampleland <- landscape.new.example()
  exampleland <- landscape.simulate(exampleland, 4)

  theta.h.mat <- landscape.theta.h(exampleland)
  theta.h.mat
}
\keyword{misc}
