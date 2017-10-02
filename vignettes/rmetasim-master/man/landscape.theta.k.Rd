\name{landscape.theta.k}
\alias{landscape.theta.k}
\title{Calculate theta using the number of alleles}

\description{
  Calculate theta using number of alleles from a landscape.}
\usage{
 landscape.theta.k(rland)
}
\arguments{
  \item{rland}{the Rmetasim landscape object}
}
\details{Uses routines in the package 'ape'}
\value{
A matrix with num loci columns and num populations rows.  Each element
reflects the estimated theta for that population x locus combination  }
\seealso{
landscape.theta.h, landscape.theta.s
}
\examples{
  exampleland <- landscape.new.example()
  exampleland <- landscape.simulate(exampleland, 4)

  theta.k.mat <- landscape.theta.k(exampleland)
  theta.k.mat
}
\keyword{misc}
