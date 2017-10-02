\name{landscape.theta.s}
\alias{landscape.theta.s}
\title{Calculate theta using segregating sites}

\description{Calculate theta from a landscape based upon
  the number of seqregating sites.}

\usage{
  landscape.theta.s(rland)
}
\arguments{
  \item{rland}{the Rmetasim landscape object}
}
\details{Uses routines in the package 'ape'}
  
\value{
A matrix with num loci columns and num populations rows.  Each element
reflects the estimated theta for that population x locus combination  }
\seealso{
theta.k.landscape, theta.h.landscape
}
\examples{
  exampleland <- landscape.new.example()
  exampleland <- landscape.simulate(exampleland, 4)
  theta.s.mat <- landscape.theta.s(exampleland)
  theta.s.mat
}
\keyword{misc}
