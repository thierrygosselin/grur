\name{landscape.amova.locus}
\alias{landscape.amova.locus}

\title{ uses functions in ade4 to calcuate phi-st for a particular locus }
\description{
Runs an amova on a locus.  Does not include information about sequence similarity or ssr size in analysis.  
}
\usage{
landscape.amova.locus(rland, l = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{l}{locus number}
  \item{rland}{ landscape object }
}
\details{
Should be the same as Weir and Cockerham's \eqn{\theta}
}
\value{
list of amova results for a locus
}

\seealso{landscape.amova, landscape.amova.pairwise}
\keyword{misc }% __ONLY ONE__ keyword per line
