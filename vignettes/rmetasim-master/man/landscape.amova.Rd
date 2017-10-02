\name{landscape.amova}
\alias{landscape.amova}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ calcuates phi-st for every locus in the landscape }
\description{
 calcuates \eqn{\phi_{ST}} for every locus in the landscape 
}
\usage{
landscape.amova(rland, np = NULL, ns = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{rland}{ landscape object}
  \item{np}{ max number of pops to include }
  \item{ns}{ max number of samples to collect}
}
\value{
vector of length equal to the number of loci
}

\seealso{ \code{\link{landscape.amova.locus}},\code{\link{landscape.amova.pairwise}}}
\keyword{ misc}% __ONLY ONE__ keyword per line
