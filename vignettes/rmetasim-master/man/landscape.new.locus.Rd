\name{landscape.new.locus}
\alias{landscape.new.locus}
\title{Add a locus}

\description{
  Add a locus to a Rmetasim landscape object}
\usage{
## must be called AFTER integer, switch, and float params have been created
landscape.new.locus(rland,type=0,ploidy=1,
             mutationrate=0,transmission=1,numalleles=2,
             allelesize=50,frequencies=NULL, states = NULL)
}
\arguments{
  \item{rland}{partially created landscape object, required}
  \item{type}{(default=0) type of locus, 0=Infinite Allele mutation model (Integer), 1=Stepwise mutation model (Integer) state, 2=DNA base (variable length string state)}
  \item{ploidy}{(default=1) locus ploidy, 1 or 2}
  \item{mutationrate}{(default=0) probability of mutation per generation, less than or equal to 1}
  \item{transmission}{(default=1) 1=uniparental inheritance, 0=biparental inheritance}
  \item{numalleles}{(default=2) number of different alleles at the time of creation}
  \item{allelesize}{(default=50) length of DNA strings if type=2}
  \item{frequencies}{(default=NULL) vector of frequencies for each
  allele, must be numalleles long and add up to 1, if NULL frequencies
  are equally distributed}
\item{states}{States for each of the alleles specified.  If NULL
  sequential numeric states used.}
}

\examples{
  exampleland <- landscape.new.empty()
  exampleland <- landscape.new.intparam(exampleland, s=2, h=2)
  exampleland <- landscape.new.floatparam(exampleland)
  exampleland <- landscape.new.switchparam(exampleland)

  exampleland <- landscape.new.locus(exampleland,type=2,ploidy=2,
                          mutationrate=.001,numalleles=5,allelesize=100)

  exampleland$loci

  rm(exampleland)
}
\keyword{misc}
