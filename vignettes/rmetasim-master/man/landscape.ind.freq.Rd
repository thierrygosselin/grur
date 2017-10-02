\name{landscape.ind.freq}
\alias{landscape.ind.freq}
\title{Return a matrix containing genotypes in 'frequency per
  individual' format}

\description{
  Converts a landscape into a format that can be used in creating genind objects
}

\usage{
  landscape.ind.freq(Rland, include.states=TRUE)
}
\arguments{
  \item{Rland}{the Rmetasim landscape object}
    \item{include.states}{deprecated}

}

\value{
  matrix
}

\details{

  This function returns a frequency by individual matrix for the
  landscape given in Rland

  Genotypes are stored in a landscape object as columns
  representing allelic position. If the locus is haploid, there is a
  single column and if diploid, two.  The allele indices are entered in
  each column for individual genotypes.

  The 'frequency per individual' format has a column for every
  locusXallele combination; genotypes are represented as the frequency
  of particular alleles in an individuals genotype (for example freq = 1
  in a column diploid means that the individual is homozygous for that
  column, 0.5 in two columns represents a heterozygote for the alleles
  represented by those two columns).  This is essentially the genind
  format from adegenet.

}

\seealso{
  landscape.freq.locnames, landscape.make.genind
  }

\examples{
  exampleland <- landscape.new.example() 
  exampleland <- landscape.simulate(exampleland, 25)
  head(landscape.ind.freq(exampleland))
}
\keyword{misc}
