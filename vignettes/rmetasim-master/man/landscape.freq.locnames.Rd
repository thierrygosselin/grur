\name{landscape.freq.locnames}
\alias{landscape.freq.locnames}
\title{Return a vector of locus IDs concatenated to Allele ids}

\description{
  Extracts locus and allele ids from landscape
}

\usage{
  landscape.freq.locnames(Rland)
}
\arguments{
  \item{Rland}{the Rmetasim landscape object}
}

\value{
  character vector
}

\details{

  This function returns a vector of reasonable column names for a
  frequency by individual matrix from an Rland object. It is mainly
  useful in creating genind objects from landscapes.

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
  landscape.ind.freq, landscape.make.genind
  }

\examples{
  exampleland <- landscape.new.example() 
  exampleland <- landscape.simulate(exampleland, 4)
  landscape.freq.locnames(exampleland)
  rm(exampleland)
}
\keyword{misc}
