\name{landscape.coalinput}
\alias{landscape.coalinput}
\title{Add loci and individuals based upon output from SimCoal 2.0}
\author{Mark Bravington and Allan Strand}
\description{
  Take rmetasim object and replaces the locus and individual data based on the results of a SimCoal run stored in Arlequin format files}

\usage{
  ## must be called AFTER integer, switch, and float params have been created
   landscape.coalinput(rland, npp=200, arlseq = NULL, arlms = NULL, seqsitemut=1e-7,
                               msmut = 5e-4, mut.rates = NULL)
}
\arguments{
  \item{rland}{partially created landscape object, required}
  \item{npp}{number per population.  Scalar or vector of length equal to number of populations.  If scalar, value replicated}
  \item{arlseq}{name of the Arlequin format file containing a single locus of haploid sequence data for any number of populations}
\item{arlms}{name of the Arlequin format file containing a single locus of diploid microsatellite data for any number of populations}
  \item{seqsitemut}{mutation rate for sequence data}
  \item{msmut}{mutation rate for diploid genotypic data}
  \item{mut.rates}{alternative means to specify mutation rates.  Legal
    values are either NULL or a vector of rates equal to the number of
    loci to simulate.  If NULL, SSR loci are assigned msmut as a
    mutation rate and sequence-based loci, seqsitemut.  If a vector,
    ovrrides msmut and seqsitemut}
}

\details{This function provides part of an interface between R and
SimCoal, an environment for simulating sequences and microsatellite
genotypes from coalescent trees.  SimCoal can be used to simulate a
standing crop of alleles and their relationships under a wide range of
demographies.  It returns haplotypes and genotypes of individuals in
Arlequin format files.

         If either 'arlseq' or 'arlms' are set to NULL, their
corresponding data will not be included in the landscape (for example
if arlseq=NULL, only diploid genotypes will be imported

         The genotypes in the Arlequin files are used to create
         rland\$loci objects based upon their frequencies and states.
         These rland\$loci sub-objects are then used to populate the
         rland\$individuals sub-object.
                                                          
         The number of populations in the Arlequin files should be the
         same among genetic locus types (sequence versus
         microsatellite) and the rland\$intparam\$habitats parameter.
         The per-population frequency data will be used in creating
         individuals
  }

\value{an rmetasim object with new loci and individuals}
\examples{
  exampleland <- landscape.new.empty()
  exampleland <- landscape.new.intparam(exampleland, s=2, h=2)
  exampleland <- landscape.new.floatparam(exampleland)
  exampleland <- landscape.new.switchparam(exampleland)

#  exampleland <- landscape.coalinput(exampleland)
#  exampleland$loci

}
\keyword{misc}
