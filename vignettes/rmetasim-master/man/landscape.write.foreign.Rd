\name{landscape.write.foreign}
\alias{landscape.write.foreign}
\title{Save a landscape to a file in a foreign format}

\description{
  Save a Rmetasim landscape object to a file in a suite of output
  formats.  This function used to be supported with backend C++
  functions. These functions have been removed from the codebase.
  Currently this function only supports writing genepop-format files to disk.
}
\usage{
  landscape.write.foreign(Rland,fn="foreign.genepop",fmt="genepop",...)
}
\arguments{
  \item{Rland}{the Rmetasim landscape object}
  \item{fn}{the path and name of the file to save the landscape
    to}
  \item{fmt}{the output format for the landscape:  Can take the
    following values:"genepop"}
  \item{...}{arguments passed to individual functions that actually do
  the work}
}
\details{This is a wrapper for R-side functions that write landscape
  data to disk.

  It is probably not necessary, because of the new
  landscape.make.genind() function.  However, there is a single format
  implemented right now: GenePop.  In addition to the parameters named
  above, a character string called 'title' can be passed to this
  function and it will give a title at the top of the resulting genepop
  file.

  Former versions of this function had a "numi" parameter to specify how
  many individuals to sample per population.  Now use landscape.sample()
  before running this one.
}

\value{None, run only for the side-effect of writing to disk}
\examples{
  ## Needs write access to the current directory, files created!!
  exampleland <- landscape.new.example()
  landscape.write.foreign(exampleland, fn="exampleland.genepop", fmt="genepop")
  rm(exampleland)
}
\keyword{misc}
