\name{landscape.mig.matrix}
\alias{landscape.mig.matrix}
\title{Creates a Migration Matrix for All Life Stages}

\description{Creates a binary matrix representing the migration between
  a set of 'h' populations containing 's' life stages each. This matrix
  can be based on a given migration model or on a custom matrix 
}

\usage{
  landscape.mig.matrix(h=3,s=2,mig.model="island",first.rep.s=s,
  h.dim=NULL, distance.fun=NULL, distance.factor=1, R.custom = NULL, ...)
}

\arguments{

  \item{h}{habitats (default=3), the number of different subpopulations
  within the landscape}

  \item{s}{stages (default=2), the number of stages in the life cycle of
  the organism}

  \item{mig.model}{migration model (default="island"), the migration
  model to be used to make the matrix. Choices are '"island",
  "stepping.stone.linear", "stepping.stone.circular", "twoD",
  "twoDwDiagonal", "distance","custom"'.  See details.}

\item{first.rep.s}{first reproductive life stage (default=s), the life
stage at which the organism starts to reproduce}

  \item{h.dim}{rectangular arrangement of populations
  (default=NULL). vector of length 2 showing the distribution of
  populations in rows and columns when the model of evolution is equal
  to "twoD" or "twoDwDiagonal".}

  \item{distance.fun}{function to calculate migration (default=NULL), an
  user created function that uses the distance between each population
  to calculate the migration rate between those two populations if the
  migration model is equal to "distance".}

  \item{distance.factor}{distance factor (default=1), the distance
  between each adjacent population if the migration model is equal to
  "distance"}

  \item{R.custom}{custom migration matrix (default=NULL), migration
matrix with 'h' by 'h' dimensions to be used to create the larger 'h*s'
by 'h*s' matrix if the migration model is equal to "custom"}

\item{...}{additional arguments passed to 'distance.fun'}
}

\details{ 
This function can work on three different ways:
\enumerate{
  \item{With a given migration model}{This will take in consideration
    one of the predefined migration models to create the migration
    matrix.}
    \itemize{
      \item{"island"}{Migration occurs among all the populations in the
	model.}
      \item{"stepping.stone.linear"}{The populations are distributed
	linearly and migration only occurs between the adjacent populations.}
      \item{"stepping.stone.circular"}{Similar to "stepping.stone.linear",
	but the populations are distributed in a circle so there is migration
	between the first and the last population.}
      \item{"twoD"}{The populations are distributed in two
	dimensions. It is necessary to provide the "h.dim" term in order to
	determine the distribution of the populations in rows and columns
	respectively. Migration only occurs between populations that are
	adjacent to each other}
      \item{"twoDwDiagonal"}{Similar to "twoD", but within a square formed
	by four populations (two rows and two columns) there is migration in
	the diagonal}
    }
  
    \item{With a custom migration matrix}{This requires the user to provide
      the "R.custom" argument. In this case the function will expand the
      migration pattern given on "R.custom" to encompass all life
      stages. For the function to work this way the "mig.model" term must be
      equal to "custom".}  
    \item{With a distance functions}{This requires a function that shows
	how migration changes with changing distance. The "distance.fun" is
	very versatile and the use of "..." allows the functions to accept
	extra terms. The "distance.factor" term allows the user to change the
	distance between the populations to facilitate the use of distance
	functions that work on greater or smaller scales. It is necessary to
	provide the "h.dim". It is necessary to provide the "h.dim" term in order to
	determine the distribution of the populations in rows and columns
	respectively. It is possible to have a linear distribution of
	populations if one of the terms of "h.dim" is equal to 1. For the
	function to work this way the "mig.model" term must be equal to "distance".}
}}
\value{
  \item{R}{Matrix containing the final result from the function call.
  This should be a "h*s" by "h*s" matrix indicating what life stages
  from what populations migrate to the first life stage of what
  populations. When the "mig.model" is equal to distance this matrix
  will indicate the rate of migration between the populations instead of
  if it just occurs or not.}
  \item{h}{the number of different subpopulations}
  \item{s}{the number of stages in the life cycle of the organism}
  \item{mig.model}{the migration model used to make the matrix}
  \item{first.rep.s}{the life stage at which the organism starts to
    reproduce}
  \item{R.int}{A "h" by "h" matrix indicating the migration pattern. If
    "mig.model" is equal to custom, "R.int" will be equal to
    "R.custom".}
}

\author{Artur Veloso and Allan Strand}

%\seealso{\code(plot.mig.matrix) in the package "rmetasim-viz"}

\examples{

#Circular stepping stone migration model
landscape.mig.matrix(s=3,h=4,mig.model="stepping.stone.linear",first.rep.s=2)

#Two dimensions with diagonal migration model
landscape.mig.matrix(h=18,h.dim=c(3,6),s=2,mig.model="twoDwDiagonal")

#Using a custom migration matrix
R.custom <- matrix(c(0, 0, 1, 0,
                     1, 0, 1, 0,
                     1, 0, 0, 0,
                     1, 0, 1, 0), ncol=4,nrow=4,byrow=TRUE)
landscape.mig.matrix(s=3,h=4,first.rep.s=2,mig.model="custom",R.custom=R.custom)

#Using a distance function. Notice that the distance function requires
#the argument "lambda" that can be given in the "make.mig.matrix"
#function call.

my.dist <- function(distance,lambda) {exp(-distance*lambda)}
landscape.mig.matrix(h=18,h.dim=c(3,6),s=2,mig.model="distance",distance.fun=my.dist,lambda=1)
}