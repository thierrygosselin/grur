\name{landscape.new.epoch}
\alias{landscape.new.epoch}
\title{Create an Epoch}

\description{
  Create an epoch for a Rmetasim landscape object}
\usage{
  ## must be called AFTER integer, switch, and float params have
  ##   been created and after the demography has been created
  ## S, R, and M matricies must be square matricies of size X by X
  ##   where X = rland$intparam$stages*rland$intparam$habitats

landscape.new.epoch(rland,S=NULL,R=NULL,M=NULL,
epochprob=1,startgen=0,extinct=NULL,carry=NULL,
localprob=NULL)
}
\arguments{
  \item{rland}{partially created landscape object, required}
  \item{S}{(default=NULL) Survivablity matrix for epoch, NULL gives no movement between subpopulations (0 matrix)}
  \item{R}{(default=NULL) female Reproduction matrix for epoch, NULL gives no dispersal between subpopulations (0 matrix)}
  \item{M}{(default=NULL) Male reporduction matrix for epoch, NULL gives
    no sperm or pollen movement between subpopulations (0 matrix)}
  \item{epochprob}{(default=1) probability of choosing this epoch randomly if randepoch==1}
  \item{startgen}{(default=0) generation in which this epoch starts}
  \item{extinct}{(default=NULL) vector of extinction probabilities per
    generation for each subpopulation, must be rland$intparam$habitats
    in length, passing NULL gives a 0\% probability of extinction to each subpopulation}
    \item{carry}{(default=NULL) vector of carrying capacities for each subpopulation, must be rland$intparam$habitats in length, passing NULL gives a 1000 individual carrying capacity to each subpopulation}
    \item{localprob}{(default=NULL) vector of probabilites for choosing local demographies, must be length(rland$demography$localdem) in length, passing NULL gives each demography an equal probability}
}

  \examples{
  
  exampleS <- matrix(c(0.1, 0, 0.5, 0.3), nrow = 2)
  exampleR <- matrix(c(0, 1.1, 0, 0), nrow = 2)
  exampleM <- matrix(c(0, 0, 0, 1), nrow = 2)
  
  exampleland <- landscape.new.empty()
  exampleland <- landscape.new.intparam(exampleland, s=2, h=2)
  exampleland <- landscape.new.floatparam(exampleland)
  exampleland <- landscape.new.switchparam(exampleland)
  exampleland <- landscape.new.local.demo(exampleland,exampleS,exampleR,exampleM)

  ## nonsense matricies
  exampleS <- matrix(c(rep(0,4),
                rep(1,4),
                rep(0,4),
                rep(1,4)), nrow = 4)
  exampleR <- matrix(c(rep(0.5,4),
                rep(0,4),
                rep(0.5,4),
                rep(0,4)), nrow = 4)
  exampleM <- matrix(c(rep(0,4),
                rep(.25,4),
                rep(0,4),
                rep(0,4)), nrow = 4)

  ## defaults
  exampleland<- landscape.new.epoch(exampleland,exampleS,exampleR,exampleM)

  exampleland$demography$epochs[[1]]

  rm(exampleS)
  rm(exampleR)
  rm(exampleM)
  rm(exampleland)
}
\keyword{misc}
