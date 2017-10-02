\name{landscape.modify.epoch}
\alias{landscape.modify.epoch}
\title{Modifies one of the landscape's epochs}

\description{
This function updates the deomgraphic parameters in a landscape for a
particular epoch
  }
\usage{
landscape.modify.epoch(rland,epoch=1,S=NULL,R=NULL,M=NULL,
   epochprob=NULL,startgen=NULL,extinct=NULL,
   carry=NULL,localprob=NULL)
}
\arguments{
  \item{rland}{landscape object, required}
  \item{epoch}{the epoch to modify, default 1}
  \item{S}{(default=NULL) Survivablity matrix for epoch, NULL leaves unchanged}
  \item{R}{(default=NULL) female Reproduction matrix for epoch, NULL
    leaves unchanged}
  \item{M}{(default=NULL) Male reporduction matrix for epoch, NULL
    leaves unchanged}
  \item{epochprob}{(default=NULL) probability of choosing this epoch, NULL
    leaves unchanged}
  \item{startgen}{(default=NULL) generation in which this epoch starts,
    NULL leaves unchanged}
  \item{extinct}{(default=NULL) vector of extinction probabilities per
    generation for each subpopulation, NULL leaves unchanged}
  \item{carry}{(default=NULL) vector of carrying capacities for each
    subpopulation, must be rland$intparam$habitats in length, NULL
    leaves unchanged}
  \item{localprob}{(default=NULL) vector of probabilites for choosing
    local demographies, must be length(rland$demography$localdem) in
    length, NULL
    leaves unchanged}

}
\keyword{misc}
