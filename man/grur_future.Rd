% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/internal.R
\name{grur_future}
\alias{grur_future}
\title{grur parallel function}
\usage{
grur_future(
  .x,
  .f,
  flat.future = c("int", "chr", "dfr", "dfc", "walk", "drop"),
  split.vec = FALSE,
  split.with = NULL,
  split.chunks = 4L,
  parallel.core = parallel::detectCores() - 1,
  ...
)
}
\description{
Updating grur to use future
}
\keyword{internal}
