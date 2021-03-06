% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{DR_statistic_time}
\alias{DR_statistic_time}
\title{DR statistic with sampling times}
\usage{
DR_statistic_time(x, sample_times)
}
\arguments{
\item{x}{an ape::phylo object}

\item{sample_times}{a vector of ages}
}
\value{
matrix of rates for each tip, by sample time
}
\description{
Computes the Jetz DR rate, but sample at certain times in the psat
}
