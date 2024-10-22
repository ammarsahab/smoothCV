% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/smoothCV.R
\name{fcCV}
\alias{fcCV}
\title{Forward chaining cross validation}
\usage{
fcCV(
  fullset,
  initialn,
  folds,
  type,
  start,
  end,
  dist,
  localoptim = F,
  alphrange,
  betarange
)
}
\arguments{
\item{fullset}{A set of univariate time series data. Can be a vector (type double) or a data.table.}

\item{initialn}{Number of observations in the initial training set.}

\item{folds}{Number of folds to divide the rest of the data into.}

\item{type}{Type of smoothing. Can be \code{"SMA"}, \code{"DMA"}, \code{"SES"}, \code{"DES"}.}

\item{start}{(if type="SMA" or type="DMA") Starting value for m (see \code{sma.dt} and \code{dma.dt}).}

\item{end}{(if type="SMA" or type="DMA") Maximum value for m. Must be less than or equal to the length of the training set for SMA,
or less than or equal to the length of the training set divided by two for DMA.}

\item{dist}{(if type="SMA" or type="DMA") Distance between successive values for m.
The vector of m values will be constructed as \code{seq(start,end,dist)}.}

\item{localoptim}{Whether to optimize in each training set. If true, a range for paramter values will not need to be provided}

\item{alphrange}{(if type="SES" or type="DES") A vector of parameters for the level component. Can be created using \code{seq}, or by manually specifying a vector.}

\item{betarange}{(if type="DES") A vector of parameters for the trend component. Can be created using \code{seq}, or by manually specifying a vector.
Not used if type="SES"}
}
\value{
A data.table containing all parameter combinations during each iteration (fold) with their respective MSE and MAPE values.
Can be grouped by parameter values to obtain the mean error for each parameter value.
}
\description{
The algorithm starts by picking a number of observations as the initial training set.
The remaining observations will be split into k folds. The first fold will be used as the initial test set.
Grid search is then committed using MA.Grid or ES.Grid.
The first fold is the incorporated to the training set. The second fold will be used as the next test set.
The algorithm repeats itself k times.
}
\examples{
fcCV(fullset=crudenow$Close, initialn=36, folds=12, type="DES", localoptim=F,
alphrange=seq(0.1,1,0.1), betarange=seq(0.1,1,0.1))
fcCV(fullset=crudenow$Close, initialn=36, folds=12, type="SMA", localoptim=F,
start=2, end=30, dist=3)
}
