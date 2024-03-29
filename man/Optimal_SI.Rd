% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Lib_SpectralIndices.R
\name{optimal_SI}
\alias{optimal_SI}
\title{This function computes the correlation between a set of vegetation properties
and all combinations of spectral bands corresponding to a given type of spectral index}
\usage{
optimal_SI(
  ReflMat,
  Spectral_Bands,
  BPvars,
  ExpressIndex,
  Permutations = FALSE,
  nbCPU = 1
)
}
\arguments{
\item{ReflMat}{numeric. reflectance matrix with ncol = nbands and nrows = nsamples}

\item{Spectral_Bands}{numeric. vector containing central wavelength for each spectral band of ReflMat}

\item{BPvars}{numeric. Biophysical properties as matrix with ncol = nBPvars and nrows = nsamples}

\item{ExpressIndex}{character. expression corresponding to the spectral index to be explored ()}

\item{Permutations}{boolean. either compute all permutations, or all combinations}

\item{nbCPU}{numeric. number of CPU for multithreading}
}
\value{
list. includes band Combinations and corresponding correlation between biophysical variables and SI
}
\description{
This function computes the correlation between a set of vegetation properties
and all combinations of spectral bands corresponding to a given type of spectral index
}
