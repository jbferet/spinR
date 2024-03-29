% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Lib_SpectralIndices.R
\name{functSI}
\alias{functSI}
\title{this function computes the spectral index defined by idx for Refl}
\usage{
functSI(idx, Refl, S2Bands)
}
\arguments{
\item{idx}{character. name of the index}

\item{Refl}{dataframe reflectance data corresponding to S2 sensor}

\item{S2Bands}{dataframe central wavelength for S2 bands}
}
\value{
numeric. spectral index corresponding to reflectance data
}
\description{
this function computes the spectral index defined by idx for Refl
}
