% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Lib_SpectralIndices_fun.R
\name{nBR_RAW}
\alias{nBR_RAW}
\title{this function computes nBR_RAW spectral index based on a reflectance data.table
where columns are named after S2 bands: B2, B3, B4, B5, B6, B7, B8, B8A, B11, B12}
\usage{
nBR_RAW(Refl, S2Bands = NULL)
}
\arguments{
\item{Refl}{data.table. each column is a spectral band (B2, B3, B4, B5, B6, B7, B8, B8A, B11, B12)}

\item{S2Bands}{data.table. wavelength}
}
\value{
numeric. nBR_RAW spectral index
}
\description{
this function computes nBR_RAW spectral index based on a reflectance data.table
where columns are named after S2 bands: B2, B3, B4, B5, B6, B7, B8, B8A, B11, B12
}
