#' this function computes ARI1 spectral index based on a reflectance data.table
#' where columns are named after S2 bands: B2, B3, B4, B5, B6, B7, B8, B8A, B11, B12
#'
#' @param Refl data.table. each column is a spectral band (B2, B3, B4, B5, B6, B7, B8, B8A, B11, B12)
#' @param S2Bands data.table. wavelength
#'
#' @return numeric. ARI1 spectral index
#' @export
#'
ARI1 <- function(Refl, S2Bands=NULL){
  SI <- (1/Refl$B3)-(1/Refl$B5)
  return(SI)
}

#' this function computes ARI2 spectral index based on a reflectance data.table
#' where columns are named after S2 bands: B2, B3, B4, B5, B6, B7, B8, B8A, B11, B12
#'
#' @param Refl data.table. each column is a spectral band (B2, B3, B4, B5, B6, B7, B8, B8A, B11, B12)
#' @param S2Bands data.table. wavelength
#'
#' @return numeric. ARI2 spectral index
#' @export
#'
ARI2 <- function(Refl, S2Bands=NULL){
  SI <- (Refl$B8/Refl$B2)-(Refl$B8/Refl$B3)
  return(SI)
}

#' this function computes ARVI spectral index based on a reflectance data.table
#' where columns are named after S2 bands: B2, B3, B4, B5, B6, B7, B8, B8A, B11, B12
#'
#' @param Refl data.table. each column is a spectral band (B2, B3, B4, B5, B6, B7, B8, B8A, B11, B12)
#' @param S2Bands data.table. wavelength
#'
#' @return numeric. ARVI spectral index
#' @export
#'
ARVI <- function(Refl, S2Bands=NULL){
  SI <- (Refl$B8-(2*Refl$B4-Refl$B2))/(Refl$B8+(2*Refl$B4-Refl$B2))
  return(SI)
}

#' this function computes BAI spectral index based on a reflectance data.table
#' where columns are named after S2 bands: B2, B3, B4, B5, B6, B7, B8, B8A, B11, B12
#'
#' @param Refl data.table. each column is a spectral band (B2, B3, B4, B5, B6, B7, B8, B8A, B11, B12)
#' @param S2Bands data.table. wavelength
#'
#' @return numeric. BAI spectral index
#' @export
BAI <- function(Refl, S2Bands=NULL){
  SI <- (1/((0.1-Refl$B4)**2+(0.06-Refl$B8)**2))
  return(SI)
}

#' this function computes BAIS2 spectral index based on a reflectance data.table
#' where columns are named after S2 bands: B2, B3, B4, B5, B6, B7, B8, B8A, B11, B12
#'
#' @param Refl data.table. each column is a spectral band (B2, B3, B4, B5, B6, B7, B8, B8A, B11, B12)
#' @param S2Bands data.table. wavelength
#'
#' @return numeric. BAIS2 spectral index
#' @export
#'
BAIS2 <- function(Refl, S2Bands=NULL){
  SI <-  (1-((Refl$B6 * Refl$B7 * Refl$B8A)/Refl$B4)**0.5) *
    ((Refl$B12 - Refl$B8A)/((Refl$B12 + Refl$B8A)**0.5)+1)
  return(SI)
}

#' this function computes CCCI spectral index based on a reflectance data.table
#' where columns are named after S2 bands: B2, B3, B4, B5, B6, B7, B8, B8A, B11, B12
#'
#' @param Refl data.table. each column is a spectral band (B2, B3, B4, B5, B6, B7, B8, B8A, B11, B12)
#' @param S2Bands data.table. wavelength
#'
#' @return numeric. CCCI spectral index
#' @export
#'
CCCI <- function(Refl, S2Bands=NULL){
  SI <- ((Refl$B8 - Refl$B5) / (Refl$B8 + Refl$B5)) /
    ((Refl$B8 - Refl$B4) / (Refl$B8 + Refl$B4))
  return(SI)
}

#' this function computes CHL_RE spectral index based on a reflectance data.table
#' where columns are named after S2 bands: B2, B3, B4, B5, B6, B7, B8, B8A, B11, B12
#'
#' @param Refl data.table. each column is a spectral band (B2, B3, B4, B5, B6, B7, B8, B8A, B11, B12)
#' @param S2Bands data.table. wavelength
#'
#' @return numeric. CHL_RE spectral index
#' @export
#'
CHL_RE <- function(Refl, S2Bands=NULL){
  SI <- Refl$B5/Refl$B8
  return(SI)
}

#' this function computes CRI1 spectral index based on a reflectance data.table
#' where columns are named after S2 bands: B2, B3, B4, B5, B6, B7, B8, B8A, B11, B12
#'
#' @param Refl data.table. each column is a spectral band (B2, B3, B4, B5, B6, B7, B8, B8A, B11, B12)
#' @param S2Bands data.table. wavelength
#'
#' @return numeric. CRI1 spectral index
#' @export
#'
CRI1 <- function(Refl, S2Bands=NULL){
  SI <- (1/Refl$B2)-(1/Refl$B3)
  return(SI)
}

#' this function computes CRI2 spectral index based on a reflectance data.table
#' where columns are named after S2 bands: B2, B3, B4, B5, B6, B7, B8, B8A, B11, B12
#'
#' @param Refl data.table. each column is a spectral band (B2, B3, B4, B5, B6, B7, B8, B8A, B11, B12)
#' @param S2Bands data.table. wavelength
#'
#' @return numeric. CRI2 spectral index
#' @export
#'
CRI2 <- function(Refl, S2Bands=NULL){
  SI <- (1/Refl$B2)-(1/Refl$B5)
  return(SI)
}

#' this function computes EVI spectral index based on a reflectance data.table
#' where columns are named after S2 bands: B2, B3, B4, B5, B6, B7, B8, B8A, B11, B12
#'
#' @param Refl data.table. each column is a spectral band (B2, B3, B4, B5, B6, B7, B8, B8A, B11, B12)
#' @param S2Bands data.table. wavelength
#'
#' @return numeric. EVI spectral index
#' @export
#'
EVI <- function(Refl, S2Bands=NULL){
  SI <- 2.5*(Refl$B8-Refl$B4)/((Refl$B8+6*Refl$B4-7.5*Refl$B2+1))
  return(SI)
}

#' this function computes EVI2 spectral index based on a reflectance data.table
#' where columns are named after S2 bands: B2, B3, B4, B5, B6, B7, B8, B8A, B11, B12
#'
#' @param Refl data.table. each column is a spectral band (B2, B3, B4, B5, B6, B7, B8, B8A, B11, B12)
#' @param S2Bands data.table. wavelength
#'
#' @return numeric. EVI2 spectral index
#' @export
#'
EVI2 <- function(Refl, S2Bands=NULL){
  SI <- 2.5*(Refl$B8-Refl$B4)/(Refl$B8+2.4*Refl$B4+1)
  return(SI)
}

#' this function computes GRVI1 spectral index based on a reflectance data.table
#' where columns are named after S2 bands: B2, B3, B4, B5, B6, B7, B8, B8A, B11, B12
#'
#' @param Refl data.table. each column is a spectral band (B2, B3, B4, B5, B6, B7, B8, B8A, B11, B12)
#' @param S2Bands data.table. wavelength
#'
#' @return numeric. GRVI1 spectral index
#' @export
#'
GRVI1 <- function(Refl, S2Bands=NULL){
  SI <- (Refl$B4-Refl$B3)/(Refl$B4+Refl$B3)
  return(SI)
}

#' this function computes GRVI1 spectral index based on a reflectance data.table
#' where columns are named after S2 bands: B2, B3, B4, B5, B6, B7, B8, B8A, B11, B12
#'
#' @param Refl data.table. each column is a spectral band (B2, B3, B4, B5, B6, B7, B8, B8A, B11, B12)
#' @param S2Bands data.table. wavelength
#'
#' @return numeric. GRVI1 spectral index
#' @export
#'
GRVI1 <- function(Refl, S2Bands=NULL){
  SI <- (Refl$B4-Refl$B3)/(Refl$B4+Refl$B3)
  return(SI)
}

#' this function computes GNDVI spectral index based on a reflectance data.table
#' where columns are named after S2 bands: B2, B3, B4, B5, B6, B7, B8, B8A, B11, B12
#'
#' @param Refl data.table. each column is a spectral band (B2, B3, B4, B5, B6, B7, B8, B8A, B11, B12)
#' @param S2Bands data.table. wavelength
#'
#' @return numeric. GNDVI spectral index
#' @export
#'
GNDVI <- function(Refl, S2Bands=NULL){
  SI <- (Refl$B8-Refl$B3)/(Refl$B8+Refl$B3)
  return(SI)
}

#' this function computes IRECI spectral index based on a reflectance data.table
#' where columns are named after S2 bands: B2, B3, B4, B5, B6, B7, B8, B8A, B11, B12
#'
#' @param Refl data.table. each column is a spectral band (B2, B3, B4, B5, B6, B7, B8, B8A, B11, B12)
#' @param S2Bands data.table. wavelength
#'
#' @return numeric. IRECI spectral index
#' @export
#'
IRECI <- function(Refl, S2Bands=NULL){
  SI <- (Refl$B7-Refl$B4)/(Refl$B5/(Refl$B6))
  return(SI)
}

#' this function computes LAI_SAVI spectral index based on a reflectance data.table
#' where columns are named after S2 bands: B2, B3, B4, B5, B6, B7, B8, B8A, B11, B12
#'
#' @param Refl data.table. each column is a spectral band (B2, B3, B4, B5, B6, B7, B8, B8A, B11, B12)
#' @param S2Bands data.table. wavelength
#'
#' @return numeric. LAI_SAVI spectral index
#' @export
#'
LAI_SAVI <- function(Refl, S2Bands=NULL){
  SI <- -log(0.371 + 1.5 * (Refl$B8 - Refl$B4) / (Refl$B8+ Refl$B4+ 0.5)) / 2.4
  return(SI)
}

#' this function computes MCARI spectral index based on a reflectance data.table
#' where columns are named after S2 bands: B2, B3, B4, B5, B6, B7, B8, B8A, B11, B12
#'
#' @param Refl data.table. each column is a spectral band (B2, B3, B4, B5, B6, B7, B8, B8A, B11, B12)
#' @param S2Bands data.table. wavelength
#'
#' @return numeric. MCARI spectral index
#' @export
#'
MCARI <- function(Refl, S2Bands=NULL){
  SI  <- (1-0.2*(Refl$B5-Refl$B3)/(Refl$B5-Refl$B4))
  return(SI)
}

#' this function computes mNDVI705 spectral index based on a reflectance data.table
#' where columns are named after S2 bands: B2, B3, B4, B5, B6, B7, B8, B8A, B11, B12
#'
#' @param Refl data.table. each column is a spectral band (B2, B3, B4, B5, B6, B7, B8, B8A, B11, B12)
#' @param S2Bands data.table. wavelength
#'
#' @return numeric. mNDVI705 spectral index
#' @export
#'
mNDVI705 <- function(Refl, S2Bands=NULL){
  SI <- (Refl$B6-Refl$B5)/(Refl$B6+Refl$B5-2*Refl$B2)
  return(SI)
}

#' this function computes MSAVI2 spectral index based on a reflectance data.table
#' where columns are named after S2 bands: B2, B3, B4, B5, B6, B7, B8, B8A, B11, B12
#'
#' @param Refl data.table. each column is a spectral band (B2, B3, B4, B5, B6, B7, B8, B8A, B11, B12)
#' @param S2Bands data.table. wavelength
#'
#' @return numeric. MSAVI2 spectral index
#' @export
#'
MSAVI2 <- function(Refl, S2Bands=NULL){
  SI <- ((Refl$B8+1)-0.5*sqrt(((2*Refl$B8)-1)**2+8*Refl$B4))
  return(SI)
}

#' this function computes MSI spectral index based on a reflectance data.table
#' where columns are named after S2 bands: B2, B3, B4, B5, B6, B7, B8, B8A, B11, B12
#'
#' @param Refl data.table. each column is a spectral band (B2, B3, B4, B5, B6, B7, B8, B8A, B11, B12)
#' @param S2Bands data.table. wavelength
#'
#' @return numeric. MSI spectral index
#' @export
#'
MSI <- function(Refl, S2Bands=NULL){
  SI <- Refl$B11/Refl$B8
  return(SI)
}

#' this function computes mSR705 spectral index based on a reflectance data.table
#' where columns are named after S2 bands: B2, B3, B4, B5, B6, B7, B8, B8A, B11, B12
#'
#' @param Refl data.table. each column is a spectral band (B2, B3, B4, B5, B6, B7, B8, B8A, B11, B12)
#' @param S2Bands data.table. wavelength
#'
#' @return numeric. mSR705 spectral index
#' @export
#'
mSR705 <- function(Refl, S2Bands=NULL){
  SI <- (Refl$B6-Refl$B2)/(Refl$B5-Refl$B2)
  return(SI)
}

#' this function computes MTCI spectral index based on a reflectance data.table
#' where columns are named after S2 bands: B2, B3, B4, B5, B6, B7, B8, B8A, B11, B12
#'
#' @param Refl data.table. each column is a spectral band (B2, B3, B4, B5, B6, B7, B8, B8A, B11, B12)
#' @param S2Bands data.table. wavelength
#'
#' @return numeric. MTCI spectral index
#' @export
#'
MTCI <- function(Refl, S2Bands=NULL){
  SI <- (Refl$B6-Refl$B5)/(Refl$B5+Refl$B4)
  return(SI)
}

#' this function computes nBR_RAW spectral index based on a reflectance data.table
#' where columns are named after S2 bands: B2, B3, B4, B5, B6, B7, B8, B8A, B11, B12
#'
#' @param Refl data.table. each column is a spectral band (B2, B3, B4, B5, B6, B7, B8, B8A, B11, B12)
#' @param S2Bands data.table. wavelength
#'
#' @return numeric. nBR_RAW spectral index
#' @export
#'
nBR_RAW <- function(Refl, S2Bands=NULL){
  SI <- (Refl$B8-Refl$B12)/(Refl$B8+Refl$B12)
  return(SI)
}

#' this function computes NDI_45 spectral index based on a reflectance data.table
#' where columns are named after S2 bands: B2, B3, B4, B5, B6, B7, B8, B8A, B11, B12
#'
#' @param Refl data.table. each column is a spectral band (B2, B3, B4, B5, B6, B7, B8, B8A, B11, B12)
#' @param S2Bands data.table. wavelength
#'
#' @return numeric. NDI_45 spectral index
#' @export
#'
NDI_45 <- function(Refl, S2Bands=NULL){
  SI <- (Refl$B5-Refl$B4)/(Refl$B5+Refl$B4)
  return(SI)
}

#' this function computes NDII spectral index based on a reflectance data.table
#' where columns are named after S2 bands: B2, B3, B4, B5, B6, B7, B8, B8A, B11, B12
#'
#' @param Refl data.table. each column is a spectral band (B2, B3, B4, B5, B6, B7, B8, B8A, B11, B12)
#' @param S2Bands data.table. wavelength
#'
#' @return numeric. NDII spectral index
#' @export
#'
NDII <- function(Refl, S2Bands=NULL){
  SI <- (Refl$B8-Refl$B11)/(Refl$B8+Refl$B11)
  return(SI)
}

#' this function computes NDSI spectral index based on a reflectance data.table
#' where columns are named after S2 bands: B2, B3, B4, B5, B6, B7, B8, B8A, B11, B12
#'
#' @param Refl data.table. each column is a spectral band (B2, B3, B4, B5, B6, B7, B8, B8A, B11, B12)
#' @param S2Bands data.table. wavelength
#'
#' @return numeric. NDSI spectral index
#' @export
#'
NDSI <- function(Refl, S2Bands=NULL){
  SI <- (Refl$B3 - Refl$B11) / (Refl$B3 + Refl$B11)
  return(SI)
}

#' this function computes NDVI spectral index based on a reflectance data.table
#' where columns are named after S2 bands: B2, B3, B4, B5, B6, B7, B8, B8A, B11, B12
#'
#' @param Refl data.table. each column is a spectral band (B2, B3, B4, B5, B6, B7, B8, B8A, B11, B12)
#' @param S2Bands data.table. wavelength
#'
#' @return numeric. NDVI spectral index
#' @export
#'
NDVI <- function(Refl, S2Bands=NULL){
  SI <- (Refl$B8-Refl$B4)/(Refl$B8+Refl$B4)
  return(SI)
}

#' this function computes NDVI_G spectral index based on a reflectance data.table
#' where columns are named after S2 bands: B2, B3, B4, B5, B6, B7, B8, B8A, B11, B12
#'
#' @param Refl data.table. each column is a spectral band (B2, B3, B4, B5, B6, B7, B8, B8A, B11, B12)
#' @param S2Bands data.table. wavelength
#'
#' @return numeric. NDVI_G spectral index
#' @export
#'
NDVI_G <- function(Refl, S2Bands=NULL){
  SI <- Refl$B3*(Refl$B8-Refl$B4)/(Refl$B8+Refl$B4)
  return(SI)
}

#' this function computes NDVI705 spectral index based on a reflectance data.table
#' where columns are named after S2 bands: B2, B3, B4, B5, B6, B7, B8, B8A, B11, B12
#'
#' @param Refl data.table. each column is a spectral band (B2, B3, B4, B5, B6, B7, B8, B8A, B11, B12)
#' @param S2Bands data.table. wavelength
#'
#' @return numeric. NDVI705 spectral index
#' @export
#'
NDVI705 <- function(Refl, S2Bands=NULL){
  SI <- (Refl$B6-Refl$B5)/(Refl$B6+Refl$B5)
  return(SI)
}

#' this function computes NDWI1 spectral index based on a reflectance data.table
#' where columns are named after S2 bands: B2, B3, B4, B5, B6, B7, B8, B8A, B11, B12
#'
#' @param Refl data.table. each column is a spectral band (B2, B3, B4, B5, B6, B7, B8, B8A, B11, B12)
#' @param S2Bands data.table. wavelength
#'
#' @return numeric. NDWI1 spectral index
#' @export
#'
NDWI1 <- function(Refl, S2Bands=NULL){
  SI <- (Refl$B8A-Refl$B11)/(Refl$B8A+Refl$B11)
  return(SI)
}

#' this function computes NDWI2 spectral index based on a reflectance data.table
#' where columns are named after S2 bands: B2, B3, B4, B5, B6, B7, B8, B8A, B11, B12
#'
#' @param Refl data.table. each column is a spectral band (B2, B3, B4, B5, B6, B7, B8, B8A, B11, B12)
#' @param S2Bands data.table. wavelength
#'
#' @return numeric. NDWI2 spectral index
#' @export
#'
NDWI2 <- function(Refl, S2Bands=NULL){
  SI <- (Refl$B8A-Refl$B12)/(Refl$B8A+Refl$B12)
  return(SI)
}

#' this function computes PSRI spectral index based on a reflectance data.table
#' where columns are named after S2 bands: B2, B3, B4, B5, B6, B7, B8, B8A, B11, B12
#'
#' @param Refl data.table. each column is a spectral band (B2, B3, B4, B5, B6, B7, B8, B8A, B11, B12)
#' @param S2Bands data.table. wavelength
#'
#' @return numeric. PSRI spectral index
#' @export
#'
PSRI <- function(Refl, S2Bands=NULL){
  SI <- (Refl$B4-Refl$B2)/(Refl$B5)
  return(SI)
}

#' this function computes PSRI_NIR spectral index based on a reflectance data.table
#' where columns are named after S2 bands: B2, B3, B4, B5, B6, B7, B8, B8A, B11, B12
#'
#' @param Refl data.table. each column is a spectral band (B2, B3, B4, B5, B6, B7, B8, B8A, B11, B12)
#' @param S2Bands data.table. wavelength
#'
#' @return numeric. PSRI_NIR spectral index
#' @export
#'
PSRI_NIR <- function(Refl, S2Bands=NULL){
  SI <- (Refl$B4-Refl$B2)/(Refl$B8)
  return(SI)
}

#' this function computes RE_NDVI spectral index based on a reflectance data.table
#' where columns are named after S2 bands: B2, B3, B4, B5, B6, B7, B8, B8A, B11, B12
#'
#' @param Refl data.table. each column is a spectral band (B2, B3, B4, B5, B6, B7, B8, B8A, B11, B12)
#' @param S2Bands data.table. wavelength
#'
#' @return numeric. RE_NDVI spectral index
#' @export
#'
RE_NDVI <- function(Refl, S2Bands=NULL){
  SI <- (Refl$B8-Refl$B6)/(Refl$B8+Refl$B6)
  return(SI)
}

#' this function computes RE_NDWI spectral index based on a reflectance data.table
#' where columns are named after S2 bands: B2, B3, B4, B5, B6, B7, B8, B8A, B11, B12
#'
#' @param Refl data.table. each column is a spectral band (B2, B3, B4, B5, B6, B7, B8, B8A, B11, B12)
#' @param S2Bands data.table. wavelength
#'
#' @return numeric. RE_NDWI spectral index
#' @export
#'
RE_NDWI <- function(Refl, S2Bands=NULL){
  SI <- (Refl$B4-Refl$B6)/(Refl$B4+Refl$B6)
  return(SI)
}

#' this function computes S2REP spectral index based on a reflectance data.table
#' where columns are named after S2 bands: B2, B3, B4, B5, B6, B7, B8, B8A, B11, B12
#'
#' @param Refl data.table. each column is a spectral band (B2, B3, B4, B5, B6, B7, B8, B8A, B11, B12)
#' @param S2Bands data.table. wavelength
#'
#' @return numeric. S2REP spectral index
#' @export
#'
S2REP <- function(Refl, S2Bands=NULL){
  SI <- 705+35*(((0.5*(Refl$B7+Refl$B4))-Refl$B6)/(Refl$B7-Refl$B6))
  return(SI)
}

#' this function computes SAVI spectral index based on a reflectance data.table
#' where columns are named after S2 bands: B2, B3, B4, B5, B6, B7, B8, B8A, B11, B12
#'
#' @param Refl data.table. each column is a spectral band (B2, B3, B4, B5, B6, B7, B8, B8A, B11, B12)
#' @param S2Bands data.table. wavelength
#'
#' @return numeric. SAVI spectral index
#' @export
#'
SAVI <- function(Refl, S2Bands=NULL){
  SI <- 1.5*(Refl$B8-Refl$B4)/(Refl$B8+Refl$B4+0.5)
  return(SI)
}

#' this function computes SIPI spectral index based on a reflectance data.table
#' where columns are named after S2 bands: B2, B3, B4, B5, B6, B7, B8, B8A, B11, B12
#'
#' @param Refl data.table. each column is a spectral band (B2, B3, B4, B5, B6, B7, B8, B8A, B11, B12)
#' @param S2Bands data.table. wavelength
#'
#' @return numeric. SIPI spectral index
#' @export
#'
SIPI <- function(Refl, S2Bands=NULL){
  SI <- (Refl$B8-Refl$B2)/(Refl$B8-Refl$B4)
  return(SI)
}

#' this function computes SR spectral index based on a reflectance data.table
#' where columns are named after S2 bands: B2, B3, B4, B5, B6, B7, B8, B8A, B11, B12
#'
#' @param Refl data.table. each column is a spectral band (B2, B3, B4, B5, B6, B7, B8, B8A, B11, B12)
#' @param S2Bands data.table. wavelength
#'
#' @return numeric. SR spectral index
#' @export
#'
SR <- function(Refl, S2Bands=NULL){
  SI <- Refl$B8/Refl$B4
  return(SI)
}

#' this function computes CR_SWIR spectral index based on a reflectance data.table
#' where columns are named after S2 bands: B2, B3, B4, B5, B6, B7, B8, B8A, B11, B12
#'
#' @param Refl data.table. each column is a spectral band (B2, B3, B4, B5, B6, B7, B8, B8A, B11, B12)
#' @param S2Bands data.table. wavelength
#'
#' @return numeric. CR_SWIR spectral index
#' @export
#'
CR_SWIR <- function(Refl, S2Bands){
  SI <- Refl$B11/(Refl$B8A+(S2Bands$B11-S2Bands$B8A)*(Refl$B12-Refl$B8A)/(S2Bands$B12-S2Bands$B8A))
  return(SI)
}

#' this function computes CR_RE spectral index based on a reflectance data.table
#' where columns are named after S2 bands: B2, B3, B4, B5, B6, B7, B8, B8A, B11, B12
#'
#' @param Refl data.table. each column is a spectral band (B2, B3, B4, B5, B6, B7, B8, B8A, B11, B12)
#' @param S2Bands data.table. wavelength
#'
#' @return numeric. CR_RE spectral index
#' @export
#'
CR_RE <- function(Refl, S2Bands){
  SI <- Refl$B5/(Refl$B4+(S2Bands$B5-S2Bands$B4)*(Refl$B6-Refl$B4)/(S2Bands$B6-S2Bands$B4))
  return(SI)
}
