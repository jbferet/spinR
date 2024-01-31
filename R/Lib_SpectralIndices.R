# ============================================================================== =
# prosail
# Lib_SpectralIndices.R
# ============================================================================== =
# PROGRAMMERS:
# Jean-Baptiste FERET <jb.feret@teledetection.fr>
# Florian de BOISSIEU <fdeboiss@gmail.com>
# Copyright 2019/11 Jean-Baptiste FERET
# ============================================================================== =
# This Library includes aims at computing spectral indices from reflectance data
# ============================================================================== =

#' This function computes Area under curve for continuum removed reflectances
#' See Malenovský et al. (2013) for details
#' http://dx.doi.org/10.1016/j.rse.2012.12.015
#'
#' @param Refl RasterBrick, RasterStack or list. Raster bands in the order of SensorBands.
#' @param SensorBands numeric. vector containing central wavelength for each spectral band in the image
#' @param AUCminmax list. wavelengths of lower and upper boundaries ('CRmin' and 'CRmax')
#' @param ReflFactor numeric. multiplying factor used to write reflectance in image (==10000 for S2)
#'
#' @return AUCval raster
#' @export
AUC <- function(Refl, SensorBands, AUCminmax, ReflFactor=1){

  AUCbands <- list()
  AUCbands[['CRmin']] <- SensorBands[get_closest_bands(SensorBands,AUCminmax[['CRmin']])]
  AUCbands[['CRmax']] <- SensorBands[get_closest_bands(SensorBands,AUCminmax[['CRmax']])]
  Bands <- get_closest_bands(SensorBands,AUCbands)
  for (i in Bands[['CRmin']]:Bands[['CRmax']]){
    if (is.na(match(i,Bands))){
      AUCbands[[paste('B',i,sep = '')]] <- SensorBands[i]
    }
  }
  # compute continuum removal for all spectral bands
  CR <- CR_WL(Refl = Refl, SensorBands = SensorBands,
              CRbands = AUCbands, ReflFactor=ReflFactor)

  WL <- sort(unlist(AUCbands),decreasing = F)
  AUCval <- 0.5*(1-CR[[1]])*(WL[2]-WL[1])
  for (i in 2:length(CR)){
    AUCval <- AUCval+0.5*(2-CR[[i-1]]-CR[[i]])*(WL[i+1]-WL[i])
  }
  AUCval <- AUCval+0.5*(1-CR[[length(CR)]])*(WL[i+2]-WL[i+1])
  return(AUCval)
}

#' This function computes continuum removal value for a spectral band of interest,
#' based on lower and upper wavelengths corresponding to boundaries of the continuum
#'
#' @param WLmin numeric. wavelength of the spectral band corresponding to minimum boundary
#' @param WLmax numeric. wavelength of the spectral band corresponding to maximum boundary
#' @param WLtarget numeric. wavelength of the spectral band for which CR is computed
#' @param boundaries list. raster objects corresponding to minimum and maximum wavelengths
#' @param target list. raster object corresponding target wavelength
#' @param p list. in case progress bar with parallel computing
#'
#' @return CR list. raster object corresponding to continuum removed value
#' @export
compute_CR <- function(WLmin, WLmax, WLtarget, boundaries, target, p = NULL){
  CR <- target/(boundaries[[1]]+(WLtarget-WLmin)*(boundaries[[2]]-boundaries[[1]])/(WLmax-WLmin))
  if (!is.null(p)){p()}
  return(CR)
}

#' this function produces a spectral index from an expression defining a spectral index
#'
#' @param listBands list. list of spectral bands defined in the 'ExpressIndex' variable
#' @param Refl RasterBrick, RasterStack or list. Raster bands in the order of SensorBands.
#' @param SensorBands numeric. wavelength in nanometers of the spectral bands of Refl.
#' @param ExpressIndex  character. expression corresponding to the spectral index to compute
#' @param ReflFactor numeric. multiplying factor used to write reflectance in image (==10000 for S2)
#' @param NameIndex character. name for the index to be computed, provided in the raster layer
#' @param p list. in case progress bar with parallel computing
#'
#' @return numeric. band numbers of original sensor corresponding to S2
#' @importFrom raster subset stack
#' @export
compute_SI_fromExp <- function(listBands, Refl, SensorBands, ExpressIndex , ReflFactor=1,
                                                  NameIndex = NULL, p = NULL){

  # define which bands to be used in the spectral index
  Bands <- get_closest_bands(SensorBands,listBands)

  # manage depending on data type
  ClassInputData <- class(Refl)[1]
  if (ClassInputData=='RasterBrick' | ClassInputData=='RasterStack' | ClassInputData=='stars'){
    if (ClassInputData=='stars'){
      Refl <- Refl[Bands]
    } else {
      Refl <- raster::subset(Refl, Bands)
    }
  } else if(ClassInputData=='matrix'){
    if (dim(Refl)[2] == length(SensorBands)){
      Refl <- Refl[,Bands]
    } else if (dim(Refl)[1] == length(SensorBands)){
      Refl <- t(Refl[Bands,])
    } else {
      stop('The dimensions of the reflectance matrix do not match with those of SensorBands')
    }
    # convert matrix into list
    Refl <- lapply(seq_len(ncol(Refl)), function(i) Refl[,i])
  } else if(is.list(Refl)){
    Refl <- raster::stack(Refl[Bands]) # checks that all rasters have same crs/extent
  } else {
    stop('Refl is expected to be a RasterStack, RasterBrick, Stars object, list of rasters or matrix')
  }
  # apply reflectance factor
  if (!ReflFactor==1){
    Refl <- Refl/ReflFactor
  }
  # add band names
  names(Refl) <- gsub(pattern = 'B',replacement = 'Band',x = names(Bands))

  # get number of bands involved in the computation of the spectral index
  nbBands <- unique(as.numeric(gsub(pattern = 'B',
                                    replacement = '',
                                    x =  unlist(regmatches(ExpressIndex,
                                                           gregexpr("B[[:digit:]]+",
                                                                    ExpressIndex))))))
  sortBand <- sort(nbBands,index.return=T,decreasing = T)

  # produce an expression corresponding to the computation of the spectral index
  matches <- unique(unlist(regmatches(ExpressIndex, gregexpr("B[[:digit:]]+", ExpressIndex))))[sortBand$ix]
  replaces <- paste("Refl[['Band",gsub(pattern = 'B',replacement = '',x = matches),"']]",sep = '')
  ExpressIndex_Final <- ExpressIndex
  for (bb in 1:length(matches)){
    ExpressIndex_Final <- gsub(pattern = matches[bb], replacement = replaces[bb], x = ExpressIndex_Final)
  }
  SI <- eval(parse(text = ExpressIndex_Final))
  if (! ClassInputData=='matrix' && !is.null(NameIndex)){
    names(SI) <- NameIndex
  } else if (ClassInputData=='matrix' && !is.null(NameIndex)){
    SI <- data.frame(SI)
    colnames(SI) <- NameIndex
  }
  if (!is.null(p)){p()}
  return(SI)
}

#' This function computes the correlation between a set of vegetation properties
#' and a spectral index
#'
#' @param listBands numeric. bands used in the expression
#' @param Refl numeric. reflectance matrix with ncol = nbands and nrows = nsamples
#' @param SensorBands numeric. vector containing central wavelength for each spectral band of ReflMat
#' @param ExpressIndex character. expression corresponding to the spectral index to be explored ()
#' @param BPvars numeric. Biophysical properties as matrix with ncol = nBPvars and nrows = nsamples
#' @param ReflFactor numeric. multplying factor for reflectance
#' @param NameIndex character. name of spectral index
#'
#' @return list. includes band Combinations and corresponding correlation between biophysical variables and SI
#' @importFrom Rfast correls
#' @export

compute_SI_fromExp_Corr <- function(listBands, Refl, SensorBands, ExpressIndex, BPvars,
                                    ReflFactor = 1, NameIndex = 'SI'){
  SIopt <- compute_SI_fromExp(listBands = listBands,
                              Refl = Refl,
                              SensorBands = SensorBands,
                              ExpressIndex = ExpressIndex,
                              ReflFactor = ReflFactor,
                              NameIndex = NameIndex)
  corr_BP_SI <- Rfast::correls(x = BPvars,y = SIopt$SI,type = 'pearson')[,1]
  return(corr_BP_SI)
}

#' This function computes the correlation between a set of vegetation properties
#' and a continuum removed spectral index
#'
#' @param listBands numeric. bands used in the expression
#' @param Refl numeric. reflectance matrix with ncol = nbands and nrows = nsamples
#' @param SensorBands numeric. vector containing central wavelength for each spectral band of ReflMat
#' @param BPvars numeric. Biophysical properties as matrix with ncol = nBPvars and nrows = nsamples
#'
#' @return list. includes band Combinations and corresponding correlation between biophysical variables and SI
#' @importFrom Rfast correls
#' @export

compute_SI_CR_Corr <- function(listBands, Refl, SensorBands, BPvars){

  if (dim(Refl)[1]==length(SensorBands)){
    Refl <- t(Refl)
  }
  bands <- get_closest_bands(SensorBands,listBands)
  SIopt <- compute_CR(WLmin = as.numeric(listBands[1]), WLmax = as.numeric(listBands[3]),
                      WLtarget = as.numeric(listBands[2]),
                      boundaries = list(Refl[,bands[1]], Refl[,bands[3]]), target = Refl[,bands[2]])
  corr_BP_SI <- Rfast::correls(x = BPvars,y = SIopt,type = 'pearson')[,1]
  return(corr_BP_SI)
}

#' this function aims at computing spectral indices from Sensor reflectance data in raster object
#' it computes the spectral indices based on their computation with Sentinel-2
#' and assumes that the bands of the S2 data follow this order
#' wavelength	= {496.6, 560.0, 664.5, 703.9, 740.2, 782.5, 835.1, 864.8, 1613.7, 2202.4}
#' Full description of the indices can be found here:
#' https://www.sentinel-hub.com/eotaxonomy/indices
#'
#' @param Refl RasterBrick, RasterStack or list. Raster bands in the order of SensorBands.
#' @param SensorBands numeric. wavelength in nanometers of the spectral bands of Refl.
#' @param Sel_Indices  list. list of spectral indices to be computed
#' @param nbCPU  numeric. number of CPU if parallel processing required
#' @param StackOut logical. If TRUE returns a stack, otherwise a list of rasters.
#' @param ReflFactor numeric. multiplying factor used to write reflectance in image (==10000 for S2)
#' @param Offset numeric. offset (when Refl between 0 and 1) to be applied on reflectance. Useful to avoid zero values
#' @param S2Bands numeric. wavelength of the spectral bands corresponding to S2 (default = S2A)
#'
#' @return list. includes
#' - SpectralIndices: List of spectral indices computed from the reflectance initially provided
#' - listIndices: list of spectral indices computable with the function
#' @importFrom methods is
#' @importFrom raster stack brick subset
#' @importFrom future.apply future_lapply
#' @importFrom future plan multisession sequential
#' @export

compute_S2SI_Raster <- function(Refl, SensorBands, Sel_Indices='ALL', nbCPU = 1,
                                StackOut = T, ReflFactor = 1, Offset = 0,
                                S2Bands = data.frame('B2'=492.7, 'B3'=559.8, 'B4'=664.6,
                                                     'B5'=704.1, 'B6'=740.5, 'B7' = 782.8,
                                                     'B8' = 832.8, 'B8A' = 864.7,
                                                     'B11' = 1613.7, 'B12' = 2202.4)){
  SpectralIndices <- list()
  Sen2S2 <- get_closest_bands(SensorBands,S2Bands)
  ClassRaster <- class(Refl)[1]
  if (ClassRaster=='RasterBrick' | ClassRaster=='RasterStack' | ClassRaster=='stars' | ClassRaster=='SpatRaster'){
    if (ClassRaster=='stars'){
      Refl <- Refl[Sen2S2]
    } else if (ClassRaster=='SpatRaster'){
      Refl <- Refl[[Sen2S2]]
    } else {
      Refl <- raster::subset(Refl, Sen2S2)
    }
  } else if(is.list(Refl)){
    Refl <- raster::stack(Refl[Sen2S2]) # checks that all rasters have same crs/extent
  } else {
    stop('Refl is expected to be a RasterStack, RasterBrick, Stars object or a list of rasters')
  }
  # if !ReflFactor == 1 then apply a reflectance factor
  if (!ReflFactor==1){
    Refl <- Refl/ReflFactor
  }
  # if offset
  if (!Offset==0){
    Refl <- Refl+Offset
  }
  names(Refl) <- names(Sen2S2)
  IndexAll <- list()
  # # set zero vaues to >0 in order to avoid problems
  # SelZero <- which(Refl==0)
  # Refl[SelZero] <- 0.005
  # if (dim(Refl)[1]==length(SensorBands)){
  #   Refl <- t(Refl)
  # }

  # inelegant but meeeeh
  listIndices <- list('ARI1','ARI2','ARVI','BAI','BAIS2','CCCI','CHL_RE','CRI1','CRI2','EVI','EVI2',
                      'GRVI1','GNDVI','IRECI','LAI_SAVI','MCARI','mNDVI705','MSAVI2',
                      'MSI','mSR705','MTCI','nBR_RAW','NDI_45','NDII','NDSI','NDVI','NDVI_G',
                      'NDVI705','NDWI1','NDWI2','PSRI','PSRI_NIR','RE_NDVI','RE_NDWI','S2REP',
                      'SAVI','SIPI','SR','CR_SWIR','CR_RE')
  if (Sel_Indices[1]=='ALL'){
    Sel_Indices <- listIndices
  }

  functSI <- function(idx, Refl, S2Bands){
    fsi <- get(idx)
    SpectralIndices <- fsi(Refl, S2Bands)
    return(SpectralIndices)
  }

  if (nbCPU>1){
    plan(multisession, workers = nbCPU) ## Parallelize using four cores
    SpectralIndices <- future_lapply(Sel_Indices,
                                     FUN = functSI,
                                     Refl = Refl,
                                     S2Bands = S2Bands)
    plan(sequential)
  } else {
    SpectralIndices <- lapply(Sel_Indices,
                              FUN = functSI,
                              Refl = Refl,
                              S2Bands = S2Bands)
  }
  names(SpectralIndices) <- Sel_Indices
  if(StackOut)
    SpectralIndices <- raster::stack(SpectralIndices)

  res <- list('SpectralIndices'=SpectralIndices,'listIndices'=listIndices)
  return(res)
}

#' this function aims at computing spectral indices from Sensor reflectance data.
#' it computes the spectral indices based on their computation with Sentinel-2
#' and assumes that the bands of the S2 data follow this order
#' wavelength	= {496.6, 560.0, 664.5, 703.9, 740.2, 782.5, 835.1, 864.8, 1613.7, 2202.4}
#' Full description of the indices can be found here:
#' https://www.sentinel-hub.com/eotaxonomy/indices
#'
#' @param Refl numeric. Reflectance dataset defined in matrix
#' @param SensorBands numeric. wavelength of the spectral bands corresponding to reflectance
#' @param Sel_Indices list. list of spectral indices to be computed
#' @param nbCPU  numeric. number of CPU if parallel processing required
#' @param S2Bands numeric. wavelength of the spectral bands corresponding to S2 (default = S2A)
#'
#' @return list. includes
#' - SpectralIndices: List of spectral indices computed from the reflectance initially provided
#' - listIndices: list of spectral indices computable with the function
#' @importFrom data.table data.table
#' @importFrom snow splitRows
#' @importFrom future plan multisession sequential
#' @importFrom future.apply future_lapply
#' @export

compute_S2SI_from_Sensor <- function(Refl, SensorBands, Sel_Indices='ALL', nbCPU = 1,
                                     S2Bands = data.frame('B2'=492.7, 'B3'=559.8, 'B4'=664.6,
                                                          'B5'=704.1, 'B6'=740.5, 'B7' = 782.8,
                                                          'B8' = 832.8, 'B8A' = 864.7,
                                                          'B11' = 1613.7, 'B12' = 2202.4)){


  SpectralIndices <- list()
  Sen2S2 <- get_closest_bands(SensorBands,S2Bands)
  # set zero vaues to >0 in order to avoid problems
  Refl <- data.table::data.table(Refl)
  names(Refl) <- names(Sen2S2)
  Refl[Refl==0] <- 0.0001
  if (dim(Refl)[1]==length(SensorBands)){
    Refl <- t(Refl)
  }

  # inelegant but meeeeh
  listIndices <- list('ARI1','ARI2','ARVI','BAI','BAIS2','CCCI','CHL_RE','CRI1','CRI2','EVI','EVI2',
                      'GRVI1','GNDVI','IRECI','LAI_SAVI','MCARI','mNDVI705','MSAVI2',
                      'MSI','mSR705','MTCI','nBR_RAW','NDI_45','NDII','NDSI','NDVI','NDVI_G',
                      'NDVI705','NDWI1','NDWI2','PSRI','PSRI_NIR','RE_NDVI','RE_NDWI','S2REP',
                      'SAVI','SIPI','SR','CR_SWIR','CR_RE')
  if (Sel_Indices[1]=='ALL'){
    Sel_Indices = listIndices
  }

  functSI <- function(idx, Refl, S2Bands){
    fsi <- get(idx)
    SpectralIndices <- fsi(Refl, S2Bands)
    return(SpectralIndices)
  }

  if (nbCPU>1){
    plan(multisession, workers = nbCPU) ## Parallelize using four cores
    SpectralIndices <- future_lapply(Sel_Indices,
                                     FUN = functSI,
                                     Refl = Refl,
                                     S2Bands = S2Bands)
    plan(sequential)
  } else {
    SpectralIndices <- lapply(Sel_Indices,
                              FUN = functSI,
                              Refl = Refl,
                              S2Bands = S2Bands)
  }
  names(SpectralIndices) <- Sel_Indices
  res <- list('SpectralIndices'=SpectralIndices,'listIndices'=listIndices)
  return(res)
}

#' This function computes spectral indices from raster data
#'
#' @param input_raster_path character. raster path / list of raster paths
#' @param input_rast_wl numeric. vector containing central wavelength for each spectral band in the image
#' @param output_dir character. main directory where SI will be written
#' @param output_rasters character. list of paths corresponding to each path
#' @param SI_list character. spectral indices to be computed (cf package spinR)
#' @param input_mask_path character. path for mask
#' @param ReflFactor numeric. multiplying factor used to write reflectance in image (==10000 for S2)
#' @param nbCPU numeric. number of CPU for parallel processing
#' @param filetype character. driver
#' @param overwrite boolean. set TRUE to overwrite existing files
#'
#' @return output_rasters list. list of path for SI rasters
#' @export

spectralindices_from_raster <- function(input_raster_path, input_rast_wl,
                                        output_dir, output_rasters = NULL,
                                        SI_list, input_mask_path =  NULL,
                                        ReflFactor = 10000, nbCPU = 1,
                                        filetype = 'COG', overwrite = TRUE){

  if (length(SI_list)==1){
    if (SI_list== 'ALL') SI_list <- list('ARI1', 'ARI2', 'ARVI', 'BAI', 'BAIS2',
                                         'CCCI', 'CHL_RE', 'CRI1', 'CRI2', 'EVI',
                                         'EVI2', 'GRVI1', 'GNDVI', 'IRECI',
                                         'LAI_SAVI', 'MCARI', 'mNDVI705',
                                         'MSAVI2', 'MSI', 'mSR705', 'MTCI',
                                         'nBR_RAW', 'NDI_45', 'NDII', 'NDSI',
                                         'NDVI', 'NDVI_G', 'NDVI705', 'NDWI1',
                                         'NDWI2', 'PSRI', 'PSRI_NIR', 'RE_NDVI',
                                         'RE_NDWI', 'S2REP', 'SAVI', 'SIPI', 'SR',
                                         'CR_SWIR', 'CR_RE')
  }
  input_args <- list('SI' = SI_list,
                     'ReflFactor' = ReflFactor,
                     'nbCPU' = nbCPU,
                     'SensorBands' = input_rast_wl)
  input_rasters <- list('img' = input_raster_path,
                        'mask' = input_mask_path)
  if (is.null(input_mask_path)) input_rasters$mask <- NULL
  # if name for individual output rasters was not provided
  if (is.null(output_rasters)) {
    output_rasters <- file.path(output_dir, SI_list)
    # add extension .tif if tif or COG
    if (filetype%in%c('GTiff', 'COG')) output_rasters <- paste0(output_rasters,
                                                                '.tif')
    output_rasters <- as.list(output_rasters)
  }
  # name rasters based on name of spectral indices
  names(output_rasters) <- SI_list
  # check if already exists and processed
  # if no overwrite, still checks if some of the files are missing
  if (overwrite==FALSE){
    fileExists <- lapply(output_rasters, file.exists)
    for (SI in names(output_rasters)){
      if (fileExists[[SI]]==TRUE) output_rasters[[SI]] <- NULL
    }
  }
  if (length(output_rasters)>0){
    output_lyrs <- 1
    unit_nrows <- FALSE # if need to process spatial units
    funct <- bigRaster::wrapperBig_SI
    bandNames <- as.list(SI_list)
    names(bandNames) <- SI_list
    output_rasters <- bigRaster::apply_bigRaster(funct = funct,
                                                 input_rasters = input_rasters,
                                                 input_args = input_args,
                                                 output_rasters = output_rasters,
                                                 output_lyrs = output_lyrs,
                                                 filetype = filetype,
                                                 bandNames = bandNames)
  }
  return(output_rasters)
}


#' This function extracts boundaries to be used to compute continuum from reflectance data
#'
#' @param Refl RasterBrick, RasterStack or list. Raster bands in the order of SensorBands.
#' @param SensorBands numeric. vector containing central wavelength for each spectral band in the image
#' @param CRbands list. list of spectral bands (central wavelength) including CRmin and CRmax
#' @param ReflFactor numeric. multiplying factor used to write reflectance in image (==10000 for S2)
#'
#' @return CRminmax list. list of rasters corresponding to minimum and maximum wavelengths
#' @export
CR_bound <- function(Refl, SensorBands, CRbands, ReflFactor=1){

  # get closest spectral bands from CR1 and CR2
  Bands <- get_closest_bands(SensorBands,list(CRbands[['CRmin']],CRbands[['CRmax']]))
  WL <- SensorBands[Bands]
  # get equation for line going from CR1 to CR2
  CRminmax <- readRasterBands(Refl = Refl, Bands = Bands, ReflFactor=ReflFactor)
  names(CRminmax) <- paste('WL_',as.character(WL),sep = '')
  return(CRminmax)
}

#' This function extracts boundaries to be used to compute continuum from reflectance data
#'
#' @param Refl RasterBrick, RasterStack or list. Raster bands in the order of SensorBands.
#' @param SensorBands numeric. vector containing central wavelength for each spectral band in the image
#' @param CRbands list. list of spectral bands (central wavelength) including CRmin and CRmax
#' @param ReflFactor numeric. multiplying factor used to write reflectance in image (==10000 for S2)
#'
#' @return outlier_IQR numeric. band numbers of original sensor corresponding to S2
#' @importFrom progress progress_bar
#' @export
CR_WL <- function(Refl, SensorBands, CRbands, ReflFactor=1){

  # Make sure CRmin and CRmax are correctly defined
  if (is.na(match('CRmin',names(CRbands))) | is.na(match('CRmax',names(CRbands)))){
    stop('Please define CRmin and CRmax (CRmin<CRmax) as spectral bands in CRbands')
  }
  if (CRbands[['CRmax']] < CRbands[['CRmin']]){
    stop('Please define CRmin < CRmax in CRbands')
  }
  # extract CRmin and CRmax
  CRminmax <- CR_bound(Refl, SensorBands, CRbands, ReflFactor=ReflFactor)
  # extract other bands and compute CR
  CRmin <- SensorBands[get_closest_bands(SensorBands,CRbands[['CRmin']])]
  CRmax <- SensorBands[get_closest_bands(SensorBands,CRbands[['CRmax']])]
  CRbands[['CRmin']] <- NULL
  CRbands[['CRmax']] <- NULL
  CR <- list()
  # initiate progress bar
  pgbarlength <- length(CRbands)
  pb <- progress_bar$new(
    format = "Computing continuum removal [:bar] :percent in :elapsedfull , estimated time remaining :eta",
    total = pgbarlength, clear = FALSE, width= 100)
  # computation for each band
  for (band in CRbands){
    pb$tick()
    bandrank <- get_closest_bands(SensorBands,band)
    raster2CR <- readRasterBands(Refl = Refl, Bands = bandrank, ReflFactor=ReflFactor)
    CR[[as.character(band)]] <- compute_CR(WLmin = CRmin, WLmax = CRmax,
                                           WLtarget = band, boundaries=CRminmax,
                                           target=raster2CR)
  }
  return(CR)
}

#' this function identifies the bands of a given sensor with closest match to its spectral characteristics
#'
#' @param SensorBands numeric. wavelength in nanometer of the sensor of interest
#' @param listBands numeric or list. Named vector or list of spectral bands corresponding to sensor
#'
#' @return numeric. band numbers of original sensor
#' @export
get_closest_bands <- function(SensorBands,listBands){
  sapply(listBands, function(x){b = which.min(abs(SensorBands-x)); names(b)=''; b})
}

#' This function computes interquartile range (IQR) criterion, which can be used
#' as a criterion for outlier detection
#'
#' @param DistVal numeric. vector of distribution of values
#' @param weightIRQ numeric. weighting factor appplied to IRQ to define lower and upper boudaries for outliers
#'
#' @return outlier_IQR numeric. band numbers of original sensor corresponding to S2
#' @importFrom stats IQR quantile
#' @export
IQR_outliers <- function(DistVal,weightIRQ = 1.5){
  iqr <- IQR(DistVal, na.rm=TRUE)
  range_IQR <- c(quantile(DistVal, 0.25,na.rm=TRUE),quantile(DistVal, 0.75,na.rm=TRUE))
  outlier_IQR <- c(range_IQR[1]-weightIRQ*iqr,range_IQR[2]+weightIRQ*iqr)
  return(outlier_IQR)
}


#' This function computes the correlation between a set of vegetation properties
#' and all combinations of spectral bands corresponding to a given type of spectral index
#'
#' @param ReflMat numeric. reflectance matrix with ncol = nbands and nrows = nsamples
#' @param Spectral_Bands numeric. vector containing central wavelength for each spectral band of ReflMat
#' @param BPvars numeric. Biophysical properties as matrix with ncol = nBPvars and nrows = nsamples
#' @param ExpressIndex character. expression corresponding to the spectral index to be explored ()
#' @param Permutations boolean. either compute all permutations, or all combinations
#' @param nbCPU numeric. number of CPU for multithreading
#'
#' @return list. includes band Combinations and corresponding correlation between biophysical variables and SI
#' @importFrom gtools combinations permutations
#' @import cli
#' @importFrom progressr progressor with_progress handlers
#' @importFrom Rfast correls
#' @importFrom snow splitRows
#' @importFrom future.apply future_lapply
#' @importFrom future plan multisession sequential
#' @export

optimal_SI <- function(ReflMat, Spectral_Bands, BPvars, ExpressIndex,Permutations = FALSE, nbCPU = 1){

  # number of bands in reflectance data
  nbBands <- length(Spectral_Bands)
  # number of bands in spectral index
  nbBandsSI <- unique(as.numeric(gsub(pattern = 'B',
                                      replacement = '',
                                      x =  unlist(regmatches(ExpressIndex,
                                                             gregexpr("B[[:digit:]]+",ExpressIndex))))))
  sortBand <- sort(nbBandsSI,index.return=T,decreasing = F)
  matches <- unique(unlist(regmatches(ExpressIndex, gregexpr("B[[:digit:]]+", ExpressIndex))))[sortBand$ix]
  if (Permutations==FALSE){
    BandCombs <- gtools::combinations(n = nbBands,r = length(matches),repeats.allowed = F)
  } else {
    BandCombs <- gtools::permutations(n = nbBands,r = length(matches),repeats.allowed = F)
  }
  Bands <- do.call(cbind,lapply(seq_len(ncol(BandCombs)), function(i) Spectral_Bands[BandCombs[,i]]))
  colnames(Bands) <- matches
  Bands <- data.frame(Bands)
  if (nbCPU==1){
    SIopt <- sub_SI_Corr(Bands = Bands,
                         Refl = ReflMat,
                         SensorBands = Spectral_Bands,
                         ExpressIndex = ExpressIndex,
                         BPvars = BPvars,
                         ReflFactor = 1,
                         p = NULL)
  } else {
    # multithread
    if ((length(Bands[[1]])/20000)>(nbCPU*10)){
      nbmultithread <- round((length(Bands[[1]])/20000))
    } else {
      nbmultithread <- nbCPU*10
    }
    Bands2 <- snow::splitRows(Bands, nbmultithread)
    plan(multisession, workers = nbCPU)
    handlers(global = TRUE)
    handlers("cli")
    with_progress({
      p <- progressr::progressor(steps = nbmultithread)
      SIopt <- future_lapply(Bands2,
                             FUN = sub_SI_Corr,
                             Refl = ReflMat,
                             SensorBands = Spectral_Bands,
                             ExpressIndex = ExpressIndex,
                             BPvars = BPvars,
                             ReflFactor = 1,
                             p = p,
                             future.packages = c("Rfast","pbapply"))
    })
    plan(sequential)
    SIopt <- do.call(rbind,SIopt)
  }
  res <- list('BandCombinations' = Bands, 'Correlation' = SIopt)
  return(res)
}

#' This function computes the correlation between a set of vegetation properties
#' and all combinations of spectral bands corresponding to a given type of spectral index
#'
#' @param ReflMat numeric. reflectance matrix with ncol = nbands and nrows = nsamples
#' @param Spectral_Bands numeric. vector containing central wavelength for each spectral band of ReflMat
#' @param BPvars numeric. Biophysical properties as matrix with ncol = nBPvars and nrows = nsamples
#' @param nbCPU numeric. number of CPU for multithreading
#'
#' @return list. includes band Combinations and corresponding correlation between biophysical variables and SI
#' @import cli
#' @importFrom gtools combinations
#' @importFrom progressr progressor with_progress handlers
#' @importFrom Rfast correls
#' @importFrom pbapply pblapply
#' @importFrom future plan multisession sequential
#' @importFrom future.apply future_lapply
#' @export

optimal_SI_CR <- function(ReflMat, Spectral_Bands, BPvars, nbCPU = 1){

  # number of bands in reflectance data
  nbBands <- length(Spectral_Bands)
  # number of bands in spectral index
  nbBandsSI <- 3
  matches <- c('B1', 'B2', 'B3')
  BandCombs <- gtools::combinations(n = nbBands,r = length(matches),repeats.allowed = F)
  Bands <- do.call(cbind,lapply(seq_len(ncol(BandCombs)), function(i) Spectral_Bands[BandCombs[,i]]))
  colnames(Bands) <- matches
  Bands <- data.frame(Bands)
  if (nbCPU==1){
    SIopt <- sub_SI_CR_Corr(Bands = Bands,
                            Refl = ReflMat,
                            SensorBands = Spectral_Bands,
                            BPvars = BPvars)
  } else {
    Bands2 <- snow::splitRows(Bands, nbCPU*50)
    # multithread
    plan(multisession, workers = nbCPU)
    handlers(global = TRUE)
    handlers("cli")
    with_progress({
      p <- progressr::progressor(steps = nbCPU*50)
      SIopt <- future_lapply(Bands2,
                             FUN = sub_SI_CR_Corr,
                             Refl = ReflMat,
                             SensorBands = Spectral_Bands,
                             BPvars = BPvars, p,
                             future.packages = c("Rfast","pbapply"))
    })
    plan(sequential)
    SIopt <- do.call(rbind,SIopt)
  }
  res <- list('BandCombinations' = Bands, 'Correlation' = SIopt)
  return(res)
}

#' This function selects bands from a raster or stars object
#'
#' @param Refl RasterBrick, RasterStack or list. Raster bands in the order of SensorBands.
#' @param Bands numeric. rank of bands to be read in Refl
#' @param ReflFactor numeric. multiplying factor used to write reflectance in image (==10000 for S2)
#'
#' @return Robj list. R object (default is raster, stars if Refl is stars object)
#' @importFrom raster subset stack
#' @export
readRasterBands <- function(Refl, Bands, ReflFactor=1){

  # get equation for line going from CR1 to CR2
  ClassRaster <- class(Refl)
  if (ClassRaster=='RasterBrick' | ClassRaster=='RasterStack' | ClassRaster=='stars'){
    if (ClassRaster=='stars'){
      Robj <- Refl[Bands]
    } else {
      Robj <- raster::subset(Refl, Bands)
    }
  } else if(is.list(Refl)){
    Robj <- raster::stack(Refl[Bands]) # checks that all rasters have same crs/extent
  } else {
    stop('Refl is expected to be a RasterStack, RasterBrick, Stars object or a list of rasters')
  }
  # if !ReflFactor == 1 then apply a reflectance factor
  if (!ReflFactor==1){
    Robj <- Robj/ReflFactor
  }
  return(Robj)
}

#' This function computes the correlation between a set of vegetation properties
#' and a set of combinations of spectral bands corresponding to a given type of spectral index
#'
#' @param Bands numeric.
#' @param Refl numeric. reflectance matrix with ncol = nbands and nrows = nsamples
#' @param SensorBands numeric. vector containing central wavelength for each spectral band of ReflMat
#' @param ExpressIndex  character. expression corresponding to the spectral index to compute
#' @param BPvars numeric. Biophysical properties as matrix with ncol = nBPvars and nrows = nsamples
#' @param ReflFactor numeric. multiplying factor used to write reflectance in image (==10000 for S2)
#' @param p function.
#'
#' @return SIopt list. includes band Combinations and corresponding correlation between biophysical variables and SI
#' @importFrom pbapply pblapply
#' @export

sub_SI_Corr <- function(Bands, Refl, SensorBands, ExpressIndex, BPvars, ReflFactor = 1, p=NULL){

  Bands_SI <- lapply(seq_len(nrow(Bands)), function(i) Bands[i,])
  SIopt <- pblapply(Bands_SI,
                    FUN = compute_SI_fromExp_Corr,
                    Refl = Refl,
                    SensorBands = SensorBands,
                    ExpressIndex = ExpressIndex,
                    BPvars = BPvars,
                    ReflFactor = ReflFactor, NameIndex = 'SI')

  SIopt <- do.call(rbind,SIopt)
  if (!is.null(p)) p()
  return(SIopt)
}

#' This function computes the correlation between a set of vegetation properties
#' and a set of combinations of spectral bands corresponding to a given type of spectral index
#'
#' @param Bands numeric.
#' @param Refl numeric. reflectance matrix with ncol = nbands and nrows = nsamples
#' @param SensorBands numeric. vector containing central wavelength for each spectral band of ReflMat
#' @param BPvars numeric. Biophysical properties as matrix with ncol = nBPvars and nrows = nsamples
#' @param p function.
#'
#' @return SIopt list. includes band Combinations and corresponding correlation between biophysical variables and SI
#' @importFrom pbapply pblapply
#' @export

sub_SI_CR_Corr <- function(Bands, Refl, SensorBands, BPvars, p=NULL){
  Bands_SI <- lapply(seq_len(nrow(Bands)), function(i) Bands[i,])
  SIopt <- pblapply(Bands_SI,
                    FUN = compute_SI_CR_Corr,
                    Refl = Refl,
                    SensorBands = SensorBands,
                    BPvars = BPvars)
  SIopt <- do.call(rbind,SIopt)
  if (!is.null(p)){p()}
  return(SIopt)
}

# compute_lasrc_indices <- function(lasrc_dir){
#   # writes indices in files
#   require('raster')
#   source('R/functions.R')
#   LASRCBands <- c('band2'=496.6, 'band3'=560.0, 'band4'=664.5, 'band5'=703.9, 'band6'=740.2,
#                   'band7' = 782.5, 'band8' = 835.1, 'band11' = 1613.7, 'band12' = 2202.4, 'band8a' = 864.8)
#
#   bands_lasrc = data.frame(
#     name = c(sprintf('band%d', c(2:8, 11:12)), 'band8a'),
#     # res = c(rep(10, 3), rep(20, 3), 10, rep(20, 3)),
#     s2name = c(sprintf('B%02d', c(2:8, 11:12)), 'B08A'),
#     stringsAsFactors=F)
#
#   r = brick(stack(band_LASRC_file(lasrc_dir, bands_lasrc$name)))
#   indices = ComputeSpectralIndices_Raster(r, LASRCBands)
#   for(i in names(indices$SpectralIndices)){
#     file = file.path(lasrc_dir, paste0(basename(lasrc_dir), '_', i, '.tif'))
#     if(!file.exists(file)){
#       writeRaster(indices$SpectralIndices[[i]], filename = file)
#       cat(sprintf('File written: %s\n', file))
#     }
#   }
# }
