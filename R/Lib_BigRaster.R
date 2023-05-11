#' this function aims at computing spectral indices from big raster data
#' it computes the spectral indices based on clsest bands with Sentinel-2
#' and assumes that the bands of the S2 data follow this order
#' wavelength	= {496.6, 560.0, 664.5, 703.9, 740.2, 782.5, 835.1, 864.8, 1613.7, 2202.4}
#' Full description of the indices can be found here:
#' https://www.sentinel-hub.com/eotaxonomy/indices
#'
#' @param raster_path character. path for the big raster
#' @param PathOut character. path for output directory
#' @param MaskRaster character. path for mask file corresponding to raster_path
#' @param SensorBands numeric. wavelength in nanometers of the spectral bands of raster_path
#' @param Sel_Indices  list. list of spectral indices to be computed
#' @param ReflFactor numeric. multiplying factor used to write reflectance in image (==10000 for S2)
#' @param Offset numeric. offset (when Refl between 0 and 1) to be applied on reflectance. Useful to avoid zero values
#' @param S2Bands numeric. wavelength of the spectral bands corresponding to S2 (default = S2A)
#'
#' @return list. list of rasters corresponding to spectral indices
#' @importFrom methods is
#' @importFrom progress progress_bar
#' @importFrom raster raster brick blockSize readStart readStop getValues writeStart writeStop writeValues
#' @export

compute_S2SI_BigRaster <- function(raster_path, PathOut, MaskRaster = FALSE,
                                   SensorBands, Sel_Indices = 'ALL',
                                   ReflFactor = 1, Offset = 0,
                                   S2Bands = data.frame('B2'=492.7, 'B3'=559.8, 'B4'=664.6,
                                                        'B5'=704.1, 'B6'=740.5, 'B7' = 782.8,
                                                        'B8' = 832.8, 'B8A' = 864.7,
                                                        'B11' = 1613.7, 'B12' = 2202.4)){
  SpectralIndices <- list()
  Sen2S2 <- get_closest_bands(SensorBands,S2Bands)
  # read by chunk to avoid memory problem
  blk <- blockSize(brick(raster_path))
  # reflectance file
  r_in <- readStart(brick(raster_path))
  # mask file
  r_inmask <- FALSE
  if (!MaskRaster==FALSE){
    if (file.exists(MaskRaster)){
      r_inmask <- readStart(raster(MaskRaster))
    } else if (!file.exists(MaskRaster)){
      message('WARNING: Mask file does not exist:')
      print(MaskRaster)
      message('Processing all image')
    }
  }
  listIndices <- list('ARI1','ARI2','ARVI','BAI','BAIS2','CCCI','CHL_RE','CRI1','CRI2','EVI','EVI2',
                      'GRVI1','GNDVI','IRECI','LAI_SAVI','MCARI','mNDVI705','MSAVI2',
                      'MSI','mSR705','MTCI','nBR_RAW','NDI_45','NDII','NDSI','NDVI','NDVI_G',
                      'NDVI705','NDWI1','NDWI2','PSRI','PSRI_NIR','RE_NDVI','RE_NDWI','S2REP',
                      'SAVI','SIPI','SR','CR_SWIR','CR_RE')
  if (Sel_Indices[1]=='ALL'){
    Sel_Indices <- listIndices
  }
  # initiate progress bar
  pgbarlength <- blk$n
  pb <- progress_bar$new(
    format = "Computing SI [:bar] :percent in :elapsedfull , estimated time remaining :eta",
    total = pgbarlength, clear = FALSE, width= 100)
  # output files
  SIpath <- r_out <- list()
  for (SI in Sel_Indices){
    SIpath[[SI]] <- file.path(PathOut,paste(basename(raster_path),SI,sep = '_'))
    r_out[[SI]] <- writeStart(raster(raster_path), filename = SIpath[[SI]],format = "ENVI", overwrite = TRUE)
  }
  SIpath <- lapply(SIpath,FUN = paste, '.envi',sep = '')
  # loop over blocks
  for (i in seq_along(blk$row)) {
    # read values for block
    # format is a matrix with rows the cells values and columns the layers
    BlockVal <- getValues(r_in, row = blk$row[i], nrows = blk$nrows[i])
    FullLength <- dim(BlockVal)[1]

    if (typeof(r_inmask)=='logical'){
      # automatically filter pixels corresponding to negative values
      SelectPixels <- which(BlockVal[,1]>0)
      BlockVal <- BlockVal[SelectPixels,Sen2S2]
    } else if (typeof(r_inmask)=='S4'){
      MaskVal <- getValues(r_inmask, row = blk$row[i], nrows = blk$nrows[i])
      SelectPixels <- which(MaskVal ==1)
      BlockVal <- BlockVal[SelectPixels,Sen2S2]
    }
    names(BlockVal) <- names(Sen2S2)

    for (SI in Sel_Indices){
      SIval <- NA*vector(length = FullLength)
      if (length(SelectPixels)>0){
        BlockVal <- BlockVal/ReflFactor
        # compute spectral index
        SIvaltmp <- compute_S2SI_from_Sensor(Refl = BlockVal,
                                             SensorBands = S2Bands,
                                             Sel_Indices = SI)
        SIval[SelectPixels] <- c(SIvaltmp$SpectralIndices[[SI]])
      }
      r_out[[SI]] <- writeValues(r_out[[SI]], SIval, blk$row[i],format = "ENVI", overwrite = TRUE)
    }
    pb$tick()
  }
  # close files
  r_in <- readStop(r_in)
  if (typeof(r_inmask)=='S4'){
    r_inmask <- readStop(r_inmask)
  }
  for (SI in Sel_Indices){
    r_out[[SI]] <- writeStop(r_out[[SI]])
    # write biophysical variable name in headers
    HDR <- read_ENVI_header(get_HDR_name(SIpath[[SI]]))
    HDR$`band names` <- paste('{',SI,'}',sep = '')
    write_ENVI_header(HDR, get_HDR_name(SIpath[[SI]]))
  }
  gc()
  return(SIpath)
}

#' this function aims at stacking individual rasters into a unique raster file
#' when the raster size is important.
#' the file is then written in ENVI BIL as it is written line per line
#'
#' @param list_rasters list. list of raster files to be stacked
#' @param stack_file character. path for output file to be stacked
#' @param MaskRaster character. path for mask file corresponding to raster_path
#' @param names_rasters character. name for each of the rasters in list_rasters, if not already named
#'
#' @return list. list of rasters corresponding to spectral indices
#' @importFrom methods is
#' @importFrom utils read.table
#' @importFrom raster raster brick blockSize readStart readStop getValues writeStart writeStop writeValues
#' @export

stack_BigRaster <- function(list_rasters, stack_file, MaskRaster = FALSE,
                            names_rasters = FALSE){

  # define number of rasters to stack
  nbBands <- length(list_rasters)
  # read list_rasters line by line
  # blk <- blockSize(brick(list_rasters[[1]]),chunksize = dim(raster(list_rasters[[1]]))[2])
  blk <- blockSize(brick(list_rasters[[1]]))
  # initiate reading each file
  r_in <- list()
  if (names_rasters == FALSE) names_rasters <- names(list_rasters)
  names(list_rasters) <- names_rasters
  for (band in names_rasters){
    r_in[[band]] <- readStart(raster(list_rasters[[band]]))
  }
  # mask file
  r_inmask <- FALSE
  if (!MaskRaster==FALSE){
    if (file.exists(MaskRaster)){
      r_inmask <- readStart(raster(MaskRaster))
    } else if (!file.exists(MaskRaster)){
      message('WARNING: Mask file does not exist:')
      print(MaskRaster)
      message('Processing all image')
    }
  }
  # initiate progress bar
  pgbarlength <- blk$n
  pb <- progress_bar$new(
    format = "Stacking rasters [:bar] :percent in :elapsedfull , estimated time remaining :eta",
    total = pgbarlength, clear = FALSE, width= 100)
  # output files
  r_out <- writeStart(brick(stack(list_rasters)), filename = stack_file,format = "EHdr", overwrite = TRUE)
  # loop over blocks
  for (i in seq_along(blk$row)) {
    # read values for block
    # format is a matrix with rows the cells values and columns the layers
    BlockVal <- list()
    for (band in names_rasters){
      BlockVal[[band]] <- getValues(r_in[[band]], row = blk$row[i], nrows = blk$nrows[i])
      FullLength <- dim(BlockVal[[band]])[1]
      if (typeof(r_inmask)=='S4'){
        MaskVal <- getValues(r_inmask, row = blk$row[i], nrows = blk$nrows[i])
        elim <- which(MaskVal ==0)
        BlockVal[[band]][elim] <- NA
      }
    }
    r_out <- writeValues(r_out, do.call('cbind',BlockVal), blk$row[i],format = "EHdr", overwrite = TRUE)
    pb$tick()
  }
  # close files
  for (band in names_rasters){
    r_in[[band]] <- readStop(r_in[[band]])
  }
  if (typeof(r_inmask)=='S4'){
    r_inmask <- readStop(r_inmask)
  }
  r_out <- writeStop(r_out)
  # stack_file <- paste(stack_file, '.envi',sep = '')
  # convert HDR to ENVI format
  raster::hdr(raster(stack_file), format = "ENVI")
  HDR <- read_ENVI_header(get_HDR_name(stack_file))
  HDR$`band names` <- paste('',names_rasters,'',sep = '')
  HDR$`coordinate system string` <- read.table(paste(file_path_sans_ext(stack_file), ".prj", sep = ""))
  HDR$`z plot range` <- NULL
  write_ENVI_header(HDR, get_HDR_name(stack_file))
  gc()
  return(stack_file)
}



#' Reads ENVI hdr file
#'
#' @param HDRpath Path of the hdr file
#'
#' @return list of the content of the hdr file
#' @export
read_ENVI_header <- function(HDRpath) {
  # header <- paste(header, collapse = "\n")
  if (!grepl(".hdr$", HDRpath)) {
    stop("File extension should be .hdr")
  }
  HDR <- readLines(HDRpath)
  ## check ENVI at beginning of file
  if (!grepl("ENVI", HDR[1])) {
    stop("Not an ENVI header (ENVI keyword missing)")
  } else {
    HDR <- HDR [-1]
  }
  ## remove curly braces and put multi-line key-value-pairs into one line
  HDR <- gsub("\\{([^}]*)\\}", "\\1", HDR)
  l <- grep("\\{", HDR)
  r <- grep("\\}", HDR)

  if (length(l) != length(r)) {
    stop("Error matching curly braces in header (differing numbers).")
  }

  if (any(r <= l)) {
    stop("Mismatch of curly braces in header.")
  }

  HDR[l] <- sub("\\{", "", HDR[l])
  HDR[r] <- sub("\\}", "", HDR[r])

  for (i in rev(seq_along(l))) {
    HDR <- c(
      HDR [seq_len(l [i] - 1)],
      paste(HDR [l [i]:r [i]], collapse = "\n"),
      HDR [-seq_len(r [i])]
    )
  }

  ## split key = value constructs into list with keys as names
  HDR <- sapply(HDR, split_line, "=", USE.NAMES = FALSE)
  names(HDR) <- tolower(names(HDR))

  ## process numeric values
  tmp <- names(HDR) %in% c(
    "samples", "lines", "bands", "header offset", "data type",
    "byte order", "default bands", "data ignore value",
    "wavelength", "fwhm", "data gain values"
  )
  HDR [tmp] <- lapply(HDR [tmp], function(x) {
    as.numeric(unlist(strsplit(x, ",")))
  })

  return(HDR)
}


#' get hdr name from image file name, assuming it is BIL format
#'
#' @param ImPath path of the image
#'
#' @return corresponding hdr
#' @import tools
#' @export
get_HDR_name <- function(ImPath) {
  if (file_ext(ImPath) == "") {
    ImPathHDR <- paste(ImPath, ".hdr", sep = "")
  } else if (file_ext(ImPath) == "bil") {
    ImPathHDR <- gsub(".bil", ".hdr", ImPath)
  } else if (file_ext(ImPath) == "zip") {
    ImPathHDR <- gsub(".zip", ".hdr", ImPath)
  } else {
    ImPathHDR <- paste(file_path_sans_ext(ImPath), ".hdr", sep = "")
  }

  if (!file.exists(ImPathHDR)) {
    message("WARNING : COULD NOT FIND HDR FILE")
    print(ImPathHDR)
    message("Process may stop")
  }
  return(ImPathHDR)
}


#' writes ENVI hdr file
#'
#' @param HDR content to be written
#' @param HDRpath Path of the hdr file
#'
#' @return None
#' @importFrom stringr str_count
#' @export
write_ENVI_header <- function(HDR, HDRpath) {
  h <- lapply(HDR, function(x) {
    if (length(x) > 1 || (is.character(x) && str_count(x, "\\w+") > 1)) {
      x <- paste0("{", paste(x, collapse = ","), "}")
    }
    # convert last numerics
    x <- as.character(x)
  })
  writeLines(c("ENVI", paste(names(HDR), h, sep = " = ")), con = HDRpath)
  return(invisible())
}

#' ENVI functions
#'
#' based on https://github.com/cran/hyperSpec/blob/master/R/read.ENVI.R
#' added wavelength, fwhm, ... to header reading
#' Title
#'
#' @param x character.
#' @param separator character
#' @param trim.blank boolean.
#'
#' @return list.
#' @export
split_line <- function(x, separator, trim.blank = TRUE) {
  tmp <- regexpr(separator, x)
  key <- substr(x, 1, tmp - 1)
  value <- substr(x, tmp + 1, nchar(x))
  if (trim.blank) {
    blank.pattern <- "^[[:blank:]]*([^[:blank:]]+.*[^[:blank:]]+)[[:blank:]]*$"
    key <- sub(blank.pattern, "\\1", key)
    value <- sub(blank.pattern, "\\1", value)
  }
  value <- as.list(value)
  names(value) <- key
  return(value)
}

