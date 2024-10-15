# spinR v0.11.1
## change
- use GTiff as default raster driver instead of COG

## fix
- fix Refl matrix transposition in compute_S2SI_from_Sensor when nb pixels = nb bands

# spinR v0.11.0
## fix
- moved bigRaster as suggestion
- moved wrapperBig_SI_S2perband from bigRaster to spinR

## addition
- added function listIndices_spinR describing available spectral indices

# spinR v0.10.0
## fix
- fix spectralindices_from_raster when no mask is provided

# spinR v0.9.0
## addition
- added function spectralindices_from_raster working with bigRaster

# spinR v0.8.0
## changes
- removed multiprocess
- removed Lib_bigRaster (check package bigRaster if needed)

# spinR v0.6.0
## changes
- each spectral index is in an individual function

# spinR v0.5.0
## additions
- added Lib_BigRaster


# spinR v0.4.0
## fix
- corrected computation of CR_SWIR for rasters
- lowcase first letter of functions when possible

# spinR v0.3.0
## fix
- changed function names
- added functions to find optimal combination of spectral bands based on an expression and reflectance data corresponding to a set of spectral bands (multiprocessing enabled for large number of combinations)

# spinR v0.2.0
## Fix
- corrected function to compute spectral indices: SpectralIndices instead of spectralindices
