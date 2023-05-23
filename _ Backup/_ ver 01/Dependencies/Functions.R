## ---------------------------------------------------------------------------------------------------------
## ---------------------------------------------------------------------------------------------------------
##
## C.M.I.P. Data Request and Processing
## Muffins 'n' Code
## https://github.com/jorgeassis
##
## ---------------------------------------------------------------------------------------------------------
## ---------------------------------------------------------------------------------------------------------
## Dependencies

packages.to.use <- c("tictoc", "raster", "devtools", "VoCC","rgeos",
                     "rasterVis","gridExtra","doParallel","foreach",
                     "scales","data.table","mapplots","ggplot2","repmis",
                     "sf", "rnaturalearth", "rnaturalearthdata", "viridis",
                     "ggpubr", "REdaS")

packages.to.use <- unique(packages.to.use)

for(package in packages.to.use) {
  print(package)
  if( ! package %in% rownames(installed.packages()) ) { install.packages(package ) }
  if( ! package %in% rownames(installed.packages()) & package == "VoCC" ) { devtools::install_github("JorGarMol/VoCC", dependencies = TRUE) }
  if( ! package %in% rownames(installed.packages()) ) { stop("Error on package instalation") }
  suppressWarnings( library(package, character.only = TRUE) )
}

gplot_data <- function(x, maxpixels = 500000000)  {
  x <- raster::sampleRegular(x, maxpixels, asRaster = TRUE)
  coords <- raster::xyFromCell(x, seq_len(raster::ncell(x)))
  ## Extract values
  dat <- utils::stack(as.data.frame(raster::getValues(x))) 
  names(dat) <- c('value', 'variable')
  
  dat <- dplyr::as.tbl(data.frame(coords, dat))
  
  if (!is.null(levels(x))) {
    dat <- dplyr::left_join(dat, levels(x)[[1]], 
                            by = c("value" = "ID"))
  }
  dat
}
## ---------------------------------------------------------------------------
## ---------------------------------------------------------------------------




