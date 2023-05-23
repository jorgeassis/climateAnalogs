
packages.to.use <- c(  
  "adehabitatLT",
  "FNN",
  "credentials",
  "maptools" ,
  "rnaturalearth",
  "rgeos",
  "data.table",
  "ggplot2",
  "foreach",
  "doParallel",
  "parallel",
  "bigmemory",
  "raster",
  "modEvA",
  "boot",
  "ecodist"
  , "gdata"
  , "leaflet"
  , "dismo"
  , "gbm"
  ,"leaflet.extras"
  , "biganalytics"
  , "vegan"
  , "rgdal"
  , "sdmpredictors"
  , "maptools"
  , "sf"
  , "gstat")

packages.to.use <- unique(packages.to.use)

for(package in packages.to.use) {
  print(package)
  if( ! package %in% rownames(installed.packages()) ) { install.packages( package ) }
  if( ! package %in% rownames(installed.packages()) ) { install.packages( package , type = "source", dependencies=TRUE) }
  if( ! package %in% rownames(installed.packages()) ) { stop("Error on package instalation") }
  library(package, character.only = TRUE)
}

## -----------------------

firstToupper <- function(string) {
  return(paste(toupper(substr(string, 1, 1)), substr(string, 2, nchar(string)), sep=""))
}

## -----------------------

loadRData <- function(fileName){
  #loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])
}

## -----------------------

rangeVar <- function(x){(x-min(x,na.rm=T))/(max(x,na.rm=T)-min(x,na.rm=T))}
