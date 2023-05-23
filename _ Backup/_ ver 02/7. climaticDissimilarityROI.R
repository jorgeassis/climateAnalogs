# ----------------------------------------------------------------------
# ----------------------------------------------------------------------
#
#
# ----------------------------------------------------------------------

# Area estimates of polygon inside polygon are wrong, needs to be determined inside loop.

# ---------

setwd("/Volumes/Jellyfish/Dropbox/Manuscripts/Global connectivity corridors and refugia of climatic analogs for marine biodiversity/Code")

closeAllConnections()
rm(list=(ls()[ls()!="v"]))
gc(reset=TRUE)
source("/Volumes/Jellyfish/Dropbox/theMarineDataScientist/gitRepositories/climaticAnalogs/mainFunctions.R")

# ---------

scenario <- "ssp119" #  ssp585

dataFolder <- "../Data/"
resultsSubFolder <- paste0("finalRun/",scenario,"/")
resultsFolder <- paste0("../Results/",resultsSubFolder)

# ---------

load(paste0(resultsFolder,"/dataStructure.RData"))
load(paste0(resultsFolder,"/dataStructureRaster.RData"))

analogResults <- loadRData(paste0(resultsFolder,"/analogResults.RData"))

# ---------

rasterFiles <- "/Volumes/Jellyfish/Dropbox/Manuscripts/Fish sensitivity to projected climate change/Results/Rasters/"
rasterFiles <- list.files(rasterFiles, pattern="aquamapsRecords.tif",full.names = TRUE, recursive = TRUE)

# ---------

shape <- raster(rasterFiles[1])
shape[] <- NA
cells <- cellFromXY(shape,analogResults[,c("x","y")])
cellsAnalog <- cellFromXY(shape,analogResults[,c("xAnalog","yAnalog")])

# ---------

for( rasterFileName in rasterFiles ) {
  
  rasterROI <- raster(rasterFileName)
  rasterROICells <- Which(rasterROI == 1, cells = TRUE)
  rasterROICells <- rasterROICells[rasterROICells %in% cells]
  
  rasterROICellResults <- analogResults[match(rasterROICells,cells),]
  
  rasterDistributionEEZ <- shape
  rasterDistributionEEZ[rasterROICells] <- 1
  rasterDistributionEEZ[is.na(rasterDistributionEEZ)] <- 0
  
  rasterROISigma <- shape
  rasterROISigma[rasterROICells] <- rasterROICellResults$localClimaticDissimilarity
  writeRaster(rasterROISigma,file=gsub("aquamapsRecords",paste0("climaticDissimilarity_",scenario),rasterFileName), format="GTiff", overwrite=TRUE)
  
  novelClimates <- rasterROISigma
  novelClimates[novelClimates < 2] <- 0
  novelClimates[novelClimates >= 2] <- 1
  writeRaster(rasterROISigma,file=gsub("aquamapsRecords",paste0("novelClimates_",scenario),rasterFileName), format="GTiff", overwrite=TRUE)
  
  novelClimates <- rasterROISigma
  novelClimates[novelClimates < 4] <- 0
  novelClimates[novelClimates >= 4] <- 1
  writeRaster(rasterROISigma,file=gsub("aquamapsRecords",paste0("extremeNovelClimates_",scenario),rasterFileName), format="GTiff", overwrite=TRUE)
  
  rasterAnalog <- rasterize(rasterROICellResults[,c("xAnalog","yAnalog")],shape,field=rasterROICellResults$localClimaticDissimilarity)
  rasterAnalog[!is.na(rasterAnalog)] <- 1
  rasterAnalog[is.na(rasterAnalog)] <- 0
  writeRaster(rasterAnalog,file=gsub("aquamapsRecords",paste0("climaticAnalogs_",scenario),rasterFileName), format="GTiff", overwrite=TRUE)
  
  gainLoss <- calc(stack(rasterDistributionEEZ,rasterAnalog), function(x) { ifelse( x[1] == 1 & x[2] == 1 , 0 , ifelse( x[1] == 0 & x[2] == 1 , 1 , ifelse( x[1] == 1 & x[2] == 0 , -1 , NA ) ) ) } )
  plot(gainLoss, col=c("red","black","green"))
  writeRaster(rasterAnalog,file=gsub("aquamapsRecords",paste0("gainLossClimates_",scenario),rasterFileName), format="GTiff", overwrite=TRUE)
  
}
