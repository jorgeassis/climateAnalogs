# ----------------------------------------------------------------------
# ----------------------------------------------------------------------
#
#
# ----------------------------------------------------------------------

closeAllConnections()
rm(list=(ls()[ls()!="v"]))
gc(reset=TRUE)
source("mainFunctions.R")

# ----------------------------------------------------------------------
# ----------------------------------------------------------------------
# Simple animation from focal points in map

scenario <- "ssp119" # ssp119 ssp585

dataFolder <- "../Data/"
resultsSubFolder <- paste0("finalRun/",scenario,"/")
resultsFolder <- paste0("../Results/",resultsSubFolder)

load(paste0(resultsFolder,"/dataStructure.RData"))
load(paste0(resultsFolder,"/climaticAnalogs.RData"))
load(file=paste0(dataFolder,"/Spatial/globalLandmass/landmassGlobal.RData"))

focal <- 105000

# ----------

focalPts <- data.frame(dataStructure[focal,c("x","y")])

for( i in 1:ncol(climaticAnalogs)) {
  
  focalPts <- rbind(focalPts,data.frame(dataStructure[dataStructure$cell == climaticAnalogs[focal,i],c("x","y")]))
  
}

landmassRegion <- crop(landmass,extent(min(focalPts$x) - 10,max(focalPts$x) + 10,min(focalPts$y) - 10,max(focalPts$y) + 10))

plot(landmassRegion, col="gray")

for( i in 1:ncol(climaticAnalogs)) {
  
  points(focalPts[i,], pch=16)
  readline(prompt="Press [enter] to continue")
  
}

# ----------------------------------------------------------------------
# ----------------------------------------------------------------------