# ----------------------------------------------------------------------
# ----------------------------------------------------------------------
#
#
# ----------------------------------------------------------------------

closeAllConnections()
rm(list=(ls()[ls()!="v"]))
gc(reset=TRUE)
source("mainFunctions.R")

# library(credentials)
# set_github_pat()

# ----------------------------------------------------------------------
# ----------------------------------------------------------------------

resultsFolder <- "../Results/ssp119/"

climateDataBaselineDir <- "../../../Manuscripts/Global connectivity corridors and refugia of climatic analogs for marine biodiversity/Data/Climate/Baseline/"
climateDataProjectionDir <- "../../../Manuscripts/Global connectivity corridors and refugia of climatic analogs for marine biodiversity/Data/Climate/ssp119/"

load(file="../../../Manuscripts/Global connectivity corridors and refugia of climatic analogs for marine biodiversity/Data/Spatial/EEZ/EEZ.RData")
regionOfInterest <- AfricaRegionEZZ

# -------------------
# -------------------
# Data Structure

if(! dir.exists(paste0(resultsFolder))) { dir.create(paste0(resultsFolder), recursive = T) }
if(! dir.exists(paste0(resultsFolder,"/Raster/"))) { dir.create(paste0(resultsFolder,"/Raster/"), recursive = T) }

dataStructureRaster <- raster(list.files(climateDataBaselineDir,full.names = TRUE)[1])
dataStructureRaster <- crop(dataStructureRaster,regionOfInterest)
dataStructureRaster <- mask(dataStructureRaster,regionOfInterest)
plot(dataStructureRaster)

dataStructure <- Which(!is.na(dataStructureRaster), cells=TRUE, na.rm=TRUE)
dataStructure <- data.frame(cell=dataStructure,row=rowFromCell(dataStructureRaster,dataStructure),column=colFromCell(dataStructureRaster,dataStructure),xyFromCell(dataStructureRaster,dataStructure))
head(dataStructure)
save(dataStructure, file=paste0(resultsFolder,"/dataStructure.RData"))

# -------------------
# Climate Data

temperatureBaselineConditionsFiles <- list.files(climateDataBaselineDir, full.names = TRUE, pattern = "OceanTemperature BenthicDepthMean Max")
oxygenBaselineConditionsFiles <- list.files(climateDataBaselineDir, full.names = TRUE, pattern = "DissolvedMolecularOxygen BenthicDepthMean Mean")

available.t <- min(length(temperatureBaselineConditionsFiles),length(oxygenBaselineConditionsFiles))

temperatureBaselineConditions <- array(NA, c(dim(dataStructureRaster)[1],dim(dataStructureRaster)[2],available.t))
oxygenBaselineConditions <- array(NA, c(dim(dataStructureRaster)[1],dim(dataStructureRaster)[2],available.t))

for( t in 1:available.t) {
  
  dataRaster <- raster(temperatureBaselineConditionsFiles[t])
  dataRaster <- crop(dataRaster,regionOfInterest)
  dataRaster <- mask(dataRaster,regionOfInterest)
  temperatureBaselineConditions[,,t] <- as.matrix(dataRaster)
  
  dataRaster <- raster(oxygenBaselineConditionsFiles[t])
  dataRaster <- crop(dataRaster,regionOfInterest)
  dataRaster <- mask(dataRaster,regionOfInterest)
  oxygenBaselineConditions[,,t] <- as.matrix(dataRaster)
  
}

# -----

temperatureProjectionConditionsFiles <- list.files(climateDataProjectionDir, full.names = TRUE, pattern = "OceanTemperature BenthicDepthMean")
oxygenProjectionConditionsFiles <- list.files(climateDataProjectionDir, full.names = TRUE, pattern = "DissolvedMolecularOxygen BenthicDepthMean")

available.t <- min(length(temperatureProjectionConditionsFiles),length(oxygenProjectionConditionsFiles))

temperatureProjectionConditions <- array(NA, c(dim(dataStructureRaster)[1],dim(dataStructureRaster)[2],available.t))
oxygenProjectionConditions <- array(NA, c(dim(dataStructureRaster)[1],dim(dataStructureRaster)[2],available.t))

for( t in 1:available.t) {
  
  dataRaster <- raster(temperatureProjectionConditionsFiles[t])
  dataRaster <- crop(dataRaster,regionOfInterest)
  dataRaster <- mask(dataRaster,regionOfInterest)
  temperatureProjectionConditions[,,t] <- as.matrix(dataRaster)
  
  dataRaster <- raster(oxygenProjectionConditionsFiles[t])
  dataRaster <- crop(dataRaster,regionOfInterest)
  dataRaster <- mask(dataRaster,regionOfInterest)
  oxygenProjectionConditions[,,t] <- as.matrix(dataRaster)
  
}

# -------------------
# -------------------

movingWindow <- 5 # 0.25 * 5 around focal

BaselineT <- dim(temperatureBaselineConditions)[3]
FutureT <- dim(temperatureProjectionConditions)[3]

# -------

file.remove(list.files(resultsFolder, pattern="climaticAnalogs", full.names=T))
climaticAnalogs.bm <- big.matrix(nrow=nrow(dataStructure),ncol=FutureT , backingpath=resultsFolder , backingfile = "climaticAnalogs.bin", descriptorfile = "climaticAnalogs.desc")
climaticAnalogs.bm.desc <- dget( paste0(resultsFolder,"/climaticAnalogs.desc") )
climaticAnalogsDist.bm <- big.matrix(nrow=nrow(dataStructure),ncol=FutureT , backingpath=resultsFolder , backingfile = "climaticAnalogsDist.bin", descriptorfile = "climaticAnalogsDist.desc")
climaticAnalogsDist.bm.desc <- dget( paste0(resultsFolder,"/climaticAnalogsDist.desc") )

# -------

number.Cores <- 16

Cluster <- makeCluster( number.Cores )
registerDoParallel( Cluster ) 

parallelProcess <- foreach(focal=1:nrow(dataStructure), .verbose=FALSE, .packages=c("bigmemory","raster")) %dopar% {
  
  # focal = 1
  # 16.95728 secs
  
  library(FNN)
  library(adehabitatLT)
  
  time.i <- Sys.time()
  
  # -------
  
  climaticAnalogs.bm.f <- attach.big.matrix(climaticAnalogs.bm.desc)
  climaticAnalogsDist.bm.f <- attach.big.matrix(climaticAnalogsDist.bm.desc)
  
  cell <- dataStructure[focal,"cell"]
  row <- dataStructure[focal,"row"]
  column <- dataStructure[focal,"column"]

  analog.cell <- cell

  for( t.future in 1:FutureT) {

    analog.row <- dataStructure[dataStructure$cell == analog.cell,"row"]
    analog.column <- dataStructure[dataStructure$cell == analog.cell,"column"]

    neighbours <- dataStructure[dataStructure$row <= analog.row + movingWindow - 1 & dataStructure$row >= analog.row - movingWindow - 1 &
                                dataStructure$column <= analog.column + movingWindow - 1 & dataStructure$column >= analog.column - movingWindow - 1 , ]
    
    climaticAnalogs.t <- array(NA, dim=c(BaselineT,nrow(neighbours)))

    for(t.baseline in 1:BaselineT) {
      
      pairedClimate <- array(NA, dim=c(nrow(neighbours) + 1 ,2))
      
      # first row is the baseline
      pairedClimate[1,] <- c(temperatureBaselineConditions[row,column,t.baseline],oxygenBaselineConditions[row,column,t.baseline])
      
      for(i in 1:nrow(neighbours)) {
        
        analog.row.i <- neighbours[i,"row"]
        analog.column.i <- neighbours[i,"column"]
        
        pairedClimate[i + 1,] <- c(temperatureProjectionConditions[analog.row.i,analog.column.i,t.future],oxygenProjectionConditions[analog.row.i,analog.column.i,t.future])
        
      }
  
      # ---------
      
      # Standardization [0-1]

      # pairedClimate[,1] <- rangeVar(pairedClimate[,1])
      # pairedClimate[,2] <- rangeVar(pairedClimate[,2])
      
      # Standardization [divide by standard deviation]
      
      pairedClimate[,1] <- pairedClimate[,1] / sd(pairedClimate[,1])
      pairedClimate[,2] <- pairedClimate[,2] / sd(pairedClimate[,2])
      
      # ---------
        
      # myDist <- dist(pairedClimate)
      # myDist <- as.matrix(myDist)[1,-1]

      # Euclidean nearest neighbour distance in the z-standardized PCs of interannual climatic variability, i.e. the Mahalanobian nearest neighbour. 
      myDist <- as.numeric(get.knnx(data=matrix(pairedClimate[1,], ncol=2),query=matrix(pairedClimate[-1,], ncol=2),k=1,algorithm="brute")$nn.dist) 
      
      # get.knnx(data=as.matrix(data.frame(i=c(1,2,3,4),j=c(2,2,3,4))),query=as.matrix(data.frame(i=c(1,6),j=c(2,9))),k=1,algorithm="brute")
      
      # percentile of the nearest neighbour distance on the chi distribution with degrees of freedom equaling the dimensionality of the distance measurement (PCs)
      myDist <- pchi(myDist, df = ncol(pairedClimate))

      # myDist <- mahalanobis(pairedClimate, colMeans(pairedClimate), cov(pairedClimate))
      # myDist <- rangeVar(myDist)
      # myDist <- myDist[-1]
      
      climaticAnalogs.t[t.baseline,] <- myDist
      
    }
    
    analog.cell <- neighbours[which.min(apply(climaticAnalogs.t,2,mean)),"cell"]
    
    climaticAnalogs.bm.f[focal,t.future] <- analog.cell
    climaticAnalogsDist.bm.f[focal,t.future] <- spDists( as.matrix(dataStructure[focal,c("x","y")]) , as.matrix(neighbours[which.min(apply(climaticAnalogs.t,2,mean)),c("x","y")]) , longlat = TRUE  )
    
    }
    
  # -------
  
  Sys.time() - time.i
  
  # -------
  
  return(NULL)
  
}
  
stopCluster(Cluster); rm(Cluster)
closeAllConnections()

# -------

climaticAnalogs <- attach.big.matrix(climaticAnalogs.bm.desc)
climaticAnalogsDist <- attach.big.matrix(climaticAnalogsDist.bm.desc)

climaticAnalogs <- as.matrix(climaticAnalogs)
climaticAnalogsDist <- as.matrix(climaticAnalogsDist)

plot(apply(climaticAnalogsDist,2,mean))

sum(dataStructure$cell == climaticAnalogs[,ncol(climaticAnalogs)])

head(dataStructure)
dataStructure[dataStructure$cell == 321,]

climaticAnalogs[5,]
climaticAnalogsDist[5,]

dataStructure[3,]
dataStructure[dataStructure$cell == 1883,]
dataStructure[1883,]

temperatureBaselineConditions[1,62,18]
temperatureProjectionConditions[63,18,1]

oxygenBaselineConditions[1,62,18]
oxygenProjectionConditions[63,18,80]

# --------------------
# --------------------

analogResults <- data.frame(dataStructure,
                            
                            refugia = dataStructure$cell == climaticAnalogs[,ncol(climaticAnalogs)],
                            corridor = rangeVar(sapply(dataStructure$cell, function(x) { sum(climaticAnalogs[,-ncol(climaticAnalogs)] == x) } )),
                            finalAnalogDistance = climaticAnalogsDist[,ncol(climaticAnalogsDist)],
                            meanAnalogDistance = apply(climaticAnalogsDist,1,mean),
                            maxAnalogDistance = apply(climaticAnalogsDist,1,max),
                            
                            finalTemperatureChange = apply(dataStructure[,c("row","column")],1,function(x) { temperatureProjectionConditions[x[[1]],x[[2]],dim(temperatureProjectionConditions)[3]] - temperatureBaselineConditions[x[[1]],x[[2]],dim(temperatureBaselineConditions)[3]] } ),
                            meanTemperatureChange = apply(apply(dataStructure[,c("row","column")],1,function(x) { temperatureProjectionConditions[x[[1]],x[[2]],] - temperatureBaselineConditions[x[[1]],x[[2]],] } ),2,mean),
                            maxTemperatureChange = apply(apply(dataStructure[,c("row","column")],1,function(x) { temperatureProjectionConditions[x[[1]],x[[2]],] - temperatureBaselineConditions[x[[1]],x[[2]],] } ),2,max),

                            finalOxygenChange = apply(dataStructure[,c("row","column")],1,function(x) { oxygenProjectionConditions[x[[1]],x[[2]],dim(oxygenProjectionConditions)[3]] - oxygenBaselineConditions[x[[1]],x[[2]],dim(oxygenBaselineConditions)[3]] } ),
                            meanOxygenChange = apply(apply(dataStructure[,c("row","column")],1,function(x) { oxygenProjectionConditions[x[[1]],x[[2]],] - oxygenBaselineConditions[x[[1]],x[[2]],] } ),2,mean),
                            maxOxygenChange = apply(apply(dataStructure[,c("row","column")],1,function(x) { oxygenProjectionConditions[x[[1]],x[[2]],] - oxygenBaselineConditions[x[[1]],x[[2]],] } ),2,max)

                            )

save(analogResults, file=paste0(resultsFolder,"/analogResults.RData"))

# to Raster

dataStructureRaster[] <- NA

refugia <- dataStructureRaster
refugia[analogResults$cell] <- analogResults$refugia
writeRaster(refugia,file=paste0(resultsFolder,"/Raster/refugia.tif"), format="GTiff", overwrite=TRUE)

corridor <- dataStructureRaster
corridor[analogResults$cell] <- analogResults$corridor
writeRaster(corridor,file=paste0(resultsFolder,"/Raster/corridor.tif"), format="GTiff", overwrite=TRUE)

# quatile corridor


finalAnalogDistance <- dataStructureRaster
finalAnalogDistance[analogResults$cell] <- analogResults$finalAnalogDistance
writeRaster(finalAnalogDistance,file=paste0(resultsFolder,"/Raster/finalAnalogDistance.tif"), format="GTiff", overwrite=TRUE)

meanAnalogDistance <- dataStructureRaster
meanAnalogDistance[analogResults$cell] <- analogResults$meanAnalogDistance
writeRaster(meanAnalogDistance,file=paste0(resultsFolder,"/Raster/meanAnalogDistance.tif"), format="GTiff", overwrite=TRUE)

maxAnalogDistance <- dataStructureRaster
maxAnalogDistance[analogResults$cell] <- analogResults$maxAnalogDistance
writeRaster(maxAnalogDistance,file=paste0(resultsFolder,"/Raster/maxAnalogDistance.tif"), format="GTiff", overwrite=TRUE)

finalTemperatureChange <- dataStructureRaster
finalTemperatureChange[analogResults$cell] <- analogResults$finalTemperatureChange
writeRaster(finalTemperatureChange,file=paste0(resultsFolder,"/Raster/finalTemperatureChange.tif"), format="GTiff", overwrite=TRUE)

meanTemperatureChange <- dataStructureRaster
meanTemperatureChange[analogResults$cell] <- analogResults$meanTemperatureChange
writeRaster(meanTemperatureChange,file=paste0(resultsFolder,"/Raster/meanTemperatureChange.tif"), format="GTiff", overwrite=TRUE)

maxTemperatureChange <- dataStructureRaster
maxTemperatureChange[analogResults$cell] <- analogResults$maxTemperatureChange
writeRaster(maxTemperatureChange,file=paste0(resultsFolder,"/Raster/maxTemperatureChange.tif"), format="GTiff", overwrite=TRUE)

finalOxygenChange <- dataStructureRaster
finalOxygenChange[analogResults$cell] <- analogResults$finalOxygenChange
writeRaster(finalOxygenChange,file=paste0(resultsFolder,"/Raster/finalOxygenChange.tif"), format="GTiff", overwrite=TRUE)

meanOxygenChange <- dataStructureRaster
meanOxygenChange[analogResults$cell] <- analogResults$meanOxygenChange
writeRaster(meanOxygenChange,file=paste0(resultsFolder,"/Raster/meanOxygenChange.tif"), format="GTiff", overwrite=TRUE)

maxOxygenChange <- dataStructureRaster
maxOxygenChange[analogResults$cell] <- analogResults$maxOxygenChange
writeRaster(maxOxygenChange,file=paste0(resultsFolder,"/Raster/maxOxygenChange.tif"), format="GTiff", overwrite=TRUE)

# ----------------------------------------------------------------------
# ----------------------------------------------------------------------