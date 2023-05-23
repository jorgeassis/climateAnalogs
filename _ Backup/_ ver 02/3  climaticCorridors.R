

OLD // Review based on climaticAnalogs


# ---------
# ----------------------------------------------------------------------
# ----------------------------------------------------------------------
#
#
# ----------------------------------------------------------------------

setwd("~/Projects/climaticAnalogs/Code")

closeAllConnections()
rm(list=(ls()[ls()!="v"]))
gc(reset=TRUE)
source("mainFunctions.R")

nCores <- 32

# library(credentials)
# set_github_pat()

# Main citation 1: Global Change Biology (2017) 23, 3934–3955, doi: 10.1111/gcb.13645
# Main citation 2: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8390509/

# ----------------------------------------------------------------------
# ----------------------------------------------------------------------

scenario <- "ssp119" # ssp119 ssp585

dataFolder <- "../Data/"
climateDataBaselineDir <- paste0(dataFolder,"/Climate/Baseline/")
climateDataProjectionDir <- paste0(dataFolder,"/Climate/",scenario,"/")

resultsSubFolder <- paste0("finalRun/",scenario,"/")
resultsFolder <- paste0("../Results/",resultsSubFolder)
# file.remove(list.files(resultsFolder, full.names=T, recursive=T))

# ---------

predictors <- c("DissolvedMolecularOxygen Surface Mean",
                "OceanTemperature Surface Max",
                "Silicate Surface Mean",
                "pH Surface Mean",
                "Salinity Surface Mean",
                "TotalPrimaryProductionPhyto Surface Mean")

regionOfInterest <- loadRData(paste0(dataFolder,"/Spatial/EEZ/EEZGlobal.RData"))

# -------------------

bathymetryReclass <- NULL

# bathymetryReclassValF <- 0
# bathymetryReclassValT <- -200
# bathymetryReclass <- raster("../Data/BathymetryDepthMinRes025.tif")
# bathymetryReclass[bathymetryReclass < bathymetryReclassValT] <- NA
# bathymetryReclass[bathymetryReclass > bathymetryReclassValF] <- NA
# bathymetryReclass[bathymetryReclass <= bathymetryReclassValF & bathymetryReclass >= bathymetryReclassValT] <- 1

# -------------------
# -------------------

if(! dir.exists(paste0(resultsFolder))) { dir.create(paste0(resultsFolder), recursive = T) }
if(! dir.exists(paste0(resultsFolder,"/Raster/"))) { dir.create(paste0(resultsFolder,"/Raster/"), recursive = T) }

# -------------------
# -------------------
# Data Structure

# regionOfInterest <- regionOfInterest[12,]

regionOfInterestBuffedVal <- 5

if( class(regionOfInterest) == "SpatialPolygonsDataFrame" ) {
  
  regionOfInterest <- st_as_sf(regionOfInterest)
  sf::sf_use_s2(FALSE)
  regionOfInterestBuffered <- st_buffer( regionOfInterest, regionOfInterestBuffedVal)
  
}

dataStructureRaster <- raster(list.files(climateDataBaselineDir,full.names = TRUE)[1])
dataStructureRasterValues <- dataStructureRaster
dataStructureRaster[!is.na(dataStructureRaster)] <- 1

if( ! is.null(bathymetryReclass) ) { dataStructureRaster <- mask(dataStructureRaster,bathymetryReclass) }

dataStructureRaster <- crop(dataStructureRaster,regionOfInterestBuffered)
dataStructureRaster <- mask(dataStructureRaster,regionOfInterestBuffered)
dataStructureRasterROI <- mask(dataStructureRaster,as_Spatial(st_collection_extract(regionOfInterest)))

plot(dataStructureRaster)
plot(dataStructureRasterROI)

dataStructureROI <- intersect( Which(!is.na(dataStructureRasterROI), cells=TRUE ) , Which(!is.na(dataStructureRasterValues), cells=TRUE ) )
dataStructure <- intersect( Which(!is.na(dataStructureRaster), cells=TRUE ) , Which(!is.na(dataStructureRasterValues), cells=TRUE ) )

dataStructure <- data.frame(cell=dataStructure,
                            row=rowFromCell(dataStructureRaster,dataStructure),
                            column=colFromCell(dataStructureRaster,dataStructure),
                            xyFromCell(dataStructureRaster,dataStructure),
                            roi=as.numeric(sapply(dataStructure,function(x) { x %in% dataStructureROI })))

head(dataStructure)
sum(dataStructure$roi == 1)
nrow(dataStructure)

# -------------------------------------------------
# -------------------------------------------------
# Climate Data

climaticBaselineConditionsFiles <- list.files(climateDataBaselineDir, full.names = TRUE)
available.t <- sapply(predictors,function(x) { which(grepl(x,climaticBaselineConditionsFiles)) })
available.t <- min(sapply(1:length(available.t), function(x) { length(available.t[[x]]) }))
climaticBaselineConditionsFiles <- climaticBaselineConditionsFiles[unlist(sapply(predictors,function(x) { which(grepl(x,climaticBaselineConditionsFiles))[1:available.t] }))]

climaticBaselineConditions <- array(NA, c(dim(dataStructureRaster)[1],dim(dataStructureRaster)[2],length(predictors),available.t))

for( pred in 1:length(predictors)) {
  
  climaticBaselineConditionsFiles.p <- climaticBaselineConditionsFiles[grep(predictors[pred],climaticBaselineConditionsFiles)]
  
  for( t in 1:available.t) {
    
    dataRaster <- raster(climaticBaselineConditionsFiles.p[t])
    dataRaster <- crop(dataRaster,dataStructureRaster)
    dataRaster <- mask(dataRaster,dataStructureRaster)
    
    climaticBaselineConditions[,,pred,t] <- raster::as.matrix(dataRaster)
    
  }
}

# -----
# Test object

plot(raster(climaticBaselineConditions[,,pred,t]))
getNaValues <- apply(dataStructure[dataStructure$roi == 1 , c("row","column")],1,function(x) is.na(climaticBaselineConditions[x[1],x[2],1,1]) )
which(getNaValues)
getNaValues <- apply(dataStructure[dataStructure$roi == 0 , c("row","column")],1,function(x) is.na(climaticBaselineConditions[x[1],x[2],1,1]) )
which(getNaValues)

# -----

climaticProjectionConditionsFiles <- list.files(climateDataProjectionDir, full.names = TRUE)
available.t <- sapply(predictors,function(x) { list(which(grepl(x,climaticProjectionConditionsFiles))) })
available.t <- min(sapply(1:length(available.t), function(x) { length(available.t[[x]]) }))
climaticProjectionConditionsFiles <- climaticProjectionConditionsFiles[unlist(sapply(predictors,function(x) { which(grepl(x,climaticProjectionConditionsFiles))[1:available.t] }))]

climaticProjectionConditions <- array(NA, c(dim(dataStructureRaster)[1],dim(dataStructureRaster)[2],length(predictors),available.t))

for( pred in 1:length(predictors)) {
  
  climaticProjectionConditionsFiles.p <- climaticProjectionConditionsFiles[grep(predictors[pred],climaticProjectionConditionsFiles)]
  
  for( t in 1:available.t) {
    
    dataRaster <- raster(climaticProjectionConditionsFiles.p[t])
    dataRaster <- crop(dataRaster,dataStructureRaster)
    dataRaster <- mask(dataRaster,dataStructureRaster)
    
    climaticProjectionConditions[,,pred,t] <- raster::as.matrix(dataRaster)
    
  }
}

plot(raster(climaticProjectionConditions[,,pred,t]))
getNaValues <- apply(dataStructure[dataStructure$roi == 1 , c("row","column")],1,function(x) is.na(climaticProjectionConditions[x[1],x[2],1,1]) )
which(getNaValues)
getNaValues <- apply(dataStructure[dataStructure$roi == 0 , c("row","column")],1,function(x) is.na(climaticProjectionConditions[x[1],x[2],1,1]) )
which(getNaValues)

# -------------------
# -------------------

movingWindow <- 10 # 0.25 * 25 around focal

BaselineT <- dim(climaticBaselineConditions)[4]
FutureT <- dim(climaticProjectionConditions)[4]
FutureTLocal <- (FutureT-9):FutureT
nCells <- sum(dataStructure$roi == 1)

# -------------------
# -------------------

file.remove(list.files(resultsFolder, pattern="\\.desc", full.names=T))
file.remove(list.files(resultsFolder, pattern="\\.bin", full.names=T))

climaticAnalogs.bm <- big.matrix(nrow=nCells,ncol=FutureT , backingpath=resultsFolder , backingfile = "climaticAnalogs.bin", descriptorfile = "climaticAnalogs.desc")
climaticAnalogs.bm.desc <- dget( paste0(resultsFolder,"/climaticAnalogs.desc") )
climaticAnalogsSigma.bm <- big.matrix(nrow=nCells,ncol=FutureT , backingpath=resultsFolder , backingfile = "climaticAnalogsSigma.bin", descriptorfile = "climaticAnalogsSigma.desc")
climaticAnalogsSigma.bm.desc <- dget( paste0(resultsFolder,"/climaticAnalogsSigma.desc") )
climaticAnalogsDist.bm <- big.matrix(nrow=nCells,ncol=FutureT , backingpath=resultsFolder , backingfile = "climaticAnalogsDist.bin", descriptorfile = "climaticAnalogsDist.desc")
climaticAnalogsDist.bm.desc <- dget( paste0(resultsFolder,"/climaticAnalogsDist.desc") )
climaticLocalSigma.bm <- big.matrix(nrow=nCells,ncol=FutureT , backingpath=resultsFolder , backingfile = "climaticLocalSigma.bin", descriptorfile = "climaticLocalSigma.desc")
climaticLocalSigma.bm.desc <- dget( paste0(resultsFolder,"/climaticLocalSigma.desc") )

# -------------------
# -------------------

# predictorsBk <- predictors
predictors <- predictorsBk

startTime <- Sys.time()
gc(reset=TRUE)
Cluster <- makeCluster( nCores - 2 )
registerDoParallel( Cluster ) # nCells

parallelProcess <- foreach(focal.i = 1:nCells , .verbose=FALSE, .packages=c("FNN","bigmemory","raster","adehabitatLT")) %dopar% {
  
  climaticAnalogs.f <- attach.big.matrix(climaticAnalogs.bm.desc)
  climaticAnalogsSigma.f <- attach.big.matrix(climaticAnalogsSigma.bm.desc)
  climaticAnalogsDist.f <- attach.big.matrix(climaticAnalogsDist.bm.desc)
  climaticLocalSigma.f <- attach.big.matrix(climaticLocalSigma.bm.desc)
  
  # -------
  
  dataStructure.i <- dataStructure[which(dataStructure$roi == 1)[focal.i],]
  
  cell <- dataStructure.i$cell
  row <- dataStructure.i$row
  column <- dataStructure.i$column
  
  focalBaselineClimate <- as.data.frame(sapply(1:length(predictors), function(x) { climaticBaselineConditions[row,column,x,] }))
  
  if( dataStructure.i$roi == 0 ) { stop("Error :: 001")}
  if( sum(is.na(focalBaselineClimate)) > 0 ) { stop("Error :: 002")}
  
  # Standardization [around the mean and divide by standard deviation]
  
  focalBaselineClimate.mean <- apply(focalBaselineClimate,2,mean, na.rm=T)
  focalBaselineClimate.sd <- apply(focalBaselineClimate,2,sd, na.rm=T)
  focalBaselineClimate.sd[focalBaselineClimate.sd == 0] <- mean(focalBaselineClimate.sd)
  
  # focalBaselineClimate <- sweep(focalBaselineClimate,MARGIN=2,focalBaselineClimate.mean,'-')
  focalBaselineClimate <- sweep( focalBaselineClimate , MARGIN=2,focalBaselineClimate.sd,'/')
  
  # Principal component truncation rule
  trunc.SDs <- 0.1
  
  # Principal components analysis.
  PCA <- prcomp(focalBaselineClimate)
  PCs <- max(which(unlist(summary(PCA)[1]) > trunc.SDs)) # min(which(summary(PCA)$importance[3,] > 1 - trunc.SDs)) # 
  
  # First interaction starts with window centered in focal point 
  analog.cell <- cell
  
  # -------
  
  for( t.future in 1:FutureT) {
    
    analog.row <- dataStructure[dataStructure$cell == analog.cell,"row"]
    analog.column <- dataStructure[dataStructure$cell == analog.cell,"column"]
    
    futureAnalogs <- dataStructure[dataStructure$row <= analog.row + movingWindow - 1 & dataStructure$row >= analog.row - movingWindow - 1 &
                                     dataStructure$column <= analog.column + movingWindow - 1 & dataStructure$column >= analog.column - movingWindow - 1 , ]
    
    # plot(dataStructure[dataStructure$cell %in% futureAnalogs$cell,c(4,5)])
    # points(dataStructure[dataStructure$cell %in% analog.cell,c(4,5)], col="red")
    
    futureAnalogsClimate <- as.data.frame(matrix(NA,nrow=nrow(futureAnalogs), ncol=length(predictors)))
    
    for( analog.p in 1:length(predictors)) {
      futureAnalogsClimate[,analog.p] <- apply(futureAnalogs[c("row" , "column")],1, function(x) { climaticProjectionConditions[x[1],x[2],analog.p,t.future] })
    }
    if( sum(is.na(futureAnalogsClimate)) > 0 ) { stop("Error :: 003")}
    
    # Standardization [around the mean and divide by standard deviation]
    #futureAnalogsClimate <- sweep(futureAnalogsClimate,MARGIN=2,focalBaselineClimate.mean,'-')
    futureAnalogsClimate <- sweep(futureAnalogsClimate,MARGIN=2,focalBaselineClimate.sd,'/')
    
    # project the pools onto the PCs
    focalBaselineClimate.PC <- predict(PCA,focalBaselineClimate)
    futureAnalogsClimate.PC <- predict(PCA,futureAnalogsClimate)
    
    # express PC scores as standardized anomalies of reference interannual variability 
    focalBaselineClimate.PC.sd <- apply(focalBaselineClimate.PC,2,sd, na.rm=T)
    focalBaselineClimate.PC <- sweep(focalBaselineClimate.PC,MARGIN=2,focalBaselineClimate.PC.sd,'/')
    futureAnalogsClimate.PC <- sweep(futureAnalogsClimate.PC,MARGIN=2,focalBaselineClimate.PC.sd,'/')
    
    # Euclidean nearest neighbour distance in the z-standardized PCs of interannual climatic variability, i.e. the Mahalanobian nearest neighbour. 
    myDist <- get.knnx(data=as.matrix( focalBaselineClimate.PC[,1:PCs] ),query=as.matrix( futureAnalogsClimate.PC[,1:PCs] ),k=1,algorithm="brute")$nn.dist
    
    # percentile of the nearest neighbour distance on the chi distribution with degrees of freedom equaling the dimensionality of the distance measurement (PCs)
    NN.chi <- pchi( as.vector(myDist) , PCs) # rel.tol=.Machine$double.eps^0.8 
    
    # values of the chi percentiles on a standard half-normal distribution (chi distribution with one degree of freedom)
    NN.sigma <- qchi(NN.chi,1)
    
    NN.sigma[is.na(NN.sigma)] <- qchi(1-1e-16,1)
    NN.sigma[NN.sigma >= qchi(1-1e-16,1)] <- qchi(1-1e-16,1)
    
    analog.cell.list <- futureAnalogs[which(NN.sigma < quantile(NN.sigma,0.05)),"cell"]
    analog.cell <- analog.cell.list[which.min(spDistsN1( as.matrix(dataStructure[which(dataStructure$cell %in% analog.cell.list),c("x","y")]) , as.matrix(dataStructure[which(dataStructure$cell == analog.cell),c("x","y")]) , longlat = TRUE  ))]
    
    analog.sigma.val <- NN.sigma[which(futureAnalogs$cell == analog.cell)]
    
    climaticAnalogs.f[ focal.i , t.future] <- analog.cell
    climaticAnalogsSigma.f[ focal.i , t.future] <- analog.sigma.val
    climaticAnalogsDist.f[ focal.i , t.future] <- spDists( as.matrix(dataStructure[which(dataStructure$cell == cell),c("x","y")]) , as.matrix(dataStructure[which(dataStructure$cell == analog.cell),c("x","y")]) , longlat = TRUE  )
    
  }
  
  # -------------------------------------
  # -------------------------------------
  
  # local sigma
  
  for( t.future in 1:FutureT) {
    
    futureAnalogsClimate <- as.data.frame(t(sapply(1:length(predictors), function(x) { climaticProjectionConditions[row,column,x,t.future] })))
    if( sum(is.na(futureAnalogsClimate)) > 0 ) { stop("Error :: 003")}
    
    # Standardization [around the mean and divide by standard deviation]
    # futureAnalogsClimate <- sweep(futureAnalogsClimate,MARGIN=2,focalBaselineClimate.mean,'-')
    futureAnalogsClimate <- sweep(futureAnalogsClimate,MARGIN=2,focalBaselineClimate.sd,'/')
    
    # project the pools onto the PCs
    focalBaselineClimate.PC <- predict(PCA,focalBaselineClimate)
    futureAnalogsClimate.PC <- predict(PCA,futureAnalogsClimate)
    
    # express PC scores as standardized anomalies of reference interannual variability 
    focalBaselineClimate.PC.sd <- apply(focalBaselineClimate.PC,2,sd, na.rm=T)
    focalBaselineClimate.PC <- sweep(focalBaselineClimate.PC,MARGIN=2,focalBaselineClimate.PC.sd,'/')
    futureAnalogsClimate.PC <- sweep(futureAnalogsClimate.PC,MARGIN=2,focalBaselineClimate.PC.sd,'/')
    
    # Euclidean nearest neighbour distance in the z-standardized PCs of interannual climatic variability, i.e. the Mahalanobian nearest neighbour. 
    myDist <- get.knnx(data=as.matrix( focalBaselineClimate.PC[,1:PCs] ),query=matrix( futureAnalogsClimate.PC[,1:PCs] , nrow=1),k=1,algorithm="brute")$nn.dist
    
    # percentile of the nearest neighbour distance on the chi distribution with degrees of freedom equaling the dimensionality of the distance measurement (PCs)
    NN.chi <- pchi( as.vector(myDist) , PCs) # , rel.tol=.Machine$double.eps^0.8
    
    # values of the chi percentiles on a standard half-normal distribution (chi distribution with one degree of freedom)
    NN.sigma <- qchi(NN.chi,1)
    NN.sigma[is.na(NN.sigma)] <- qchi(1-1e-16,1)
    NN.sigma[NN.sigma >= qchi(1-1e-16,1)] <- qchi(1-1e-16,1)
    
    climaticLocalSigma.f[ focal.i , t.future] <- NN.sigma
    
  }
  
  # -------------------------
  # -------------------------
  
  gc(reset=TRUE)
  return(NULL)
  
}

stopCluster(Cluster); rm(Cluster)
closeAllConnections()
Sys.time() - startTime

# -------

dataStructure <- dataStructure[dataStructure$roi == 1,]
save(dataStructure, file=paste0(resultsFolder,"/dataStructure.RData"))
save(dataStructureRaster, file=paste0(resultsFolder,"/dataStructureRaster.RData"))

climaticAnalogsDist <- attach.big.matrix(climaticAnalogsDist.bm.desc)
climaticAnalogsSigma <- attach.big.matrix(climaticAnalogsSigma.bm.desc)
climaticAnalogs <- attach.big.matrix(climaticAnalogs.bm.desc)
climaticLocalSigma <- attach.big.matrix(climaticLocalSigma.bm.desc)

climaticAnalogsSigma <- climaticAnalogsSigma[]
climaticAnalogsDist <- climaticAnalogsDist[]
climaticAnalogs <- climaticAnalogs[]
climaticLocalSigma <- climaticLocalSigma[]

if(sum(climaticAnalogs[,1] == 0)) { stop("Error :: 002")}

save(climaticAnalogs, file=paste0(resultsFolder,"/climaticAnalogs.RData"))
save(climaticAnalogsDist, file=paste0(resultsFolder,"/climaticAnalogsDist.RData"))
save(climaticAnalogsSigma, file=paste0(resultsFolder,"/climaticAnalogsSigma.RData"))
save(climaticLocalSigma, file=paste0(resultsFolder,"/climaticLocalSigma.RData"))

file.remove(list.files(resultsFolder, pattern="\\.desc", full.names=T))
file.remove(list.files(resultsFolder, pattern="\\.bin", full.names=T))

# -----------------------------------------------------
# -----------------------------------------------------
# Export data

firstYearProjection <- 2020
bathymetry <- raster("../Data/Spatial/BathymetryDepthMinRes025.tif")

proporRefugiaF <- function(x) { sum(x == climaticAnalogs[which(dataStructure$cell == x),] & climaticLocalSigma[which(dataStructure$cell == x),] <= 4 ) / ncol(climaticAnalogs) }
absoluteRefugiaF <- function(x) { x == climaticAnalogs[which(dataStructure$cell == x),ncol(climaticAnalogs)] & climaticLocalSigma[which(dataStructure$cell == x),ncol(climaticAnalogs)] < 4 }
tEmergenceF <- function(x) { which( ! climaticLocalSigma[which(dataStructure$cell == x),] < 4)[1] + firstYearProjection - 1 }
corridorF <- function(x) { sum( climaticAnalogs[,-c(1,ncol(climaticAnalogs))] == x ) }

clust <- makeCluster(nCores)
clusterExport(clust, c("climaticAnalogs","dataStructure","climaticLocalSigma","firstYearProjection"))
proporRefugia <- parSapply(clust, dataStructure$cell , proporRefugiaF )
absoluteRefugia <- parSapply(clust, dataStructure$cell , absoluteRefugiaF )
tEmergence <- parSapply(clust, dataStructure$cell , tEmergenceF )
corridor <- rangeVar(parSapply(clust, dataStructure$cell , corridorF ))
stopCluster(clust)
closeAllConnections()
gc(reset=TRUE)

analogResults <- data.frame(dataStructure,
                            xAnalog = unlist( sapply( climaticAnalogs[,ncol(climaticAnalogs)] , function(x) { res <- dataStructure[ dataStructure$cell == x ,c("x")]; ifelse(length(res) != 0 , res, NA) } ) ),
                            yAnalog =  unlist( sapply( climaticAnalogs[,ncol(climaticAnalogs)] , function(x) { res <- dataStructure[ dataStructure$cell == x ,c("y")]; ifelse(length(res) != 0 , res, NA) } ) ),
                            proporRefugia = proporRefugia ,
                            absoluteRefugia = absoluteRefugia ,
                            tEmergence = tEmergence ,
                            corridor = corridor ,
                            localClimaticDissimilarity = climaticLocalSigma[,ncol(climaticAnalogsSigma)] ,
                            analogClimaticDissimilarity = climaticAnalogsSigma[,ncol(climaticAnalogsSigma)] ,
                            analogDistance = climaticAnalogsDist[,ncol(climaticAnalogsDist)]  )

analogResults$analogDistanceVert <- raster::extract(bathymetry,analogResults[,c("x","y")]) - raster::extract(bathymetry,analogResults[,c("xAnalog","yAnalog")])
analogResults$corridorQuantile95 <- analogResults$corridor >= quantile(analogResults$corridor, 0.95)
analogResults$corridorQuantile75 <- analogResults$corridor >= quantile(analogResults$corridor, 0.75)
analogResults$analogLatitudeShiftPoleward <- abs(analogResults$yAnalog) - abs(analogResults$y)

analogResults[analogResults == Inf] <- NA
analogResults[analogResults == -Inf] <- NA
save(analogResults, file=paste0(resultsFolder,"/analogResults.RData"))

# -----------------------------------------------------
# -----------------------------------------------------
# Address relative rate of change between predictors from t1 to tfinal
# https://www.nature.com/articles/s41598-021-94872-4

t.future <- 81
trunc.SDs <- 0.1
predictors.i <- 1:length(predictors)

file.remove(list.files(resultsFolder, pattern="\\.desc", full.names=T))
file.remove(list.files(resultsFolder, pattern="\\.bin", full.names=T))
dissimilarityPredictors.bm <- big.matrix(nrow=nrow(dataStructure),ncol=length(predictors)+1 , backingpath=resultsFolder , backingfile = "dissimilarityPredictors.bin", descriptorfile = "dissimilarityPredictors.desc")
dissimilarityPredictors.bm.desc <- dget( paste0(resultsFolder,"/dissimilarityPredictors.desc") )

# -----------

startTime <- Sys.time()
gc(reset=TRUE)
Cluster <- makeCluster( nCores )
registerDoParallel( Cluster )

parallelProcess <- foreach(focal.i = 1:nrow(dataStructure) , .verbose=FALSE, .packages=c("FNN","bigmemory","raster","adehabitatLT")) %dopar% {
  
  dissimilarityPredictors.f <- attach.big.matrix(dissimilarityPredictors.bm.desc)
  
  cell <- dataStructure[ focal.i ,"cell"]
  row <- dataStructure[ focal.i ,"row"]
  column <- dataStructure[ focal.i ,"column"]
  
  # Get dissimilarity with all predictors
  
  focalBaselineClimate <- as.data.frame(sapply(predictors.i, function(x) { climaticBaselineConditions[row,column,x,] }))
  focalBaselineClimate.mean <- apply(focalBaselineClimate,2,mean, na.rm=T)
  focalBaselineClimate.sd <- apply(focalBaselineClimate,2,sd, na.rm=T)
  
  # Each observation can be expressed as a conventional standardized anomaly by subtracting the mean and dividing by the standard deviation of the time series
  
  # focalBaselineClimate <- sweep(focalBaselineClimate,MARGIN=2,focalBaselineClimate.mean,'-')
  focalBaselineClimate <- sweep( focalBaselineClimate , MARGIN=2,focalBaselineClimate.sd,'/')
  focalBaselineClimate[focalBaselineClimate == Inf] <- 0
  focalBaselineClimate[is.na(focalBaselineClimate)] <- 0
  
  PCA <- prcomp(focalBaselineClimate)
  PCs <- max(which(unlist(summary(PCA)[1]) > trunc.SDs))
  
  focalfutureClimate <- as.data.frame(t(sapply(predictors.i, function(x) { climaticProjectionConditions[row,column,x,t.future-1] })))
  
  #focalfutureClimate <- sweep(focalfutureClimate,MARGIN=2,focalBaselineClimate.mean,'-')
  focalfutureClimate <- sweep(focalfutureClimate,MARGIN=2,focalBaselineClimate.sd,'/')
  focalfutureClimate[focalfutureClimate == Inf] <- 0
  focalfutureClimate[is.na(focalfutureClimate)] <- 0
  
  focalBaselineClimate.PC <- predict(PCA,focalBaselineClimate)
  futureClimate.PC <- predict(PCA,focalfutureClimate)
  
  focalBaselineClimate.PC.sd <- apply(focalBaselineClimate.PC,2,sd, na.rm=T)
  focalBaselineClimate.PC <- sweep(focalBaselineClimate.PC,MARGIN=2,focalBaselineClimate.PC.sd,'/')
  futureClimate.PC <- sweep(futureClimate.PC,MARGIN=2,focalBaselineClimate.PC.sd,'/')
  
  myDistAll <- get.knnx(data=as.matrix( focalBaselineClimate.PC[,1:PCs] ),query= matrix( futureClimate.PC[,1:PCs], nrow=1 ),k=1,algorithm="brute")$nn.dist
  
  # Get dissimilarity with all but one predictor
  
  results.i <- numeric(length(predictors))
  
  for( int in 1:length(predictors)){
    
    focalBaselineClimate <- as.data.frame(sapply(predictors.i, function(x) { climaticBaselineConditions[row,column,x,] }))
    focalfutureClimate <- as.data.frame(t(sapply(predictors.i, function(x) { climaticProjectionConditions[row,column,x,t.future-1] })))
    
    focalBaselineClimate[,int] <- focalBaselineClimate[1,int]
    focalfutureClimate[1,int] <- focalBaselineClimate[1,int]
    
    # focalBaselineClimate <- sweep(focalBaselineClimate,MARGIN=2,focalBaselineClimate.mean,'-')
    focalBaselineClimate <- sweep( focalBaselineClimate , MARGIN=2,focalBaselineClimate.sd,'/')
    focalBaselineClimate[focalBaselineClimate == Inf] <- 0
    focalBaselineClimate[is.na(focalBaselineClimate)] <- 0
    
    # focalfutureClimate <- sweep(focalfutureClimate,MARGIN=2,focalBaselineClimate.mean,'-')
    focalfutureClimate <- sweep(focalfutureClimate,MARGIN=2,focalBaselineClimate.sd,'/')
    focalfutureClimate[focalfutureClimate == Inf] <- 0
    focalfutureClimate[is.na(focalfutureClimate)] <- 0
    
    focalBaselineClimate.PC <- predict(PCA,focalBaselineClimate)
    futureClimate.PC <- predict(PCA,focalfutureClimate)
    
    focalBaselineClimate.PC.sd <- apply(focalBaselineClimate.PC,2,sd, na.rm=T)
    focalBaselineClimate.PC <- sweep(focalBaselineClimate.PC,MARGIN=2,focalBaselineClimate.PC.sd,'/')
    futureClimate.PC <- sweep(futureClimate.PC,MARGIN=2,focalBaselineClimate.PC.sd,'/')
    
    myDist <- get.knnx(data=as.matrix( focalBaselineClimate.PC[,1:PCs] ),query= matrix( futureClimate.PC[,1:PCs], nrow=1 ),k=1,algorithm="brute")$nn.dist
    
    results.i[int] <- myDistAll - myDist
    
  }
  
  if( Inf %in% results.i) { stop(paste0("! :: ", focal.i))}
  if( -Inf %in% results.i) { stop(paste0("! :: ", focal.i))}
  
  results.i <- ( (results.i - min(results.i) ) / max( results.i - min(results.i) ) ) / sum(  (results.i - min(results.i) ) / max( results.i - min(results.i) )  ) * 100
  
  for( int in 1:length(predictors)){ dissimilarityPredictors.f[focal.i,int+1] <- results.i[int] }
  
  gc(reset=TRUE)
  return(NULL)
  
}

stopCluster(Cluster); rm(Cluster)
closeAllConnections()
Sys.time() - startTime

# -----------

dissimilarityPredictors <- attach.big.matrix(dissimilarityPredictors.bm.desc)
dissimilarityPredictors <- dissimilarityPredictors[]
dissimilarityPredictors[,1] <- dataStructure$cell
colnames(dissimilarityPredictors) <- c("cell",predictors)

save(dissimilarityPredictors, file=paste0(resultsFolder,"/dissimilarityPredictors.RData"))
file.remove(list.files(resultsFolder, pattern="\\.bin", full.names=T))
file.remove(list.files(resultsFolder, pattern="\\.desc", full.names=T))

# -----------------------------------------------------
# -----------------------------------------------------
# Make pairwise matrix of analogs

file.remove(list.files(resultsFolder, pattern="squareMatrix", full.names=T))
squareMatrix.bm <- big.matrix(nrow=length(dataStructure$cell),ncol=length(dataStructure$cell) , backingpath=resultsFolder , backingfile = "squareMatrix.bin", descriptorfile = "squareMatrix.desc")
squareMatrix.bm.desc <- dget( paste0(resultsFolder,"/squareMatrix.desc") )

startTime <- Sys.time()
gc(reset=TRUE)
Cluster <- makeCluster( nCores )
registerDoParallel( Cluster )

parallelProcess <- foreach( i = 1:nrow(climaticAnalogs) , .verbose=FALSE, .packages=c("FNN","bigmemory","raster","adehabitatLT")) %dopar% {
  
  squareMatrix.bm.i <- attach.big.matrix(squareMatrix.bm.desc)
  
  climaticAnalogs.i  <- climaticAnalogs[i,]
  climaticAnalogs.i <- unique(climaticAnalogs[i,]) 
  climaticAnalogs.i <- climaticAnalogs.i[climaticAnalogs.i%in% dataStructure$cell]
  
  for(j in climaticAnalogs.i ) {
    squareMatrix.bm.i[i, which(dataStructure$cell == j) ] <- sum( climaticAnalogs.i == j) / ncol(climaticAnalogs)
  }
  
}

stopCluster(Cluster); rm(Cluster)
closeAllConnections()
Sys.time() - startTime

# -----------

squareMatrix <- attach.big.matrix(squareMatrix.bm.desc)
squareMatrix <- squareMatrix[]
save(squareMatrix, file=paste0(resultsFolder,"/squareMatrix.RData"))
file.remove(list.files(resultsFolder, pattern="\\.desc", full.names=T))
file.remove(list.files(resultsFolder, pattern="\\.bin", full.names=T))

# ----------------------------------------------------------------------
# ----------------------------------------------------------------------