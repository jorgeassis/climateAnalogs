# ----------------------------------------------------------------------
# ----------------------------------------------------------------------
#
#
# ----------------------------------------------------------------------

# Main references:

# https://onlinelibrary.wiley.com/doi/full/10.1111/gcb.13645
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8390509/
# https://www.nature.com/articles/s41598-021-94872-4
# https://www.cell.com/one-earth/pdf/S2590-3322%2821%2900602-3.pdf

# ------------------

setwd("~/Projects/climaticAnalogs/Code")
setwd("/Volumes/Jellyfish/Dropbox/Manuscripts/Global connectivity corridors and refugia of climatic analogs for marine biodiversity/Code")

closeAllConnections()
rm(list=(ls()[ls()!="v"]))
gc(reset=TRUE)
source("mainFunctions.R")
nCores <- 30

# library(credentials)
# set_github_pat()

# ----------------------------------------------------------------------
# ----------------------------------------------------------------------

scenario <- "ssp585" # ssp119 ssp585

dataFolder <- "../Data/"
tempFolder <- "../Temp/"
climateDataBaselineDir <- paste0(dataFolder,"/Climate/Baseline/")
climateDataProjectionDir <- paste0(dataFolder,"/Climate/",scenario,"/")

resultsSubFolder <- paste0("GlobalClimaticNovelty/",scenario,"/")
resultsFolder <- paste0("../Results/",resultsSubFolder)
# unlink(resultsFolder, recursive=T)

# -------------------

predictors <- c("DissolvedMolecularOxygen Surface Mean",
                "OceanTemperature Surface Max",
                "pH Surface Mean",
                "TotalPhytoplankton Surface Mean"
) # "Silicate Surface Mean", "TotalPrimaryProductionPhyto Surface Mean" "Salinity Surface Mean",


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
if(! dir.exists(paste0(tempFolder))) { dir.create(paste0(tempFolder), recursive = T) }
if(! dir.exists(paste0(resultsFolder,"/Raster/"))) { dir.create(paste0(resultsFolder,"/Raster/"), recursive = T) }

# -------------------
# -------------------
# Data Structure

dataStructureRaster <- raster(list.files(climateDataBaselineDir,full.names = TRUE)[1])
dataStructureRaster[!is.na(dataStructureRaster)] <- 1
if( ! is.null(bathymetryReclass) ) { dataStructureRaster <- mask(dataStructureRaster,bathymetryReclass) }

plot(dataStructureRaster)

dataStructure <- Which(!is.na(dataStructureRaster), cells=TRUE )
dataStructure <- data.frame(cell=dataStructure,
                            row=rowFromCell(dataStructureRaster,dataStructure),
                            column=colFromCell(dataStructureRaster,dataStructure),
                            xyFromCell(dataStructureRaster,dataStructure))

head(dataStructure)
nrow(dataStructure)
save(dataStructure, file=paste0(resultsFolder,"/dataStructure.RData"))

# -------------------------------------------------
# -------------------------------------------------
# Climate Data

file.remove(list.files(tempFolder, pattern="\\.desc", full.names=T))
file.remove(list.files(tempFolder, pattern="\\.bin", full.names=T))

# -------------

climaticBaselineConditionsFiles <- list.files(climateDataBaselineDir, full.names = TRUE)
available.t <- sapply(predictors,function(x) { which(grepl(x,climaticBaselineConditionsFiles)) })
available.t <- min(sapply(1:length(available.t), function(x) { length(available.t[[x]]) }))
climaticBaselineConditionsFiles <- climaticBaselineConditionsFiles[unlist(sapply(predictors,function(x) { which(grepl(x,climaticBaselineConditionsFiles))[1:available.t] }))]

baselineClimaticVar.CellLoc <- data.table(cell=as.numeric(sapply(dataStructure$cell, function(x) { rep(x, available.t ) } )),id=1:(length(dataStructure$cell) * available.t))
baselineClimaticVar.TimeLoc <- data.table(time=rep(1:available.t,nrow(dataStructure)), id=1:(nrow(dataStructure)*available.t))
setkey(baselineClimaticVar.CellLoc,cell)
setkey(baselineClimaticVar.TimeLoc,time)

baselineClimaticVar.bm <- big.matrix(nrow=nrow(dataStructure)*available.t, ncol=length(predictors) , backingpath=tempFolder , backingfile = "baselineClimaticVarBM.bin", descriptorfile = "baselineClimaticVarBM.desc")
baselineClimaticVar.bm.desc <- dget( paste0(tempFolder,"/baselineClimaticVarBM.desc") )
baselineClimaticVar.bm <- attach.big.matrix(baselineClimaticVar.bm.desc)

for( pred in 1:length(predictors)) {
  
  climaticBaselineConditionsFiles.p <- climaticBaselineConditionsFiles[grep(predictors[pred],climaticBaselineConditionsFiles)]
  
  for( t in 1:available.t) {
    
    dataRaster <- raster(climaticBaselineConditionsFiles.p[t])
    dataRaster <- crop(dataRaster,dataStructureRaster)
    dataRaster <- mask(dataRaster,dataStructureRaster)
    baselineClimaticVar.bm[ baselineClimaticVar.TimeLoc[ time == t, id] ,pred] <- dataRaster[dataStructure$cell]
    
  }
}

rasterPlot <- dataStructureRaster
rasterPlot[] <- NA
rasterPlot[dataStructure$cell] <- baselineClimaticVar.bm[baselineClimaticVar.TimeLoc[ time == t, id],2]
plot(rasterPlot, col = topo.colors(20))
plot(baselineClimaticVar.bm[baselineClimaticVar.CellLoc[ cell == 100000, id],4])
rm(baselineClimaticVar.bm)

# -------------

climaticProjectionConditionsFiles <- list.files(climateDataProjectionDir, full.names = TRUE)
available.t <- sapply(predictors,function(x) { list(which(grepl(x,climaticProjectionConditionsFiles))) })
available.t <- min(sapply(1:length(available.t), function(x) { length(available.t[[x]]) }))
climaticProjectionConditionsFiles <- climaticProjectionConditionsFiles[unlist(sapply(predictors,function(x) { which(grepl(x,climaticProjectionConditionsFiles))[1:available.t] }))]

futureClimaticVar.CellLoc <- data.table(cell=as.numeric(sapply(dataStructure$cell, function(x) { rep(x, available.t ) } )),id=1:(length(dataStructure$cell) * available.t))
futureClimaticVar.TimeLoc <- data.table(time=rep(1:available.t,nrow(dataStructure)), id=1:(nrow(dataStructure)*available.t))
setkey(futureClimaticVar.CellLoc,cell)
setkey(futureClimaticVar.TimeLoc,time)

futureClimaticVar.bm <- big.matrix(nrow=nrow(dataStructure)*available.t, ncol=length(predictors) , backingpath=tempFolder , backingfile = "futureClimaticVarBM.bin", descriptorfile = "futureClimaticVarBM.desc")
futureClimaticVar.bm.desc <- dget( paste0(tempFolder,"/futureClimaticVarBM.desc") )
futureClimaticVar.bm <- attach.big.matrix(futureClimaticVar.bm.desc)

for( pred in 1:length(predictors)) {
  
  climaticProjectionConditionsFiles.p <- climaticProjectionConditionsFiles[grep(predictors[pred],climaticProjectionConditionsFiles)]
  
  for( t in 1:available.t) {
    
    dataRaster <- raster(climaticProjectionConditionsFiles.p[t])
    dataRaster <- crop(dataRaster,dataStructureRaster)
    dataRaster <- mask(dataRaster,dataStructureRaster)
    futureClimaticVar.bm[ futureClimaticVar.TimeLoc[time == t,id] ,pred] <- dataRaster[dataStructure$cell]
    
  }
}

rasterPlot <- dataStructureRaster
rasterPlot[] <- NA
rasterPlot[dataStructure$cell] <- futureClimaticVar.bm[futureClimaticVar.TimeLoc[time==t,id],2]
plot(rasterPlot, col = topo.colors(20))
plot(futureClimaticVar.bm[futureClimaticVar.CellLoc[cell== 200000,id],2])
rm(futureClimaticVar.bm)

# -------------------------------------------------
# -------------------------------------------------

BaselineT <- (length(unique(baselineClimaticVar.TimeLoc[,time])) - 9):(length(unique(baselineClimaticVar.TimeLoc[,time])))
FutureT <- (length(unique(futureClimaticVar.TimeLoc[,time])) - 9):(length(unique(futureClimaticVar.TimeLoc[,time])))

# -------------------

# Global 
ROIRasterType <- "continuous" #  binomial
ROIRasterList <- list.files(climateDataBaselineDir,full.names = TRUE)[1]
ROIRasterNames <- "Global"
PCAType <- "cumulativeImportance" # cumulativeImportance [less aggressive] standardDeviation [more aggressive]

# Region based
# ROIRasterType <- "binomial" # continuous binomial
# ROIRasterList <- list.files("../Data/Rasters/", recursive = TRUE, pattern = "\\.tif", full.names = TRUE)
# ROIRasterList <- ROIRasterList[grepl("Records",ROIRasterList)]
# ROIRasterNames <- list.files("../Data/Rasters/", recursive = TRUE, pattern = "\\.tif", full.names = FALSE)
# ROIRasterNames <- substr(ROIRasterNames,1,unlist(regexpr("/",ROIRasterNames))[1]-1)
# PCAType <- "standardDeviation" # cumulativeImportance standardDeviation

for(ROIRaster.file in ROIRasterList) {
  
  # Prepare results DF
  dataStructureResult <- dataStructure
  
  # Read ROI
  ROIRaster <- raster(ROIRaster.file)
  if( ROIRasterType == "continuous") { ROIRaster <- xyFromCell(ROIRaster,Which(!is.na(ROIRaster), cells=TRUE)) }
  if( ROIRasterType == "binomial") { ROIRaster <- xyFromCell(ROIRaster,Which(ROIRaster == 1, cells=TRUE)) }
  
  Cells <- cellFromXY(dataStructureRaster,ROIRaster)
  Cells <- Cells[ Cells %in% dataStructure$cell ]
  nCells <- length(Cells)
  
  # ------------------------------------------
  # ------------------------------------------
  
  # Get ICV data [historically experienced at the focal stations] and do standardization [divide by standard deviation of baselineClimaticVarSD]
  
  startTime <- Sys.time()
  Cluster <- makeCluster( nCores )
  registerDoParallel( Cluster )
  
  baselineClimaticVarSD <- foreach(cell.i = 1:nCells , .combine = rbind , .verbose=FALSE , .packages=c("bigmemory","data.table") ) %dopar% {
    baselineClimaticVar.bm <- attach.big.matrix(baselineClimaticVar.bm.desc)
    apply(baselineClimaticVar.bm[ baselineClimaticVar.CellLoc[ cell == Cells[cell.i], id],],2,sd)
  }
  row.names(baselineClimaticVarSD) <- NULL
  
  doParallelProcess <- foreach(cell.i = 1:nCells , .verbose=FALSE , .packages=c("bigmemory","data.table") ) %dopar% {
    baselineClimaticVar.bm <- attach.big.matrix(baselineClimaticVar.bm.desc)
    baselineClimaticVar.bm[ baselineClimaticVar.CellLoc[ cell == Cells[cell.i], id],] <- sweep(baselineClimaticVar.bm[ baselineClimaticVar.CellLoc[ cell == Cells[cell.i], id] ,], MARGIN=2, baselineClimaticVarSD[cell.i,],'/')
    return(NULL)
  }
  stopCluster(Cluster); rm(Cluster); closeAllConnections(); gc(reset=TRUE)
  Sys.time() - startTime
  
  # -------
  # -------
  
  # Get baseline and future data and do standardization [divide by standard deviation of interannualClimaticVar]
  
  startTime <- Sys.time()
  Cluster <- makeCluster( nCores )
  registerDoParallel( Cluster )
  focalBaselineClimate <- foreach(cell.i = 1:nCells , .combine = rbind, .verbose=FALSE, .packages=c("bigmemory","data.table")) %dopar% {
    baselineClimaticVar.bm <- attach.big.matrix(baselineClimaticVar.bm.desc)
    apply(baselineClimaticVar.bm[ baselineClimaticVar.CellLoc[ cell == Cells[cell.i], id],][BaselineT,],2,mean)
  }
  futureAnalogsClimate <- foreach(cell.i = 1:nCells , .combine = rbind, .verbose=FALSE, .packages=c("bigmemory","data.table")) %dopar% {
    futureClimaticVar.bm <- attach.big.matrix(futureClimaticVar.bm.desc)
    apply(futureClimaticVar.bm[ futureClimaticVar.CellLoc[ cell == Cells[cell.i], id],][FutureT,],2,mean) / baselineClimaticVarSD[cell.i,]
  }
  stopCluster(Cluster); rm(Cluster); closeAllConnections(); gc(reset=TRUE)
  Sys.time() - startTime
  
  row.names(focalBaselineClimate) <- NULL
  row.names(futureAnalogsClimate) <- NULL
  if( sum(is.na(focalBaselineClimate)) > 0 | sum(is.na(futureAnalogsClimate)) > 0  ) { stop("Error :: 003")}
  
  # ------------------------------------------
  # ------------------------------------------
  
  file.remove(list.files(tempFolder, pattern="ClimateBM", full.names=T))
  
  focalBaselineClimate <- as.big.matrix(focalBaselineClimate , backingpath=tempFolder , backingfile = "focalBaselineClimateBM.bin", descriptorfile = "focalBaselineClimateBM.desc")
  focalBaselineClimate.bm.desc <- dget( paste0(tempFolder,"/focalBaselineClimateBM.desc") )
  futureAnalogsClimate <- as.big.matrix(futureAnalogsClimate , backingpath=tempFolder , backingfile = "futureAnalogsClimateBM.bin", descriptorfile = "futureAnalogsClimateBM.desc")
  futureAnalogsClimate.bm.desc <- dget( paste0(tempFolder,"/futureAnalogsClimateBM.desc") )
  
  startTime <- Sys.time()
  Cluster <- makeCluster( nCores )
  registerDoParallel( Cluster )
  
  sigmaNovelClimate <- foreach(focal.i = 1:nCells , .combine=rbind, .verbose=FALSE, .packages=c("FNN","bigmemory","raster","adehabitatLT","data.table")) %dopar% {
    
    baselineClimaticVar.bm <- attach.big.matrix(baselineClimaticVar.bm.desc)
    focalBaselineClimate.bm <- attach.big.matrix(focalBaselineClimate.bm.desc)
    futureAnalogsClimate.bm <- attach.big.matrix(futureAnalogsClimate.bm.desc)
    
    # Principal component truncation rule
    trunc.SDs <- 0.1
    
    # Principal components analysis at the focal cell
    PCA <- prcomp( baselineClimaticVar.bm[ baselineClimaticVar.CellLoc[ cell == Cells[focal.i], id] , ] )
    
    if( PCAType == "standardDeviation" ) { PCs <- max(which(unlist(summary(PCA)[1]) > trunc.SDs)) }
    if( PCAType == "cumulativeImportance" ) { PCs <- min(which(summary(PCA)$importance[3,] >= 1 - trunc.SDs)) }
    
    # -------
    
    # project the pools onto the PCs
    baselineClimaticVar.PC <- predict(PCA, baselineClimaticVar.bm[ baselineClimaticVar.CellLoc[ cell == Cells[focal.i], id] , ] )
    focalBaselineClimate.PC <- predict(PCA,focalBaselineClimate.bm[])
    futureAnalogsClimate.PC <- predict(PCA,matrix(futureAnalogsClimate.bm[focal.i,], nrow=1))
    
    # express PC scores as standardized anomalies of reference interannual variability 
    baselineClimaticVar.PC.sd <- apply(baselineClimaticVar.PC,2,sd, na.rm=T)
    focalBaselineClimate.PC <- sweep(focalBaselineClimate.PC,MARGIN=2,baselineClimaticVar.PC.sd,'/')
    futureAnalogsClimate.PC <- sweep(futureAnalogsClimate.PC,MARGIN=2,baselineClimaticVar.PC.sd,'/')
    
    # plot(focalBaselineClimate.PC[,1:2], col="red")
    # points(futureAnalogsClimate.PC[1],futureAnalogsClimate.PC[2], col="green")
    
    # Euclidean nearest neighbour distance in the z-standardized PCs of interannual climatic variability, i.e. the Mahalanobian nearest neighbour. 
    myDist <- get.knnx(data = focalBaselineClimate.PC[,1:PCs] , query = matrix( futureAnalogsClimate.PC[1:PCs], nrow=1 ),k=1,algorithm="brute")
    myAnalog <- myDist$nn.index
    myDist <- myDist$nn.dist
    
    # percentile of the nearest neighbour distance on the chi distribution with degrees of freedom equaling the dimensionality of the distance measurement (PCs)
    NN.chi <- pchi( as.vector(myDist) , PCs, rel.tol=.Machine$double.eps^0.8) 
    if( NN.chi >= (1-1e-16) ){ NN.chi <- 1-1e-16 }
    
    # values of the chi percentiles on a standard half-normal distribution (chi distribution with one degree of freedom)
    NN.sigma <- qchi(NN.chi,1)
    NN.sigma[is.na(NN.sigma)] <- qchi(1-1e-16,1)
    NN.sigma[NN.sigma >= qchi(1-1e-16,1)] <- qchi(1-1e-16,1)
    
    # --------
    
    return( data.frame( cell = Cells[focal.i], analog = Cells[myAnalog], sigma=NN.sigma) )
    
  }
  
  # ---------------
  
  sigmaDisappearClimate <- foreach(focal.i = 1:nCells , .combine=rbind, .verbose=FALSE, .packages=c("FNN","bigmemory","raster","adehabitatLT")) %dopar% {
    
    # This may need to use future climates as "baselineClimaticVar.bm"
    baselineClimaticVar.bm <- attach.big.matrix(baselineClimaticVar.bm.desc)
    focalBaselineClimate.bm <- attach.big.matrix(focalBaselineClimate.bm.desc)
    futureAnalogsClimate.bm <- attach.big.matrix(futureAnalogsClimate.bm.desc)
    
    # Principal component truncation rule
    trunc.SDs <- 0.1
    
    # Principal components analysis at the focal cell
    PCA <- prcomp( baselineClimaticVar.bm[ baselineClimaticVar.CellLoc[ cell == Cells[focal.i], id] , ] )
    
    if( PCAType == "standardDeviation" ) { PCs <- max(which(unlist(summary(PCA)[1]) > trunc.SDs)) }
    if( PCAType == "cumulativeImportance" ) { PCs <- min(which(summary(PCA)$importance[3,] >= 1 - trunc.SDs)) }
    
    # -------
    
    # project the pools onto the PCs
    baselineClimaticVar.PC <- predict(PCA, baselineClimaticVar.bm[ baselineClimaticVar.CellLoc[ cell == Cells[focal.i], id] , ] )
    focalBaselineClimate.PC <- predict(PCA,matrix(focalBaselineClimate.bm[focal.i,], nrow=1))
    futureAnalogsClimate.PC <- predict(PCA,futureAnalogsClimate.bm[])
    
    # express PC scores as standardized anomalies of reference interannual variability 
    baselineClimaticVar.PC.sd <- apply(baselineClimaticVar.PC,2,sd, na.rm=T)
    focalBaselineClimate.PC <- sweep(focalBaselineClimate.PC,MARGIN=2,baselineClimaticVar.PC.sd,'/')
    futureAnalogsClimate.PC <- sweep(futureAnalogsClimate.PC,MARGIN=2,baselineClimaticVar.PC.sd,'/')
    
    # plot(futureAnalogsClimate.PC[,1:2], col="red")
    # points(focalBaselineClimate.PC[1],focalBaselineClimate.PC[2], col="green")
    
    # Euclidean nearest neighbour distance in the z-standardized PCs of interannual climatic variability, i.e. the Mahalanobian nearest neighbour. 
    myDist <- get.knnx(data = futureAnalogsClimate.PC[,1:PCs] , query = matrix( focalBaselineClimate.PC[1:PCs], nrow=1 ),k=1,algorithm="brute")
    myAnalog <- myDist$nn.index
    myDist <- myDist$nn.dist
    
    # percentile of the nearest neighbour distance on the chi distribution with degrees of freedom equaling the dimensionality of the distance measurement (PCs)
    NN.chi <- pchi( as.vector(myDist) , PCs, rel.tol=.Machine$double.eps^0.8) 
    if( NN.chi >= (1-1e-16) ){ NN.chi <- 1-1e-16 }
    
    # values of the chi percentiles on a standard half-normal distribution (chi distribution with one degree of freedom)
    NN.sigma <- qchi(NN.chi,1)
    NN.sigma[is.na(NN.sigma)] <- qchi(1-1e-16,1)
    NN.sigma[NN.sigma >= qchi(1-1e-16,1)] <- qchi(1-1e-16,1)
    
    # --------
    
    return( data.frame(cell=Cells[focal.i],analog=Cells[myAnalog],sigma=NN.sigma) )
    
  }
  
  stopCluster(Cluster); rm(Cluster); gc(reset=TRUE)
  closeAllConnections()
  Sys.time() - startTime
  
  # ------------------------------------------
  # ------------------------------------------
  
  # Save data
  
  dataStructureResult[match(sigmaNovelClimate$cell,dataStructureResult$cell),"novelClimateSigma"] <- sigmaNovelClimate$sigma
  dataStructureResult[match(sigmaNovelClimate$cell,dataStructureResult$cell),"novelClimateAnalog"] <- sigmaNovelClimate$analog
  dataStructureResult[match(sigmaNovelClimate$cell,dataStructureResult$cell),"disappearClimateSigma"] <- sigmaDisappearClimate$sigma
  dataStructureResult[match(sigmaNovelClimate$cell,dataStructureResult$cell),"disappearClimateAnalog"] <- sigmaDisappearClimate$analog
  
  # Distance to analog
  dataStructureResult[,"xAnalog"] <- dataStructure[match(dataStructureResult$novelClimateAnalog,dataStructure$cell),"x"]
  dataStructureResult[,"yAnalog"] <- dataStructure[match(dataStructureResult$novelClimateAnalog,dataStructure$cell),"y"]
  
  distanceToAnalogF <- function(x) { ifelse(!is.na(dataStructureResult[x,"xAnalog"]) , spDistsN1( as.matrix(dataStructureResult[x,c("x","y")]) , as.matrix(dataStructureResult[x,c("xAnalog","yAnalog")]) , longlat = TRUE  ) , NA) }
  clust <- makeCluster(nCores)
  clusterExport(clust, c("dataStructureResult","spDistsN1"))
  dataStructureResult[,"analogDistance"] <- parSapply(clust, 1:nrow(dataStructureResult) , distanceToAnalogF )
  stopCluster(clust); closeAllConnections()
  
  # Refugia
  refugiaF <- function(x) { dataStructureResult[x,"cell"] == dataStructureResult[x,"novelClimateAnalog"] }
  clust <- makeCluster(nCores)
  clusterExport(clust, c("dataStructureResult"))
  dataStructureResult[,"refugia"] <- parSapply(clust, 1:nrow(dataStructureResult) , refugiaF )
  stopCluster(clust); closeAllConnections()
  
  save(dataStructureResult,file=paste0(resultsFolder,"/novelClimates",paste0(firstToupper(strsplit(ROIRasterNames," ")[[1]]), collapse = ""),".RData"))
  
  # --------
  
  rasterPlot <- rasterFromXYZ(dataStructureResult[,c("x","y","novelClimateSigma")])
  plot(rasterPlot, col = (topo.colors(20)) , main=paste0(scenario," | ", PCAType))
  
  # --------
  
  gc(reset=TRUE)
  
  # ------------------------------------------
  # ------------------------------------------
  
  # Relative contribution of predictors
  
  startTime <- Sys.time()
  Cluster <- makeCluster( nCores )
  registerDoParallel( Cluster )
  
  novelClimatePredictorContrib <- foreach(focal.i = 1:nCells , .combine=rbind, .verbose=FALSE, .packages=c("FNN","bigmemory","raster","adehabitatLT","data.table")) %dopar% {
    
    baselineClimaticVar.bm <- attach.big.matrix(baselineClimaticVar.bm.desc)
    focalBaselineClimate.bm <- attach.big.matrix(focalBaselineClimate.bm.desc)
    futureAnalogsClimate.bm <- attach.big.matrix(futureAnalogsClimate.bm.desc)
    
    # Principal component truncation rule
    trunc.SDs <- 0.15
    
    # Principal components analysis at the focal cell
    PCA <- prcomp( baselineClimaticVar.bm[ baselineClimaticVar.CellLoc[ cell == Cells[focal.i], id] , ] )
    
    if( PCAType == "standardDeviation" ) { PCs <- max(which(unlist(summary(PCA)[1]) > trunc.SDs)) }
    if( PCAType == "cumulativeImportance" ) { PCs <- min(which(summary(PCA)$importance[3,] >= 1 - trunc.SDs)) }
    
    # -------
    
    # project the pools onto the PCs
    baselineClimaticVar.PC <- predict(PCA, baselineClimaticVar.bm[ baselineClimaticVar.CellLoc[ cell == Cells[focal.i], id] , ] )
    focalBaselineClimate.PC <- predict(PCA,focalBaselineClimate.bm[])
    futureAnalogsClimate.PC <- predict(PCA,matrix(futureAnalogsClimate.bm[focal.i,], nrow=1))
    
    # express PC scores as standardized anomalies of reference interannual variability 
    baselineClimaticVar.PC.sd <- apply(baselineClimaticVar.PC,2,sd, na.rm=T)
    focalBaselineClimate.PC <- sweep(focalBaselineClimate.PC,MARGIN=2,baselineClimaticVar.PC.sd,'/')
    futureAnalogsClimate.PC <- sweep(futureAnalogsClimate.PC,MARGIN=2,baselineClimaticVar.PC.sd,'/')
    
    myDist <- get.knnx(data = focalBaselineClimate.PC[,1:PCs] , query = matrix( futureAnalogsClimate.PC[1:PCs], nrow=1 ),k=1,algorithm="brute")$nn.dist
    results.i <- numeric(length(predictors))
    
    for( int in 1:length(predictors)) {
      
      focalBaselineClimate.PC.Bk <- focalBaselineClimate.PC
      focalBaselineClimate.PC[,int] <- futureAnalogsClimate.PC[,int]
      myDist.int <- get.knnx(data = focalBaselineClimate.PC[,1:PCs] , query = matrix( futureAnalogsClimate.PC[1:PCs], nrow=1 ),k=1,algorithm="brute")$nn.dist
      results.i[int] <- myDist - myDist.int
      focalBaselineClimate.PC <- focalBaselineClimate.PC.Bk
      
    }
    
    results.i <- (results.i / sum( results.i ) ) * 100
    results.i <- data.frame( Cells[focal.i],t(results.i))
    
    # --------
    
    return( results.i )
    
  }
  
  colnames(novelClimatePredictorContrib) <- c("cell", predictors)
  
  # ---------------
  
  disappearClimatePredictorContrib <- foreach(focal.i = 1:nCells , .combine=rbind, .verbose=FALSE, .packages=c("FNN","bigmemory","raster","adehabitatLT")) %dopar% {
    
    baselineClimaticVar.bm <- attach.big.matrix(baselineClimaticVar.bm.desc)
    focalBaselineClimate.bm <- attach.big.matrix(focalBaselineClimate.bm.desc)
    futureAnalogsClimate.bm <- attach.big.matrix(futureAnalogsClimate.bm.desc)
    
    # Principal component truncation rule
    trunc.SDs <- 0.1
    
    # Principal components analysis at the focal cell
    PCA <- prcomp( baselineClimaticVar.bm[ baselineClimaticVar.CellLoc[ cell == Cells[focal.i], id] , ] )
    
    if( PCAType == "standardDeviation" ) { PCs <- max(which(unlist(summary(PCA)[1]) > trunc.SDs)) }
    if( PCAType == "cumulativeImportance" ) { PCs <- min(which(summary(PCA)$importance[3,] >= 1 - trunc.SDs)) }
    
    # -------
    
    # project the pools onto the PCs
    baselineClimaticVar.PC <- predict(PCA, baselineClimaticVar.bm[ baselineClimaticVar.CellLoc[ cell == Cells[focal.i], id] , ] )
    focalBaselineClimate.PC <- predict(PCA,matrix(focalBaselineClimate.bm[focal.i,], nrow=1))
    futureAnalogsClimate.PC <- predict(PCA,futureAnalogsClimate.bm[])
    
    # express PC scores as standardized anomalies of reference interannual variability 
    baselineClimaticVar.PC.sd <- apply(baselineClimaticVar.PC,2,sd, na.rm=T)
    focalBaselineClimate.PC <- sweep(focalBaselineClimate.PC,MARGIN=2,baselineClimaticVar.PC.sd,'/')
    futureAnalogsClimate.PC <- sweep(futureAnalogsClimate.PC,MARGIN=2,baselineClimaticVar.PC.sd,'/')
    
    # Euclidean nearest neighbour distance in the z-standardized PCs of interannual climatic variability, i.e. the Mahalanobian nearest neighbour. 
    myDist <- get.knnx(data = futureAnalogsClimate.PC[,1:PCs] , query = matrix( focalBaselineClimate.PC[1:PCs], nrow=1 ),k=1,algorithm="brute")$nn.dist
    results.i <- numeric(length(predictors))
    
    for( int in 1:length(predictors)) {
      futureAnalogsClimate.PC.Bk <- futureAnalogsClimate.PC
      futureAnalogsClimate.PC[,int] <- focalBaselineClimate.PC[,int]
      myDist.int <- get.knnx(data = futureAnalogsClimate.PC[,1:PCs] , query = matrix( focalBaselineClimate.PC[1:PCs], nrow=1 ),k=1,algorithm="brute")$nn.dist
      results.i[int] <- myDist - myDist.int
      futureAnalogsClimate.PC <- futureAnalogsClimate.PC.Bk
    }
    
    results.i <- (results.i / sum( results.i ) ) * 100
    results.i <- data.frame( Cells[focal.i],t(results.i))
    
    # --------
    
    return( results.i )
    
  }
  
  stopCluster(Cluster); rm(Cluster); gc(reset=TRUE)
  closeAllConnections()
  Sys.time() - startTime
  
  colnames(disappearClimatePredictorContrib) <- c("cell", predictors)
  
  # -----------
  
  save(novelClimatePredictorContrib, file=paste0(resultsFolder,"/novelClimatePredictorContrib.RData"))
  save(disappearClimatePredictorContrib, file=paste0(resultsFolder,"/disappearClimatePredictorContrib.RData"))
  
  # -----------
  
}

# -----------------------------------------------------
# -----------------------------------------------------

file.remove(list.files(tempFolder, pattern="\\.bin", full.names=T))
file.remove(list.files(tempFolder, pattern="\\.desc", full.names=T))