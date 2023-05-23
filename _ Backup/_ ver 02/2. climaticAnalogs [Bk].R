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

closeAllConnections()
rm(list=(ls()))
gc(reset=TRUE)
source("mainFunctions.R")
nCores <- 16

# library(credentials)
# set_github_pat()

# ----------------------------------------------------------------------
# ----------------------------------------------------------------------

scenario <- "ssp119" # ssp119 ssp585

dataFolder <- "../../Data/"
tempFolder <- "../../Temp/"
climateDataBaselineDir <- paste0(dataFolder,"/Climate/Baseline/")
climateDataProjectionDir <- paste0(dataFolder,"/Climate/",scenario,"/")

resultsSubFolder <- paste0(scenario,"/")
resultsFolder <- paste0("../../Results/",resultsSubFolder)
# unlink(resultsFolder, recursive=T)

# -------------------

predictors <- c("DissolvedMolecularOxygen Surface Mean",
                "OceanTemperature Surface Max",
                "pH Surface Mean",
                "TotalPhytoplankton Surface Mean")

# -------------------

bathymetryReclass <- NULL
calcRelativeContribution <- FALSE
overwrite <- FALSE

# Study type

# Global 
# ROIRasterType <- "continuous" # continuous binomial
# ROIRasterList <- list.files(climateDataBaselineDir,full.names = TRUE)[1]
# ROIRasterNames <- "global"
# PCAType <- "standardDeviation" # cumulativeImportance standardDeviation

# Region based
rangeMapsDir <- "../../Data/Rasters/"
ROIRasterType <- "binomial" # continuous binomial
ROIRasterList <- list.files(rangeMapsDir, full.names = TRUE)
ROIRasterNames <- list.files(rangeMapsDir, full.names = FALSE)
PCAType <- "standardDeviation" # cumulativeImportance standardDeviation

subsetSpecies <- "../../Data/AquacultureSpeciesListAll.csv"
subsetSpecies <- unique(read.csv(subsetSpecies)[,1])
subsetSpecies <- unlist(sapply(subsetSpecies, function(x) { which(grepl(x,ROIRasterNames)) } ))

ROIRasterList <- ROIRasterList[subsetSpecies]
ROIRasterNames <- ROIRasterNames[subsetSpecies]

subsetAnalysisHemisphere <- TRUE
subsetGlobalBaselineBuffer <- 50 # NULL

# ------------------------------------------------------
# ------------------------------------------------------

if(! dir.exists(paste0(resultsFolder))) { dir.create(paste0(resultsFolder), recursive = T) }
if(! dir.exists(paste0(tempFolder))) { dir.create(paste0(tempFolder), recursive = T) }

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

names(dataStructureRaster) <- "dataStructure"
save(dataStructureRaster, file=paste0(resultsFolder,"/dataStructureRaster.RData"))

# -------------------------------------------------
# -------------------------------------------------
# Climate Data

file.remove(list.files(tempFolder, pattern="\\.desc", full.names=T))
file.remove(list.files(tempFolder, pattern="\\.bin", full.names=T))

# -------------

# Get ICV data [historically experienced at the focal stations] and do standardization [divide by standard deviation]

ICVPeriod <- 2000:2009

climaticBaselineConditionsFiles <- list.files(climateDataBaselineDir, full.names = TRUE)
climaticBaselineConditionsFiles <- climaticBaselineConditionsFiles[unique(c(t(sapply(ICVPeriod, function(x) { which(grepl(x,climaticBaselineConditionsFiles)) } ))))]

baselineClimaticVar.Loc <- data.table(cell=as.numeric(sapply(dataStructure$cell, function(x) { rep(x, length(ICVPeriod) ) } )),
                                      time=as.numeric(sapply(dataStructure$cell, function(x) { 1:length(ICVPeriod) } )) )
setkey(baselineClimaticVar.Loc,cell)

baselineClimaticVar.bm <- big.matrix(nrow=nrow(baselineClimaticVar.Loc), ncol=length(predictors) , backingpath=tempFolder , backingfile = "baselineClimaticVarBM.bin", descriptorfile = "baselineClimaticVarBM.desc")
baselineClimaticVar.bm.desc <- dget( paste0(tempFolder,"/baselineClimaticVarBM.desc") )
baselineClimaticVar.bm <- attach.big.matrix(baselineClimaticVar.bm.desc)

for( pred in 1:length(predictors)) {
  
  predFiles <- climaticBaselineConditionsFiles[grep(predictors[pred],climaticBaselineConditionsFiles)]

  # ICV standard deviation per predictor
  assign(paste0("baselineClimaticVarSDPred",pred) , calc(stack(predFiles),sd) )
  assign(paste0("baselineClimaticVarMeanPred",pred) , calc(stack(predFiles),mean) )
  
  for( t in 1:length(predFiles) ) {
    
    dataRaster <- raster(predFiles[t])
    
    # Divide per ICV standard deviation
    dataRaster <- ( dataRaster -  get( paste0("baselineClimaticVarMeanPred",pred) ) ) / 1 + get( paste0("baselineClimaticVarSDPred",pred) )
    # dataRaster <- dataRaster / ( 1 + get( paste0("baselineClimaticVarSDPred",pred) ) )
    
    dataRaster <- crop(dataRaster,dataStructureRaster)
    dataRaster <- mask(dataRaster,dataStructureRaster)
    baselineClimaticVar.bm[ which(baselineClimaticVar.Loc$time == t) , pred ] <- dataRaster[ unlist(baselineClimaticVar.Loc[ time == t ,"cell"]) ]
    
  }
}

# -------------

BaselinePeriod <- 2010:2018

climaticBaselineConditionsFiles <- list.files(climateDataBaselineDir, full.names = TRUE)
climaticBaselineConditionsFiles <- climaticBaselineConditionsFiles[unique(c(t(sapply(BaselinePeriod, function(x) { which(grepl(x,climaticBaselineConditionsFiles)) } ))))]
baselineClimaticConditions.Loc <- data.table(cell=as.numeric(sapply(dataStructure$cell, function(x) { rep(x, length(BaselinePeriod) ) } )),
                                             time=as.numeric(sapply(dataStructure$cell, function(x) { 1:length(BaselinePeriod) } )) )
setkey(baselineClimaticConditions.Loc,cell)

baselineClimaticConditions.bm <- big.matrix(nrow=nrow(baselineClimaticConditions.Loc), ncol=length(predictors) , backingpath=tempFolder , backingfile = "baselineClimaticConditionsBM.bin", descriptorfile = "baselineClimaticConditionsBM.desc")
baselineClimaticConditions.bm.desc <- dget( paste0(tempFolder,"/baselineClimaticConditionsBM.desc") )
baselineClimaticConditions.bm <- attach.big.matrix(baselineClimaticConditions.bm.desc)

for( pred in 1:length(predictors)) {
  
  predFiles <- climaticBaselineConditionsFiles[grep(predictors[pred],climaticBaselineConditionsFiles)]

  for( t in 1:length(predFiles) ) {
    
    dataRaster <- raster(predFiles[t])
    
    # Divide per ICV standard deviation
    dataRaster <- ( dataRaster -  get( paste0("baselineClimaticVarMeanPred",pred) ) ) / 1 + get( paste0("baselineClimaticVarSDPred",pred) )
    # dataRaster <- dataRaster / ( 1 + get( paste0("baselineClimaticVarSDPred",pred) ) )

    dataRaster <- crop(dataRaster,dataStructureRaster)
    dataRaster <- mask(dataRaster,dataStructureRaster)
    baselineClimaticConditions.bm[ which(baselineClimaticConditions.Loc$time == t) , pred ] <- dataRaster[ unlist(baselineClimaticConditions.Loc[ time == t ,"cell"]) ]
    
  }
}

rasterPlot <- dataStructureRaster
rasterPlot[] <- NA
rasterPlot[ unlist(baselineClimaticConditions.Loc[ time == 1 ,"cell"]) ] <- baselineClimaticConditions.bm[ which(baselineClimaticConditions.Loc$time == 1),2]
plot(rasterPlot, col = topo.colors(20))
plot(baselineClimaticConditions.bm[ which(baselineClimaticConditions.Loc$cell == 1000),2])

# -------------

FuturePeriod <- 2091:2100

climaticProjectionConditionsFiles <- list.files(climateDataProjectionDir, full.names = TRUE)
climaticProjectionConditionsFiles <- climaticProjectionConditionsFiles[unique(c(t(sapply(FuturePeriod, function(x) { which(grepl(x,climaticProjectionConditionsFiles)) } ))))]
futureClimaticConditions.Loc <- data.table(cell=as.numeric(sapply(dataStructure$cell, function(x) { rep(x, length(FuturePeriod) ) } )),
                                           time=as.numeric(sapply(dataStructure$cell, function(x) { 1:length(FuturePeriod) } )) )
setkey(futureClimaticConditions.Loc,cell)

futureClimaticConditions.bm <- big.matrix(nrow=nrow(futureClimaticConditions.Loc), ncol=length(predictors) , backingpath=tempFolder , backingfile = "futureClimaticConditionsBM.bin", descriptorfile = "futureClimaticConditionsBM.desc")
futureClimaticConditions.bm.desc <- dget( paste0(tempFolder,"/futureClimaticConditionsBM.desc") )
futureClimaticConditions.bm <- attach.big.matrix(futureClimaticConditions.bm.desc)

for( pred in 1:length(predictors)) {
  
  predFiles <- climaticProjectionConditionsFiles[grep(predictors[pred],climaticProjectionConditionsFiles)]
  
  for( t in 1:length(predFiles) ) {
    
    dataRaster <- raster(predFiles[t])
    
    # Divide per ICV standard deviation
    dataRaster <- ( dataRaster -  get( paste0("baselineClimaticVarMeanPred",pred) ) ) / 1 + get( paste0("baselineClimaticVarSDPred",pred) )
    # dataRaster <- dataRaster / ( 1 + get( paste0("baselineClimaticVarSDPred",pred) ) )
    
    dataRaster <- crop(dataRaster,dataStructureRaster)
    dataRaster <- mask(dataRaster,dataStructureRaster)
    futureClimaticConditions.bm[ which(futureClimaticConditions.Loc$time == t) , pred ] <- dataRaster[ unlist(futureClimaticConditions.Loc[ time == t ,"cell"]) ]
    
  }
}

rasterPlot <- dataStructureRaster
rasterPlot[] <- NA
rasterPlot[ unlist(futureClimaticConditions.Loc[ time == 1 ,"cell"]) ] <- futureClimaticConditions.bm[ which(futureClimaticConditions.Loc$time == 1),2]
plot(rasterPlot, col = topo.colors(20))
plot(futureClimaticConditions.bm[ which(futureClimaticConditions.Loc$cell == 10000),1])

# -------------

rm(baselineClimaticVar.bm)
rm(baselineClimaticConditions.bm)
rm(futureClimaticConditions.bm)

# -------------------------------------------------
# -------------------------------------------------

globalCells <- unique(dataStructure$cell)
length(globalCells)

# ------------------------------------------------------
# ------------------------------------------------------

for( file.i in 1:length(ROIRasterList) ) { # 
  
  cat("\014")
  cat("\n")
  cat("# ---------------------------------\n")
  cat(file.i,"out of",length(ROIRasterList),"\n")
  cat("# ---------------------------------\n")
  cat("\n")
  
  ROIRasterList.i <- ROIRasterList[file.i]
  ROIRasterNames.i <- ROIRasterNames[file.i]

  # Read ROI
  
  ROIRaster <- list.files(ROIRasterList.i, pattern="\\.tif", full.names = TRUE)
  if(class(ROIRaster) == "character" & length(ROIRaster) == 1) { ROIRaster <- raster(ROIRaster) }
  if(class(ROIRaster) == "character" & length(ROIRaster) > 1) { ROIRaster <- calc(stack( ROIRaster ),max,na.rm=T) }

  if( ROIRasterType == "continuous") { ROIRaster <- xyFromCell(ROIRaster,Which(!is.na(ROIRaster), cells=TRUE)) }
  if( ROIRasterType == "binomial") { ROIRaster <- xyFromCell(ROIRaster,Which(ROIRaster == 1, cells=TRUE)) }

  resultsFolder.i <- paste0(resultsFolder,"/",ROIRasterNames.i)
  if(! dir.exists(paste0(resultsFolder.i))) { dir.create(paste0(resultsFolder.i), recursive = T) }
  if( file.exists(paste0(resultsFolder.i,"/climateDissimilarity.RData")) & ! overwrite ) { next }

  focalCells <- cellFromXY(dataStructureRaster,ROIRaster)
  focalCells <- focalCells[ which(focalCells %in% globalCells) ]

  # -----------------------
  
  if( nrow(baselineClimaticVar.Loc) != nrow(attach.big.matrix(baselineClimaticVar.bm.desc)) ) { stop("Error :: 0991")}
  if( nrow(baselineClimaticConditions.Loc) != nrow(attach.big.matrix(baselineClimaticConditions.bm.desc)) ) { stop("Error :: 0992")}
  if( nrow(futureClimaticConditions.Loc) != nrow(attach.big.matrix(futureClimaticConditions.bm.desc)) ) { stop("Error :: 0993")}

  file.remove(list.files(tempFolder, pattern="FocalBM", full.names=T))
  
  baselineClimaticConditionsFocal.Loc <- baselineClimaticConditions.Loc[which(baselineClimaticConditions.Loc$cell %in% focalCells),]
  futureClimaticConditionsFocal.Loc <- futureClimaticConditions.Loc[which(futureClimaticConditions.Loc$cell %in% focalCells),]

  baselineClimaticConditionsFocal <- as.big.matrix( attach.big.matrix(baselineClimaticConditions.bm.desc)[which(baselineClimaticConditions.Loc$cell %in% focalCells),] , backingpath=tempFolder , backingfile = "baselineClimaticConditionsFocalBM.bin", descriptorfile = "baselineClimaticConditionsFocalBM.desc")
  baselineClimaticConditionsFocal.desc <- dget( paste0(tempFolder,"/baselineClimaticConditionsFocalBM.desc") )
  
  futureClimaticConditionsFocal <- as.big.matrix( attach.big.matrix(futureClimaticConditions.bm.desc)[which(futureClimaticConditions.Loc$cell %in% focalCells),] , backingpath=tempFolder , backingfile = "futureClimaticConditionsFocalBM.bin", descriptorfile = "futureClimaticConditionsFocalBM.desc")
  futureClimaticConditionsFocal.desc <- dget( paste0(tempFolder,"/futureClimaticConditionsFocalBM.desc") )

  rm(baselineClimaticConditionsFocal); rm(futureClimaticConditionsFocal)
  
  # ------------------------------------------
  # ------------------------------------------
  
  globalCells.i <- globalCells
  
  if( subsetAnalysisHemisphere ) {
    if( min(ROIRaster[,2]) > 0 ) { 
      globalCells.i <- dataStructure[ dataStructure$y > 0, "cell"]
    }
    if( max(ROIRaster[,2]) < 0 ) { 
      globalCells.i <- dataStructure[ dataStructure$y < 0, "cell"]
    }
  }
  
  if( ! is.null(subsetGlobalBaselineBuffer) ) {
    
    dataStructure.i <- dataStructure[ dataStructure$cell %in% focalCells, ]
    dataStructure.i <- dataStructure[dataStructure$x >= min(dataStructure.i$x) - subsetGlobalBaselineBuffer &
                                     dataStructure$x <= max(dataStructure.i$x) + subsetGlobalBaselineBuffer &
                                     dataStructure$y >= min(dataStructure.i$y) - subsetGlobalBaselineBuffer &
                                     dataStructure$y <= max(dataStructure.i$y) + subsetGlobalBaselineBuffer ,  ]
    
    dataStructure.i <- dataStructure.i[dataStructure.i$cell %in% globalCells.i,]
    dataStructureRaster.i <- crop(dataStructureRaster, extent(min(dataStructure.i$x),max(dataStructure.i$x),min(dataStructure.i$y),max(dataStructure.i$y)))
    dataStructureRaster.i <- clump(dataStructureRaster.i)
    regionsUsed <- unique(extract(dataStructureRaster.i,dataStructure[ dataStructure$cell %in% focalCells,  c("x","y")  ]))
    
    globalCells.i <- dataStructure.i[which(extract(dataStructureRaster.i,dataStructure.i[ ,  c("x","y")  ]) %in% regionsUsed),"cell"]
  
  }

  # plot(dataStructure[ dataStructure$cell %in% globalCells.i,c("x","y")])
  
  # ------------------------------------------
  # ------------------------------------------
  
  time.i <- Sys.time()
  
  Cluster <- makeCluster( nCores )
  registerDoParallel( Cluster )
  
  dataStructureResult <- foreach(cell.i = focalCells , .combine=rbind, .verbose=FALSE, .packages=c("FNN","bigmemory","raster","adehabitatLT","data.table")) %dopar% {
    
    baselineClimaticVar.bm <- attach.big.matrix(baselineClimaticVar.bm.desc)
    baselineClimateGlobalCells.bm <- attach.big.matrix(baselineClimaticConditions.bm.desc)
    baselineClimateFocalCells.bm <- attach.big.matrix(baselineClimaticConditionsFocal.desc)
    futureClimateGlobalCells.bm <- attach.big.matrix(futureClimaticConditions.bm.desc)
    futureClimateFocalCells.bm <- attach.big.matrix(futureClimaticConditionsFocal.desc)
    
    # rm(baselineClimaticVar.bm); rm(baselineClimateGlobalCells.bm); rm(baselineClimateFocalCells.bm); rm(futureClimateGlobalCells.bm); rm(futureClimateFocalCells.bm)
    
    # Principal component truncation rule
    trunc.SDs <- 0.1
    
    # Principal components analysis at the focal cell
    PCA <- prcomp( baselineClimaticVar.bm[ which(baselineClimaticVar.Loc$cell %in% cell.i) , ] )
    
    if( PCAType == "standardDeviation" ) { PCs <- max(which(unlist(summary(PCA)[1]) > trunc.SDs)) }
    if( PCAType == "cumulativeImportance" ) { PCs <- min(which(summary(PCA)$importance[3,] >= 1 - trunc.SDs)) }
    if( length(predictors) > 1 & PCs == 1) { PCs <- 2 }

    # -------

    # Degree of global novelty
    
    # project the pools onto the PCs
    # baselineClimaticVar.PC <- predict(PCA, baselineClimaticVar.bm[ which(baselineClimaticVar.Loc$cell %in% cell.i) , ])
    baselineClimateFocalCells.PC <- predict(PCA,baselineClimateFocalCells.bm[])
    futureClimateGlobalCells.PC <- predict(PCA,futureClimateGlobalCells.bm[ which(futureClimaticConditions.Loc$cell %in% cell.i)  ,])
  
    # express PC scores as standardized anomalies of reference interannual variability 
    # baselineClimaticVar.PC.sd <- 1 + apply(baselineClimaticVar.PC,2,sd, na.rm=T)
    # baselineClimateFocalCells.PC <- sweep(baselineClimateFocalCells.PC,MARGIN=2,baselineClimaticVar.PC.sd,'/')
    # futureClimateGlobalCells.PC <- sweep(futureClimateGlobalCells.PC,MARGIN=2,baselineClimaticVar.PC.sd,'/')
    
    # plot(baselineClimateFocalCells.PC[,1:2], col="red")
    # points(futureClimateGlobalCells.PC[1,1],futureClimateGlobalCells.PC[1,2], col="black")
    
    # Euclidean nearest neighbour distance in the z-standardized PCs of interannual climatic variability, i.e. the Mahalanobian nearest neighbour. 
    myDist <- get.knnx(data = matrix(baselineClimateFocalCells.PC[,1:PCs], ncol=PCs ) , query = matrix( futureClimateGlobalCells.PC[,1:PCs], ncol=PCs ),k=1,algorithm="brute")
    myDist <- myDist$nn.dist
    myDist <- mean(myDist)

    # percentile of the nearest neighbour distance on the chi distribution with degrees of freedom equaling the dimensionality of the distance measurement (PCs)
    NN.chi <- pchi( as.vector(myDist) , PCs, rel.tol=.Machine$double.eps^0.8) 
    if( NN.chi >= (1-1e-16) ){ NN.chi <- 1-1e-16 }
    
    # values of the chi percentiles on a standard half-normal distribution (chi distribution with one degree of freedom)
    NN.sigma <- qchi(NN.chi,1)
    NN.sigma[is.na(NN.sigma)] <- qchi(1-1e-16,1)
    NN.sigma[NN.sigma >= qchi(1-1e-16,1)] <- qchi(1-1e-16,1)
    sigmaNovelty <- NN.sigma
    
    # --------
    # --------
  
    # get analog based on sigma novelty
    
    baselineClimateGlobalCells.PC <- predict(PCA,baselineClimateGlobalCells.bm[which(baselineClimaticConditions.Loc$cell %in% globalCells.i),])
    # baselineClimateGlobalCells.PC <- sweep(baselineClimateGlobalCells.PC,MARGIN=2,baselineClimaticVar.PC.sd,'/')
    myDist <- get.knnx(data = matrix(baselineClimateGlobalCells.PC[,1:PCs], ncol=PCs ) , query = matrix( apply(futureClimateGlobalCells.PC[,1:PCs],2,mean), ncol=PCs ),k=min(10,nrow(baselineClimateGlobalCells.PC)),algorithm="brute")
    
    myDist.dist <- myDist$nn.dist
    myDist.index <- myDist$nn.index
    
    # Get the cell id
    sigmaNoveltyAnalogPool <- baselineClimaticConditions.Loc[ which(baselineClimaticConditions.Loc$cell %in% globalCells.i), cell ]
    sigmaNoveltyAnalogPool <- dataStructure[ dataStructure$cell %in% sigmaNoveltyAnalogPool[myDist.index],c("x","y","cell")]
    
    myDistGeo <- get.knnx(data = as.matrix(sigmaNoveltyAnalogPool[,1:2], ncol=2 ) , query = as.matrix( dataStructure[ dataStructure$cell %in% cell.i,c("x","y")], ncol=2 ),k=1,algorithm="brute")
    sigmaNoveltyAnalog <- sigmaNoveltyAnalogPool[myDistGeo$nn.index,"cell"]
    
    # plot(dataStructure[ dataStructure$cell %in% sigmaNoveltyAnalogPool$cell,c("x","y")])
    # points(dataStructure[ dataStructure$cell %in% cell.i,c("x","y")], col="red")
    # points(dataStructure[ dataStructure$cell %in% sigmaNoveltyAnalog,c("x","y")], col="green")
    
    # --------
    # --------
    
    # Degree of global disappearance
  
    # project the pools onto the PCs
    # baselineClimaticVar.PC <- predict(PCA, baselineClimaticVar.bm[ which(baselineClimaticVar.Loc$cell %in% cell.i) , ] )
    futureClimateFocalCells.PC <- predict(PCA,futureClimateFocalCells.bm[])
    baselineClimateGlobalCells.PC <- predict(PCA,matrix(baselineClimateGlobalCells.bm[ which(baselineClimaticConditions.Loc$cell %in% cell.i)  ,], ncol=length(predictors)))
    
    # express PC scores as standardized anomalies of reference interannual variability 
    # baselineClimaticVar.PC.sd <- apply(baselineClimaticVar.PC,2,sd, na.rm=T)
    # futureClimateFocalCells.PC <- sweep(futureClimateFocalCells.PC,MARGIN=2,baselineClimaticVar.PC.sd,'/')
    # baselineClimateGlobalCells.PC <- sweep(baselineClimateGlobalCells.PC,MARGIN=2,baselineClimaticVar.PC.sd,'/')
    
    # plot(baselineClimateFocalCells.PC[,1:2], col="red")
    # points(futureClimateGlobalCells.PC[1],futureClimateGlobalCells.PC[2], col="black")
    
    # Euclidean nearest neighbour distance in the z-standardized PCs of interannual climatic variability, i.e. the Mahalanobian nearest neighbour. 
    myDist <- get.knnx(data = matrix(futureClimateFocalCells.PC[,1:PCs], ncol=PCs ) , query = matrix( baselineClimateGlobalCells.PC[,1:PCs], ncol=PCs ),k=1,algorithm="brute")
    myDist <- myDist$nn.dist
    myDist <- mean(myDist)
    
    # percentile of the nearest neighbour distance on the chi distribution with degrees of freedom equaling the dimensionality of the distance measurement (PCs)
    NN.chi <- pchi( as.vector(myDist) , PCs, rel.tol=.Machine$double.eps^0.8) 
    if( NN.chi >= (1-1e-16) ){ NN.chi <- 1-1e-16 }
    
    # values of the chi percentiles on a standard half-normal distribution (chi distribution with one degree of freedom)
    NN.sigma <- qchi(NN.chi,1)
    NN.sigma[is.na(NN.sigma)] <- qchi(1-1e-16,1)
    NN.sigma[NN.sigma >= qchi(1-1e-16,1)] <- qchi(1-1e-16,1)
    sigmaDisapper <- NN.sigma
    
    # --------
    # --------
    
    # get analog based on sigma disappear
    
    futureClimateGlobalCells.PC <- predict(PCA,futureClimateGlobalCells.bm[which(futureClimaticConditions.Loc$cell %in% globalCells.i),])
    #futureClimateGlobalCells.PC <- sweep(futureClimateGlobalCells.PC,MARGIN=2,baselineClimaticVar.PC.sd,'/')
    myDist <- get.knnx(data = matrix(futureClimateGlobalCells.PC[,1:PCs], ncol=PCs ) , query = matrix( apply(baselineClimateGlobalCells.PC[,1:PCs],2,mean), ncol=PCs ),k=10,algorithm="brute")
    
    myDist.dist <- myDist$nn.dist
    myDist.index <- myDist$nn.index
    
    # Get the cell id
    sigmaNoveltyAnalogPool <- futureClimaticConditions.Loc[ which(futureClimaticConditions.Loc$cell %in% globalCells.i), cell ]
    sigmaNoveltyAnalogPool <- dataStructure[ dataStructure$cell %in% sigmaNoveltyAnalogPool[myDist.index],c("x","y","cell")]
    
    myDistGeo <- get.knnx(data = as.matrix(sigmaNoveltyAnalogPool[,1:2], ncol=2 ) , query = as.matrix( dataStructure[ dataStructure$cell %in% cell.i,c("x","y")], ncol=2 ),k=1,algorithm="brute")
    sigmaNoveltyAnalog <- sigmaNoveltyAnalogPool[myDistGeo$nn.index,"cell"]
    
    #
    sigmaDisapperAnalog <- futureClimaticConditions.Loc[ which(futureClimaticConditions.Loc$cell %in% globalCells.i), cell ]
    sigmaDisapperAnalog <- sigmaDisapperAnalog[myDist$nn.index[which.min(myDist$nn.dist)]]
    
    # --------
    
    return( data.frame( cell = cell.i, analogNovelty = sigmaNoveltyAnalog, sigmaNovelty=sigmaNovelty, analogDisappearance = sigmaDisapperAnalog, sigmaDisappearance=sigmaDisapper) )
    
  }

  stopCluster(Cluster); rm(Cluster); gc(reset=TRUE)
  closeAllConnections()

  time.i - Sys.time()
    
  # ------------------------------------------
  # ------------------------------------------
  
  # Save data
  
  # dataStructureRaster[] <- NA
  # resultDataStructureRaster <- dataStructureRaster
  # resultDataStructureRaster[dataStructureResult$cell] <- dataStructureResult$sigmaNovelty
  # resultDataStructureRaster <- crop(resultDataStructureRaster,dataStructureRaster.i)
  # plot(resultDataStructureRaster)
  
  dataStructureResult[,"x"] <- dataStructure[match(dataStructureResult$cell,dataStructure$cell),"x"]
  dataStructureResult[,"y"] <- dataStructure[match(dataStructureResult$cell,dataStructure$cell),"y"]
  
  dataStructureResult[,"analogNovelty.x"] <- dataStructure[match(dataStructureResult$analogNovelty,dataStructure$cell),"x"]
  dataStructureResult[,"analogNovelty.y"] <- dataStructure[match(dataStructureResult$analogNovelty,dataStructure$cell),"y"]

  distanceToAnalogF <- function(x) { ifelse(!is.na(dataStructureResult[x,"analogNovelty.x"]) , spDistsN1( as.matrix(dataStructureResult[x,c("x","y")]) , as.matrix(dataStructureResult[x,c("analogNovelty.x","analogNovelty.y")]) , longlat = TRUE  ) , NA) }
  clust <- makeCluster(nCores)
  clusterExport(clust, c("dataStructureResult","spDistsN1"))
  dataStructureResult[,"analogNoveltyDist"] <- parSapply(clust, 1:nrow(dataStructureResult) , distanceToAnalogF )
  stopCluster(clust); closeAllConnections()
  
  dataStructureResult[,"analogDisappearance.x"] <- dataStructure[match(dataStructureResult$analogDisappearance,dataStructure$cell),"x"]
  dataStructureResult[,"analogDisappearance.y"] <- dataStructure[match(dataStructureResult$analogDisappearance,dataStructure$cell),"y"]
  
  distanceToAnalogF <- function(x) { ifelse(!is.na(dataStructureResult[x,"analogDisappearance.x"]) , spDistsN1( as.matrix(dataStructureResult[x,c("x","y")]) , as.matrix(dataStructureResult[x,c("analogDisappearance.x","analogDisappearance.y")]) , longlat = TRUE  ) , NA) }
  clust <- makeCluster(nCores)
  clusterExport(clust, c("dataStructureResult","spDistsN1"))
  dataStructureResult[,"analogDisappearanceDist"] <- parSapply(clust, 1:nrow(dataStructureResult) , distanceToAnalogF )
  stopCluster(clust); closeAllConnections()
  
  # Refugia within focal pool
  refugiaF <- function(x) { dataStructureResult[x,"cell"] %in% dataStructureResult[,"analogNovelty"] }
  clust <- makeCluster(nCores)
  clusterExport(clust, c("dataStructureResult"))
  dataStructureResult[,"refugiaWithinFocalPoll"] <- parSapply(clust, 1:nrow(dataStructureResult) , refugiaF )
  stopCluster(clust); closeAllConnections()
  
  # Refugia 
  refugiaF <- function(x) { dataStructureResult[x,"cell"] == dataStructureResult[x,"analogNovelty"] }
  clust <- makeCluster(nCores)
  clusterExport(clust, c("dataStructureResult"))
  dataStructureResult[,"refugia"] <- parSapply(clust, 1:nrow(dataStructureResult) , refugiaF )
  stopCluster(clust); closeAllConnections()
  
  save(dataStructureResult,file=paste0(resultsFolder.i,"/climateDissimilarity.RData"))
  
  # --------

  gc(reset=TRUE)

  # ------------------------------------------
  # ------------------------------------------
  
  # Relative contribution of predictors
  
  if( calcRelativeContribution ) {
        
      Cluster <- makeCluster( nCores )
      registerDoParallel( Cluster )
      
      novelClimatePredictorContrib <- foreach(focal.i = 1:focalCellsN , .combine=rbind, .verbose=FALSE, .packages=c("FNN","bigmemory","raster","adehabitatLT","data.table")) %dopar% {
        
        baselineClimaticVar.bm <- attach.big.matrix(baselineClimaticVar.bm.desc)
        focalBaselineClimate.bm <- attach.big.matrix(focalBaselineClimate.bm.desc)
        focalFutureClimate.bm <- attach.big.matrix(focalFutureClimate.bm.desc)
        
        # Principal component truncation rule
        trunc.SDs <- 0.1
        
        # Principal components analysis at the focal cell
        PCA <- prcomp( baselineClimaticVar.bm[ baselineClimaticVar.CellLoc[ cell == focalCells[focal.i], id] , ] )
        
        if( PCAType == "standardDeviation" ) { PCs <- max(which(unlist(summary(PCA)[1]) > trunc.SDs)) }
        if( PCAType == "cumulativeImportance" ) { PCs <- min(which(summary(PCA)$importance[3,] >= 1 - trunc.SDs)) }
        
        # -------
        
        # Degree of global novelty
        
        # project the pools onto the PCs
        baselineClimaticVar.PC <- predict(PCA, baselineClimaticVar.bm[ baselineClimaticVar.CellLoc[ cell == focalCells[focal.i], id] , ] )
        focalBaselineClimate.PC <- predict(PCA,focalBaselineClimate.bm[])
        focalFutureClimate.PC <- predict(PCA,matrix(focalFutureClimate.bm[focal.i,], nrow=1))
        
        # express PC scores as standardized anomalies of reference interannual variability 
        baselineClimaticVar.PC.sd <- apply(baselineClimaticVar.PC,2,sd, na.rm=T)
        focalBaselineClimate.PC <- sweep(focalBaselineClimate.PC,MARGIN=2,baselineClimaticVar.PC.sd,'/')
        focalFutureClimate.PC <- sweep(focalFutureClimate.PC,MARGIN=2,baselineClimaticVar.PC.sd,'/')
        
        # plot(focalBaselineClimate.PC[,1:2], col="red")
        # points(focalFutureClimate.PC[1],focalFutureClimate.PC[2], col="green")
        
        # Euclidean nearest neighbour distance in the z-standardized PCs of interannual climatic variability, i.e. the Mahalanobian nearest neighbour. 
        myDist <- get.knnx(data = focalBaselineClimate.PC[,1:PCs] , query = matrix( focalFutureClimate.PC[1:PCs], nrow=1 ),k=1,algorithm="brute")$nn.dist
    
        results.i <- numeric(length(predictors))
        
        for( int in 1:length(predictors)) {
          
          focalBaselineClimate.PC.Bk <- focalBaselineClimate.PC
          focalBaselineClimate.PC[,int] <- focalFutureClimate.PC[,int]
          myDist.int <- get.knnx(data = focalBaselineClimate.PC[,1:PCs] , query = matrix( focalFutureClimate.PC[1:PCs], nrow=1 ),k=1,algorithm="brute")$nn.dist
          results.i[int] <- myDist - myDist.int
          focalBaselineClimate.PC <- focalBaselineClimate.PC.Bk
          
        }
        
        results.i <- (results.i / sum( results.i ) ) * 100
        results.i <- data.frame( cell = Cells[focal.i],t(results.i))
        
        return( results.i )
        
      }
    
      colnames(novelClimatePredictorContrib) <- c("cell", predictors)
      
      # -----------
      
      disappearClimatePredictorContrib <- foreach(focal.i = 1:focalCellsN , .combine=rbind, .verbose=FALSE, .packages=c("FNN","bigmemory","raster","adehabitatLT","data.table")) %dopar% {
        
        baselineClimaticVar.bm <- attach.big.matrix(baselineClimaticVar.bm.desc)
        focalBaselineClimate.bm <- attach.big.matrix(focalBaselineClimate.bm.desc)
        focalFutureClimate.bm <- attach.big.matrix(focalFutureClimate.bm.desc)
        
        # Principal component truncation rule
        trunc.SDs <- 0.1
        
        # Principal components analysis at the focal cell
        PCA <- prcomp( baselineClimaticVar.bm[ baselineClimaticVar.CellLoc[ cell == Cells[focal.i], id] , ] )
        
        if( PCAType == "standardDeviation" ) { PCs <- max(which(unlist(summary(PCA)[1]) > trunc.SDs)) }
        if( PCAType == "cumulativeImportance" ) { PCs <- min(which(summary(PCA)$importance[3,] >= 1 - trunc.SDs)) }
        
        # -------
      
        # Degree of global disappearance
        
        # project the pools onto the PCs
        baselineClimaticVar.PC <- predict(PCA, baselineClimaticVar.bm[ baselineClimaticVar.CellLoc[ cell == focalCells[focal.i], id] , ] )
        focalBaselineClimate.PC <- predict(PCA,matrix(focalBaselineClimate.bm[focal.i,], nrow=1))
        focalFutureClimate.PC <- predict(PCA,focalFutureClimate.bm[])
        
        # express PC scores as standardized anomalies of reference interannual variability 
        baselineClimaticVar.PC.sd <- apply(baselineClimaticVar.PC,2,sd, na.rm=T)
        focalBaselineClimate.PC <- sweep(focalBaselineClimate.PC,MARGIN=2,baselineClimaticVar.PC.sd,'/')
        focalFutureClimate.PC <- sweep(focalFutureClimate.PC,MARGIN=2,baselineClimaticVar.PC.sd,'/')
        
        # plot(focalFutureClimate.PC[,1:2], col="red")
        # points(focalBaselineClimate.PC[1],focalBaselineClimate.PC[2], col="green")
        
        # Euclidean nearest neighbour distance in the z-standardized PCs of interannual climatic variability, i.e. the Mahalanobian nearest neighbour. 
        myDist <- get.knnx(data = focalFutureClimate.PC[,1:PCs] , query = matrix( focalBaselineClimate.PC[1:PCs], nrow=1 ),k=1,algorithm="brute")$nn.dist
    
        results.i <- numeric(length(predictors))
        
        for( int in 1:length(predictors)) {
          
          focalFutureClimate.PC.Bk <- focalFutureClimate.PC
          focalFutureClimate.PC[,int] <- focalBaselineClimate.PC[,int]
          myDist.int <- get.knnx(data = focalFutureClimate.PC[,1:PCs] , query = matrix( focalBaselineClimate.PC[1:PCs], nrow=1 ),k=1,algorithm="brute")$nn.dist
          results.i[int] <- myDist - myDist.int
          focalFutureClimate.PC <- focalFutureClimate.PC.Bk
          
        }
        
        results.i <- (results.i / sum( results.i ) ) * 100
        results.i <- data.frame( cell = Cells[focal.i],t(results.i))
        
        return( results.i )
        
      }
      
      stopCluster(Cluster); rm(Cluster); gc(reset=TRUE)
      closeAllConnections()
      
      colnames(disappearClimatePredictorContrib) <- c("cell", predictors)
      
      # -----------
    
      save(novelClimatePredictorContrib, file=paste0(resultsFolder.i,"/novelClimatePredictorContrib.RData"))
      save(disappearClimatePredictorContrib, file=paste0(resultsFolder.i,"/disappearClimatePredictorContrib.RData"))
      
  }
  
}

# -----------------------------------------------------
# -----------------------------------------------------

file.remove(list.files(tempFolder, pattern="\\.bin", full.names=T))
file.remove(list.files(tempFolder, pattern="\\.desc", full.names=T))
