# ----------------------------------------------------------------------
# ----------------------------------------------------------------------
#
#
# ----------------------------------------------------------------------

# Main references:
 
# https://www.pnas.org/doi/10.1073/pnas.0606292104
# https://pubmed.ncbi.nlm.nih.gov/28145063/
# https://esajournals.onlinelibrary.wiley.com/doi/10.1890/070037
# https://reader.elsevier.com/reader/sd/pii/S2590332221006023?token=0880C0F541F28588CC91307246A7FE2632A408AD7ABB7D548B270CDCFD487C5203D88F6921DF5651A9F5FE043756E658&originRegion=eu-west-1&originCreation=20230512120846
# https://www.sciencedirect.com/science/article/pii/S2590332221006023?via%3Dihub#appsec2
# https://www.nature.com/articles/s41467-019-08540-3
# https://d-nb.info/1155721160/34
# https://iopscience.iop.org/article/10.1088/1748-9326/acc2d4
# https://onlinelibrary.wiley.com/doi/10.1111/gcb.13645
# https://onlinelibrary.wiley.com/doi/full/10.1111/gcb.13645
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8390509/
# https://www.nature.com/articles/s41598-021-94872-4
# https://www.cell.com/one-earth/pdf/S2590-3322%2821%2900602-3.pdf

# ------------------

closeAllConnections()
rm(list=(ls()))
gc(reset=TRUE)
source("mainFunctions.R")
nCores <- 8

# library(credentials)
# set_github_pat() 

# ----------------------------------------------------------------------
# ----------------------------------------------------------------------

scenario <- "ssp585" # ssp119 ssp370 ssp585

dataFolder <- "/Volumes/Dropbox/Dropbox/Manuscripts/Aquaculture exposure to projected climate change/Data/"
tempFolder <- "/Volumes/Dropbox/Dropbox/Manuscripts/Aquaculture exposure to projected climate change/Temp/"
resultsFolder <- "/Volumes/Dropbox/Dropbox/Manuscripts/Aquaculture exposure to projected climate change/Results/G6/"

BaselinePeriod <- 2010:2018
FuturePeriod <- 2091:2100

# G1 
# predictors <- c("Nitrate Surface Mean","OceanTemperature Surface Max","Salinity Surface Mean")

# G2 
# predictors <- c("Phosphate Surface Min","TotalPhytoplankton Surface Mean","OceanTemperature Surface Max","Salinity Surface Mean")

# G3
# predictors <- c("TotalPhytoplankton Surface Mean","OceanTemperature Surface Max","Salinity Surface Mean")

# G4
# predictors <- c("DissolvedMolecularOxygen Surface Mean","OceanTemperature Surface Max","Salinity Surface Mean")

# G5
# predictors <- c("DissolvedMolecularOxygen Surface Mean","TotalPhytoplankton Surface Mean","OceanTemperature Surface Max","Salinity Surface Mean")

# G6
predictors <- c("TotalPhytoplankton Surface Mean","OceanTemperature Surface Max")

bathymetryReclass <- NULL
calcRelativeContribution <- TRUE
calculateMESS <- TRUE
overwrite <- FALSE

# -------------------
# Study type

# Global 
# ROIRasterType <- "continuous" # continuous binomial
# ROIRasterList <- list.files(climateDataBaselineDir,full.names = TRUE)[1]
# ROIRasterNames <- "global"
# PCAType <- "standardDeviation" # cumulativeImportance standardDeviation

# Region based
rangeMapsDir <- "/Volumes/Dropbox/Dropbox/Data/Biodiversity Data/Range Maps [0.25]/"
ROIRasterType <- "binomial" # continuous binomial
ROIRasterList <- list.files(rangeMapsDir, full.names = TRUE)
ROIRasterNames <- list.files(rangeMapsDir, full.names = FALSE)
PCAType <- "standardDeviation" # cumulativeImportance standardDeviation

subsetSpecies <- "/Volumes/Dropbox/Dropbox/Manuscripts/Aquaculture exposure to projected climate change/Data/SpList_G6.csv"
subsetSpecies <- unique(read.csv(subsetSpecies, header=FALSE)[,1])
subsetSpecies <- trimws(subsetSpecies)
subsetSpecies <- unlist(sapply(subsetSpecies, function(x) { which(grepl(x,ROIRasterNames)) } ))

# Get spaces

ROIRasterList <- ROIRasterList[subsetSpecies]
ROIRasterNames <- ROIRasterNames[subsetSpecies]

subsetAnalysisHemisphere <- TRUE
subsetGlobalBaselineBuffer <- 50 # NULL

# ------------------------------------------------------
# ------------------------------------------------------

climateDataBaselineDir <- paste0(dataFolder,"/Climate/Baseline/")
climateDataProjectionDir <- paste0(dataFolder,"/Climate/",scenario,"/")
resultsFolder <- paste0(resultsFolder,scenario,"/")

if(! dir.exists(paste0(resultsFolder))) { dir.create(paste0(resultsFolder), recursive = T) }

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

save(dataStructure, file=paste0(resultsFolder,"/dataStructure.RData"))
names(dataStructureRaster) <- "dataStructure"
save(dataStructureRaster, file=paste0(resultsFolder,"/dataStructureRaster.RData"))

globalCells <- unique(dataStructure$cell)
length(globalCells)

# -------------------------------------------------
# -------------------------------------------------
# Climate Data

climaticBaselineConditionsFiles <- list.files(climateDataBaselineDir, full.names = TRUE)
climaticBaselineConditionsFiles <- climaticBaselineConditionsFiles[unique(c(t(sapply(BaselinePeriod, function(x) { which(grepl(x,climaticBaselineConditionsFiles)) } ))))]

climaticBaselineConditions <- matrix(NA, ncol=length(predictors), nrow=nrow(dataStructure))
colnames(climaticBaselineConditions) <- predictors
climaticBaselineConditions <- data.table(dataStructure[,c("cell","x","y")],climaticBaselineConditions)
setkey(climaticBaselineConditions,cell)

rasterStackBaseline <- list()

for( pred in 1:length(predictors)) {
  
  predFiles <- climaticBaselineConditionsFiles[grep(predictors[pred],climaticBaselineConditionsFiles)]
  dataRaster <- calc(stack(predFiles),mean)
  dataRaster <- crop(dataRaster,dataStructureRaster)
  dataRaster <- mask(dataRaster,dataStructureRaster)
  rasterStackBaseline <- c(rasterStackBaseline,list(dataRaster))

  dataRasterSD <- calc(stack(predFiles),sd)
  dataRasterSD <- crop(dataRasterSD,dataStructureRaster)
  dataRasterSD <- mask(dataRasterSD,dataStructureRaster)
  
  assign(paste0("dataRaster",pred) , dataRaster )
  assign(paste0("dataRasterSD",pred) , dataRasterSD )
  
  dataRaster <- ( get(paste0("dataRaster",pred)) - cellStats( get(paste0("dataRaster",pred)) ,mean,na.rm=T) ) / 1 + get(paste0("dataRasterSD",pred))
  climaticBaselineConditions[,predictors[pred]] <- dataRaster[climaticBaselineConditions$cell]
  
}

rasterStackBaseline <- stack(rasterStackBaseline)

rasterPlot <- dataStructureRaster
rasterPlot[] <- NA
rasterPlot[ climaticBaselineConditions$cell ] <- unlist(climaticBaselineConditions[ , 5])
plot(rasterPlot, col = topo.colors(20))

# -------------

climaticProjectionConditionsFiles <- list.files(climateDataProjectionDir, full.names = TRUE)
climaticProjectionConditionsFiles <- climaticProjectionConditionsFiles[unique(c(t(sapply(FuturePeriod, function(x) { which(grepl(x,climaticProjectionConditionsFiles)) } ))))]

climaticProjectionConditions <- matrix(NA, ncol=length(predictors), nrow=nrow(dataStructure))
colnames(climaticProjectionConditions) <- predictors
climaticProjectionConditions <- data.table(dataStructure[,c("cell","x","y")],climaticProjectionConditions)
setkey(climaticProjectionConditions,cell)

rasterStackProjection <- list()

for( pred in 1:length(predictors)) {
  
  predFiles <- climaticProjectionConditionsFiles[grep(predictors[pred],climaticProjectionConditionsFiles)]
  dataRaster <- calc(stack(predFiles),mean)
  dataRaster <- crop(dataRaster,dataStructureRaster)
  dataRaster <- mask(dataRaster,dataStructureRaster)
  rasterStackProjection <- c(rasterStackProjection,list(dataRaster))
  
  dataRaster <- ( dataRaster - cellStats( get(paste0("dataRaster",pred)) ,mean,na.rm=T) ) / 1 + get(paste0("dataRasterSD",pred))
  climaticProjectionConditions[,predictors[pred]] <- dataRaster[climaticProjectionConditions$cell]
  
}

rasterStackProjection <- stack(rasterStackProjection)

rasterPlot <- dataStructureRaster
rasterPlot[] <- NA
rasterPlot[ climaticProjectionConditions$cell ] <- unlist(climaticProjectionConditions[ , 5])
plot(rasterPlot, col = topo.colors(20))

# ----------

rasterPlot <- dataStructureRaster
rasterPlot[] <- NA
rasterPlot[ climaticProjectionConditions$cell ] <- unlist(climaticProjectionConditions[ , 5]) - unlist(climaticBaselineConditions[ , 5])
plot(rasterPlot, col = topo.colors(20))

# ------------------------------------------------------
# ------------------------------------------------------

for( file.i in 1:length(ROIRasterList) ) {
  
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

  # ---------------------------
  
  ROIRasterMask <- ROIRaster
  ROIRasterMask[ROIRasterMask == 0] <- NA
  rasterStackBaseline.Focal <- mask(rasterStackBaseline,ROIRasterMask)
  rasterStackProjection.Focal <- mask(rasterStackProjection,ROIRasterMask)
  
  # ---------------------------

  if( ROIRasterType == "continuous") { ROIRaster <- xyFromCell(ROIRaster,Which(!is.na(ROIRaster), cells=TRUE)) }
  if( ROIRasterType == "binomial") { ROIRaster <- xyFromCell(ROIRaster,Which(ROIRaster == 1, cells=TRUE)) }

  resultsFolder.i <- paste0(resultsFolder,"/",ROIRasterNames.i)
  if(! dir.exists(paste0(resultsFolder.i))) { dir.create(paste0(resultsFolder.i), recursive = T) }
  if( file.exists(paste0(resultsFolder.i,"/climateDissimilarity.RData")) & ! overwrite ) { next }

  focalCells <- cellFromXY(dataStructureRaster,ROIRaster)
  focalCells <- focalCells[ which(focalCells %in% globalCells) ]
 
  if( length(focalCells) == 1 ) { next }
  
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

  climaticBaselineConditions.Focal <- climaticBaselineConditions[cell %in% focalCells,]
  climaticProjectionConditions.Focal <- climaticProjectionConditions[cell %in% focalCells,]

  climaticBaselineConditions.All <- climaticBaselineConditions[cell %in% globalCells.i,]
  climaticProjectionConditions.All <- climaticProjectionConditions[cell %in% globalCells.i,]

  # Principal component truncation rule
  trunc.SDs <- 0.01
  
  # Principal components analysis at the focal cell
  PCA <- prcomp( climaticBaselineConditions.Focal[,-c(1:3)] )
  
  if( PCAType == "standardDeviation" ) { PCs <- max(which(unlist(summary(PCA)[1]) > trunc.SDs)) }
  if( PCAType == "cumulativeImportance" ) { PCs <- min(which(summary(PCA)$importance[3,] >= 1 - trunc.SDs)) }
  if( length(predictors) > 1 & PCs == 1 | PCs == -Inf ) { PCs <- 2 }
  
  # project the pools onto the PCs
  climaticBaselineConditions.Focal.PC <- predict(PCA, climaticBaselineConditions.Focal[,-c(1:3)] )
  climaticProjectionConditions.Focal.PC <- predict(PCA,climaticProjectionConditions.Focal[,-c(1:3)] )
  climaticBaselineConditions.All.PC <- predict(PCA, climaticBaselineConditions.All[,-c(1:3)] )
  climaticProjectionConditions.All.PC <- predict(PCA,climaticProjectionConditions.All[,-c(1:3)] )

  # ------------------------------------------
  # ------------------------------------------
  
  Cluster <- makeCluster( nCores )
  registerDoParallel( Cluster )
  
  dataStructureResult <- foreach(cell.i = 1:length(focalCells) , .combine=rbind, .verbose=FALSE, .packages=c("FNN","bigmemory","raster","adehabitatLT","data.table")) %dopar% {
    
    # -------

    # Degree of global novelty

    # Euclidean nearest neighbour distance in the z-standardized PCs of interannual climatic variability, i.e. the Mahalanobian nearest neighbour. 
    myDist <- get.knnx(data = climaticBaselineConditions.Focal.PC[,1:PCs] , query = matrix( climaticProjectionConditions.Focal.PC[cell.i,1:PCs], ncol=PCs ),k=1,algorithm="brute")
    myDist <- myDist$nn.dist

    # percentile of the nearest neighbour distance on the chi distribution with degrees of freedom equaling the dimensionality of the distance measurement (PCs)
    NN.chi <- pchi( as.vector(myDist) , PCs, rel.tol=.Machine$double.eps^0.8) 
    if( NN.chi >= (1-1e-16) ){ NN.chi <- 1-1e-16 }
    
    # values of the chi percentiles on a standard half-normal distribution (chi distribution with one degree of freedom)
    NN.sigma <- qchi(NN.chi,1)
    NN.sigma[is.na(NN.sigma)] <- qchi(1-1e-16,1)
    NN.sigma[NN.sigma >= qchi(1-1e-16,1)] <- qchi(1-1e-16,1)
    sigmaNovelty <- NN.sigma
    
    # --------

    # get analog based on sigma novelty
    
    myDist <- get.knnx(data = climaticBaselineConditions.All.PC[,1:PCs] , query =  matrix( climaticProjectionConditions.Focal.PC[cell.i,1:PCs], ncol=PCs ),k=100,algorithm="brute")
    myDist.index <- myDist$nn.index
    myDist.dist <- myDist$nn.dist
    
    myDist.index <- myDist.index[which( myDist.dist <= quantile(myDist.dist,0.05))]
    myDist.dist <- myDist.dist[which( myDist.dist <= quantile(myDist.dist,0.05))]
    myDistGeo <- get.knnx(data = as.matrix(climaticBaselineConditions.All[myDist.index,.(x,y)], ncol=2 ) , query = as.matrix( climaticProjectionConditions.Focal[cell.i,.(x,y)], ncol=2 ),k=1,algorithm="brute")
    
    sigmaNoveltyAnalog <- climaticBaselineConditions.All[myDist.index[myDistGeo$nn.index],cell]
    
    # --------
    # --------
  
    # Degree of global disappearance
    
    # Euclidean nearest neighbour distance in the z-standardized PCs of interannual climatic variability, i.e. the Mahalanobian nearest neighbour. 
    myDist <- get.knnx(data = climaticProjectionConditions.Focal.PC[,1:PCs] , query = matrix( climaticBaselineConditions.Focal.PC[cell.i,1:PCs], ncol=PCs ),k=1,algorithm="brute")
    myDist <- myDist$nn.dist
    
    # percentile of the nearest neighbour distance on the chi distribution with degrees of freedom equaling the dimensionality of the distance measurement (PCs)
    NN.chi <- pchi( as.vector(myDist) , PCs, rel.tol=.Machine$double.eps^0.8) 
    if( NN.chi >= (1-1e-16) ){ NN.chi <- 1-1e-16 }
    
    # values of the chi percentiles on a standard half-normal distribution (chi distribution with one degree of freedom)
    NN.sigma <- qchi(NN.chi,1)
    NN.sigma[is.na(NN.sigma)] <- qchi(1-1e-16,1)
    NN.sigma[NN.sigma >= qchi(1-1e-16,1)] <- qchi(1-1e-16,1)
    sigmaDisapper <- NN.sigma
    
    # --------
    
    # get analog based on sigma disappearance
    
    myDist <- get.knnx(data = climaticProjectionConditions.All.PC[,1:PCs] , query =  matrix( climaticBaselineConditions.Focal.PC[cell.i,1:PCs], ncol=PCs ),k=100,algorithm="brute")
    myDist.index <- myDist$nn.index
    myDist.dist <- myDist$nn.dist
    
    myDist.index <- myDist.index[which( myDist.dist <= quantile(myDist.dist,0.05))]
    myDist.dist <- myDist.dist[which( myDist.dist <= quantile(myDist.dist,0.05))]
    myDistGeo <- get.knnx(data = as.matrix(climaticProjectionConditions.All[myDist.index,.(x,y)], ncol=2 ) , query = as.matrix( climaticBaselineConditions.Focal[cell.i,.(x,y)], ncol=2 ),k=1,algorithm="brute")
    
    sigmaDisapperAnalog <- climaticProjectionConditions.All[myDist.index[myDistGeo$nn.index],cell]
    
    # --------
    # --------
    
    return( data.frame( cell = climaticProjectionConditions.Focal[cell.i,cell], analogNovelty = sigmaNoveltyAnalog, sigmaNovelty=sigmaNovelty, analogDisappearance = sigmaDisapperAnalog, sigmaDisappearance=sigmaDisapper) )
    
  }
  
  stopCluster(Cluster); rm(Cluster); gc(reset=TRUE)
  closeAllConnections()
  
  # ------------------------------------------
  # ------------------------------------------
  
  # hist(dataStructureResult$sigmaNovelty)
  # hist(dataStructureResult$sigmaDisappearance)
  
  # dataStructureRaster[] <- NA
  # resultDataStructureRaster <- dataStructureRaster
  # resultDataStructureRaster[dataStructureResult$cell] <- dataStructureResult$sigmaNovelty
  # resultDataStructureRaster <- crop(resultDataStructureRaster,dataStructureRaster.i)
  # plot(resultDataStructureRaster)
  
  # plot(resultDataStructureRaster > 4)
  # plot(resultDataStructureRaster > 2)

  # ------------------------------------------
  # ------------------------------------------
  
  # Save data

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
  
  # --------

  # Calculate MESS
  
  if( calculateMESS ) {
    
    messCalc <- mess(rasterStackProjection.Focal, as.data.frame(rasterStackBaseline.Focal, na.rm=T), full=T)
    messCalcMESS <- messCalc$rmess
    messCalcMESS[messCalcMESS == Inf ] <- NA
    dataStructureResult[,"messDissimilarity"] <- extract(messCalcMESS ,dataStructureResult[,c("x","y")])
    
  }
  
  # --------
  
  save(dataStructureResult,file=paste0(resultsFolder.i,"/climateDissimilarity.RData"))
  gc(reset=TRUE)

  # ------------------------------------------
  # ------------------------------------------
  
  # Relative contribution of predictors
  
  if( calcRelativeContribution ) {
        
    # -----------
    
    dataStructureResult.pred.novelty <- matrix(NA,ncol=length(predictors),nrow=nrow(dataStructureResult))
    dataStructureResult.pred.disappearance <- matrix(NA,ncol=length(predictors),nrow=nrow(dataStructureResult))
    
    for( pred in predictors ) {
      
      climaticBaselineConditions.Focal <- climaticBaselineConditions[cell %in% focalCells,]
      climaticProjectionConditions.Focal <- climaticProjectionConditions[cell %in% focalCells,]
      
      climaticBaselineConditions.All <- climaticBaselineConditions[cell %in% globalCells.i,]
      climaticProjectionConditions.All <- climaticProjectionConditions[cell %in% globalCells.i,]
      
      climaticProjectionConditions.Focal[,which(colnames(climaticProjectionConditions.Focal) == pred) := climaticBaselineConditions.Focal[,get(pred)] ] 
      climaticProjectionConditions.All[,which(colnames(climaticProjectionConditions.All) == pred) := climaticBaselineConditions.All[,get(pred)] ] 

      # Principal component truncation rule
      trunc.SDs <- 0.01
      
      # Principal components analysis at the focal cell
      PCA <- prcomp( climaticBaselineConditions.Focal[,-c(1:3)] )
      
      if( PCAType == "standardDeviation" ) { PCs <- max(which(unlist(summary(PCA)[1]) > trunc.SDs)) }
      if( PCAType == "cumulativeImportance" ) { PCs <- min(which(summary(PCA)$importance[3,] >= 1 - trunc.SDs)) }
      if( length(predictors) > 1 & PCs == 1) { PCs <- 2 }
      
      # project the pools onto the PCs
      climaticBaselineConditions.Focal.PC <- predict(PCA, climaticBaselineConditions.Focal[,-c(1:3)] )
      climaticProjectionConditions.Focal.PC <- predict(PCA,climaticProjectionConditions.Focal[,-c(1:3)] )
      climaticBaselineConditions.All.PC <- predict(PCA, climaticBaselineConditions.All[,-c(1:3)] )
      climaticProjectionConditions.All.PC <- predict(PCA,climaticProjectionConditions.All[,-c(1:3)] )
      
      Cluster <- makeCluster( nCores )
      registerDoParallel( Cluster )
      
      dataStructureResult.pred.i <- foreach(cell.i = 1:length(focalCells) , .combine=rbind, .verbose=FALSE, .packages=c("FNN","bigmemory","raster","adehabitatLT","data.table")) %dopar% {
        
        # -------
        
        # Degree of global novelty
        
        # Euclidean nearest neighbour distance in the z-standardized PCs of interannual climatic variability, i.e. the Mahalanobian nearest neighbour. 
        myDist <- get.knnx(data = climaticBaselineConditions.Focal.PC[,1:PCs] , query = matrix( climaticProjectionConditions.Focal.PC[cell.i,1:PCs], ncol=PCs ),k=1,algorithm="brute")
        myDist <- myDist$nn.dist
        
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
        
        # Degree of global disappearance
        
        # Euclidean nearest neighbour distance in the z-standardized PCs of interannual climatic variability, i.e. the Mahalanobian nearest neighbour. 
        myDist <- get.knnx(data = climaticProjectionConditions.Focal.PC[,1:PCs] , query = matrix( climaticBaselineConditions.Focal.PC[cell.i,1:PCs], ncol=PCs ),k=1,algorithm="brute")
        myDist <- myDist$nn.dist
        
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
        
        return( data.frame( sigmaNovelty=sigmaNovelty, sigmaDisappearance=sigmaDisapper) )
        
      }
      
      stopCluster(Cluster); rm(Cluster); gc(reset=TRUE)
      closeAllConnections()
      
      dataStructureResult.pred.novelty[,which(predictors == pred)] <- dataStructureResult.pred.i$sigmaNovelty - dataStructureResult$sigmaNovelty
      dataStructureResult.pred.disappearance[,which(predictors == pred)] <- dataStructureResult.pred.i$sigmaDisappearance - dataStructureResult$sigmaDisappearance
      
    }
  
    # -----------
  
    save(dataStructureResult.pred.novelty, file=paste0(resultsFolder.i,"/climateDissimilarityPredictorContribNovelty.RData"))
    save(dataStructureResult.pred.disappearance, file=paste0(resultsFolder.i,"/climateDissimilarityPredictorContribDisapper.RData"))
    
  }
  
}

# -----------------------------------------------------
# -----------------------------------------------------