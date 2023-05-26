# ----------------------------------------------------------------------
# ----------------------------------------------------------------------
#
#
# ----------------------------------------------------------------------

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

predictors <- c("DissolvedMolecularOxygen Surface Mean","OceanTemperature Surface Max","Salinity Surface Mean")

# -------------------
# Study type

rangeMapsDir <- "/Volumes/Dropbox/Dropbox/Data/Biodiversity Data/Range Maps [0.25]/"
ROIRasterType <- "binomial" # continuous binomial
ROIRasterList <- list.files(rangeMapsDir, full.names = TRUE)
ROIRasterNames <- list.files(rangeMapsDir, full.names = FALSE)
PCAType <- "standardDeviation" # cumulativeImportance standardDeviation

subsetSpecies <- "Dicentrarchus labrax"
subsetSpecies <- unlist(sapply(subsetSpecies, function(x) { which(grepl(x,ROIRasterNames)) } ))
ROIRasterList <- ROIRasterList[subsetSpecies]
ROIRasterNames <- ROIRasterNames[subsetSpecies]
subsetAnalysisHemisphere <- TRUE
subsetGlobalBaselineBuffer <- 50 # NULL

# ------------------------------------------------------
# ------------------------------------------------------

climateDataBaselineDir <- paste0(dataFolder,"/Climate/Baseline/")
climateDataProjectionDir <- paste0(dataFolder,"/Climate/",scenario,"/")
resultsFolder <- paste0(resultsFolder,scenario,"/")

# -------------------
# -------------------
# Data Structure

dataStructureRaster <- raster(list.files(climateDataBaselineDir,full.names = TRUE)[1])
dataStructureRaster[!is.na(dataStructureRaster)] <- 1

dataStructure <- Which(!is.na(dataStructureRaster), cells=TRUE )
dataStructure <- data.frame(cell=dataStructure,
                            row=rowFromCell(dataStructureRaster,dataStructure),
                            column=colFromCell(dataStructureRaster,dataStructure),
                            xyFromCell(dataStructureRaster,dataStructure))

names(dataStructureRaster) <- "dataStructure"

globalCells <- unique(dataStructure$cell)

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

# ------------------------------------------------------
# ------------------------------------------------------

# Read ROI

ROIRaster <- list.files(ROIRasterList, pattern="\\.tif", full.names = TRUE)
if(class(ROIRaster) == "character" & length(ROIRaster) == 1) { ROIRaster <- raster(ROIRaster) }
if(class(ROIRaster) == "character" & length(ROIRaster) > 1) { ROIRaster <- calc(stack( ROIRaster ),max,na.rm=T) }
plot(ROIRaster)

maskRegion <- shapefile("../../../Manuscripts/Aquaculture exposure to projected climate change/Data/fineTuneRangeMaps/Salmo salar.shp")
ROIRaster <- mask(ROIRaster,maskRegion)
plot(ROIRaster)

# ---------------------------

ROIRasterMask <- ROIRaster
ROIRasterMask[ROIRasterMask == 0] <- NA
rasterStackBaseline.Focal <- mask(rasterStackBaseline,ROIRasterMask)
rasterStackProjection.Focal <- mask(rasterStackProjection,ROIRasterMask)

# ---------------------------

if( ROIRasterType == "continuous") { ROIRaster <- xyFromCell(ROIRaster,Which(!is.na(ROIRaster), cells=TRUE)) }
if( ROIRasterType == "binomial") { ROIRaster <- xyFromCell(ROIRaster,Which(ROIRaster == 1, cells=TRUE)) }

focalCells <- cellFromXY(dataStructureRaster,ROIRaster)
focalCells <- focalCells[ which(focalCells %in% globalCells) ]

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

dataStructureRaster[] <- NA
resultDataStructureRaster <- dataStructureRaster
resultDataStructureRaster[dataStructureResult$cell] <- dataStructureResult$sigmaNovelty
resultDataStructureRaster <- crop(resultDataStructureRaster,dataStructureRaster.i)
plot(resultDataStructureRaster)

dataStructureResult[,"x"] <- dataStructure[match(dataStructureResult$cell,dataStructure$cell),"x"]
dataStructureResult[,"y"] <- dataStructure[match(dataStructureResult$cell,dataStructure$cell),"y"]

dataStructureResult[,"analogNovelty.x"] <- dataStructure[match(dataStructureResult$analogNovelty,dataStructure$cell),"x"]
dataStructureResult[,"analogNovelty.y"] <- dataStructure[match(dataStructureResult$analogNovelty,dataStructure$cell),"y"]

# -------------

which(dataStructureResult$sigmaNovelty < 2 & 
        dataStructureResult$x > 0 & 
        dataStructureResult$x < 10 &
        dataStructureResult$y > 47.5 )

plot(resultDataStructureRaster)
points(dataStructureResult[3093,c("x","y")])
points(dataStructureResult[3093,c("analogNovelty.x","analogNovelty.y")], col="red")

# 3093 - 3321 // 8.292361
# 400 - 2018 // 1.030699

cell.i <- 400
myDist <- get.knnx(data = climaticBaselineConditions.Focal.PC[,1:PCs] , query = matrix( climaticProjectionConditions.Focal.PC[cell.i,1:PCs], ncol=PCs ),k=1,algorithm="brute")
myDist$nn.index

myDist <- get.knnx(data = climaticBaselineConditions.Focal.PC[,1:PCs] , query = matrix( climaticProjectionConditions.Focal.PC[cell.i,1:PCs], ncol=PCs ),k=1,algorithm="brute")
myDist <- myDist$nn.dist
NN.chi <- pchi( as.vector(myDist) , PCs, rel.tol=.Machine$double.eps^0.8) 
if( NN.chi >= (1-1e-16) ){ NN.chi <- 1-1e-16 }
NN.sigma <- qchi(NN.chi,1)
NN.sigma[is.na(NN.sigma)] <- qchi(1-1e-16,1)
NN.sigma[NN.sigma >= qchi(1-1e-16,1)] <- qchi(1-1e-16,1)
NN.sigma

# -------------

dataStructureResult[3093,]

A <- climaticBaselineConditions.Focal.PC[3321,1:2]
A1 <- climaticProjectionConditions.Focal.PC[3093,1:2]

B <- climaticBaselineConditions.Focal.PC[2018,1:2] + c(-3,+2)
B1 <- climaticProjectionConditions.Focal.PC[400,1:2] + c(+4,+2)

ggplot() + 
  geom_point(data=data.frame(climaticBaselineConditions.Focal.PC), aes(x=PC1, y=PC2), color="gray") +
  
  geom_line(data=data.frame(x=c(A[1],A1[1]),y=c(A[2],A1[2]),grp=1), linetype=3, size=0.5, aes(x=x, y=y, group = grp)) +
  geom_point(data=data.frame(x=A[1],y=A[2]), aes(x=x, y=y), color="black", size=2) +
  geom_point(data=data.frame(x=A1[1],y=A1[2]), aes(x=x, y=y), color="red", size=2) +
  geom_point(data=data.frame(x=A[1],y=A[2]), aes(x=x, y=y), color="black", pch=21, size=5) +
  geom_point(data=data.frame(x=A1[1],y=A1[2]), aes(x=x, y=y), color="red", pch=21, size=5) +
  
  geom_line(data=data.frame(x=c(B[1],B1[1]),y=c(B[2],B1[2]),grp=1), linetype=3, size=0.5, aes(x=x, y=y, group = grp)) +
  geom_point(data=data.frame(x=B[1],y=B[2]), aes(x=x, y=y), color="black", size=2) +
  geom_point(data=data.frame(x=B1[1],y=B1[2]), aes(x=x, y=y), color="red", size=2) +
  geom_point(data=data.frame(x=B[1],y=B[2]), aes(x=x, y=y), color="black", pch=21, size=5) +
  geom_point(data=data.frame(x=B1[1],y=B1[2]), aes(x=x, y=y), color="red", pch=21, size=5) +
  
  theme_light() +
  
  annotate(geom="text", x=A[1],A[2], label="A",size=4.25,family="Helvetica", color = "#444444",hjust = -1) +
  annotate(geom="text", x=A1[1],A1[2], label="A'",size=4.25,family="Helvetica", color = "#444444",hjust = -1) +
  annotate(geom="text", x=B[1],B[2], label="B",size=4.25,family="Helvetica", color = "#444444",hjust = -1) +
  annotate(geom="text", x=B1[1],B1[2], label="B'",size=4.25,family="Helvetica", color = "#444444",hjust = -1) 

# ------------------------------------------
# ------------------------------------------

# hist(dataStructureResult$sigmaNovelty)
# hist(dataStructureResult$sigmaDisappearance)

which(dataStructureResult$sigmaNovelty < 2 & 
        dataStructureResult$x < -25 )

plot(resultDataStructureRaster)
points(dataStructureResult[4916,c("x","y")])
points(dataStructureResult[4916,c("analogNovelty.x","analogNovelty.y")], col="red")


# 3093 - 3321 // 8.292361
# 400 - 2018 // 1.030699

themeMap <- 
  theme_minimal() +
  theme(
    text = element_text(family = "Helvetica", color = "#22211d"),
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major = element_line(color = "black", size = 0),
    panel.grid.minor = element_blank(),
    plot.background = element_rect(fill = "#FFFFFF", color = NA), # F3F3F3
    panel.background = element_rect(fill = "#FFFFFF", color = NA), # F3F3F3
    legend.background = element_rect(fill = "#FFFFFF", color = NA),
    legend.box.background = element_rect(fill='#FFFFFF'),
    panel.border = element_blank(),
    legend.position = "none"
  )


"#575757"
"#CDCDCD"

sf_use_s2(FALSE)
mapExtent <- c(-35, 50,20, 75)
names(mapExtent) <- c("xmin","xmax","ymin","ymax")
# worldMapCoordRef <- paste0("+proj=laea +lat_0=",mapExtent[4]-mapExtent[3]," +lon_0=10 +x_0=4321000 +y_0=3210000 +ellps=GRS80 +units=m +no_defs")
worldMapCoordRef <- "+proj=eqdc +lat_1=10 +lat_2=50 +lon_0=20"

worldMap <- ne_countries(scale = 10, returnclass = "sp")
worldMap <- crop(worldMap,mapExtent)
plot(worldMap)

library( rnaturalearth )
library( h3js )
library( h3jsr )
library( sf )
library(maptools)

rasterMap <- resultDataStructureRaster
rasterMap <- crop(rasterMap, extent(mapExtent))
resolutionH3 <- 3
rasterMapDF <- data.frame(xyFromCell(rasterMap, Which( !is.na(rasterMap) , cells=T)),val=rasterMap[Which( !is.na(rasterMap) , cells=T)])

cl <- makeCluster( detectCores() / 2 )
clusterExport(cl, c("resolutionH3","rasterMapDF") )
hexAddress <- parApply(cl, data.frame(rasterMapDF), 1, function(x) { h3js::h3_geo_to_h3(x[[2]], x[[1]], res = resolutionH3) } )
hexAddressHexagon <- parLapply(cl, unique(hexAddress), function(x) { mean(rasterMapDF[rasterMapDF$hex == x , "val"],na.rm=T) } )
stopCluster(cl) 
rasterMapDF <- data.frame(rasterMapDF,hex=hexAddress)
cl <- makeCluster( detectCores() / 2 )
clusterExport(cl, c("resolutionH3","rasterMapDF") )
hexAddressHexagon <- parLapply(cl, unique(hexAddress), function(x) { mean(rasterMapDF[rasterMapDF$hex == x , "val"],na.rm=T) } )
stopCluster(cl)
closeAllConnections()

rasterMapDF <- data.frame(hex=unique(rasterMapDF$hex),val=unlist(hexAddressHexagon))
rasterMapDF.polygons <- h3jsr::cell_to_polygon(input =rasterMapDF$hex, simple = FALSE)
rasterMapDF.polygons$hex <- as.character(rasterMapDF$hex)
rasterMapDF.polygons$value <- as.numeric(as.character(rasterMapDF$val))
minLegend <- min(rasterMapDF.polygons$value)
maxLegend <- max(rasterMapDF.polygons$value)
nColors <- 7
colorBreaks <- seq( min(rasterMapDF.polygons$value),max(rasterMapDF.polygons$value), length.out=nColors)
colorBreaks[1] <- min(rasterMapDF.polygons$value)
colorBreaks[nColors] <- max(rasterMapDF.polygons$value)
myColors <- c("#E4FAFF","#E8EF15","#ec7a06","#e31515", "#450751") # Blue Yellow Red Purple
myColors <- colorRampPalette(myColors)(nColors)
hexagons <- as_Spatial(rasterMapDF.polygons)

worldMapP <-  st_transform(st_as_sf(worldMap),worldMapCoordRef)
rasterMapDF.polygonsP <-  st_transform(rasterMapDF.polygons,worldMapCoordRef)

dataStructureResult[3093,c("x","y")]
dataStructureResult[3093,c("analogNovelty.x","analogNovelty.y")]
dataStructureResult[400,c("x","y")]
dataStructureResult[400,c("analogNovelty.x","analogNovelty.y")]

ggplot() +
  geom_sf(data = worldMapP , fill = "#CDCDCD", colour = "#CDCDCD" , size=0.25 ) +
  geom_sf(data = rasterMapDF.polygonsP, aes(fill=value), colour ="black", size = 0.05) + # round signif
  scale_colour_gradientn(colours = myColors, breaks= colorBreaks, aesthetics = "fill", labels=signif(colorBreaks, digits = 3), limits=c(minLegend,maxLegend) ) +
  themeMap + coord_sf(crs = worldMapCoordRef)

