# ----------------------------------------------------------------------
# ----------------------------------------------------------------------
#
#
# ----------------------------------------------------------------------

closeAllConnections()
rm(list=(ls()[ls()!="v"]))
gc(reset=TRUE)
source("mainFunctions.R")
nCores <- 8

rgdal::setCPLConfigOption("GDAL_PAM_ENABLED", "FALSE")

# ----------------------------------------------------------------------
# ----------------------------------------------------------------------

scenario <- "ssp119" # ssp119 ssp585
producePDF <- FALSE

dataFolder <- "../../Data/"
climateDataBaselineDir <- paste0(dataFolder,"/Climate/Baseline/")
climateDataProjectionDir <- paste0(dataFolder,"/Climate/",scenario,"/")

resultsSubFolder <- paste0(scenario,"/")
resultsFolder <- paste0("../../Results/",resultsSubFolder)

# -----------

dataStructureRaster <- loadRData(paste0(resultsFolder,"/dataStructureRaster.RData"))
dataStructure <- loadRData(paste0(resultsFolder,"/dataStructure.RData"))

resultsFolder <- list.files(resultsFolder, full.names = T)
resultsFolder <- resultsFolder[!grepl("\\.RData",resultsFolder)]

for( resultsFolder.i in resultsFolder) {
  
  if( ! file.exists(paste0(resultsFolder.i,"/climateDissimilarity.RData")) ) { next }
  
  analogResults <- loadRData(paste0(resultsFolder.i,"/climateDissimilarity.RData"))
  if(! dir.exists(paste0(resultsFolder.i,"/Raster/"))) { dir.create(paste0(resultsFolder.i,"/Raster/"), recursive = T) }
  
  # ------------
  
  dataStructureRaster[] <- NA
  resultDataStructureRaster <- dataStructureRaster
  resultDataStructureRaster[analogResults$cell] <- analogResults$sigmaNovelty
  writeRaster(resultDataStructureRaster,file=paste0(resultsFolder.i,"/Raster/sigmaNovelty.tif"), format="GTiff", overwrite=TRUE)
  
  resultDataStructureRaster <- dataStructureRaster
  resultDataStructureRaster[analogResults$cell] <- analogResults$analogNovelty
  writeRaster(resultDataStructureRaster,file=paste0(resultsFolder.i,"/Raster/analogNovelty.tif"), format="GTiff", overwrite=TRUE)
  
  resultDataStructureRaster <- dataStructureRaster
  resultDataStructureRaster[analogResults$cell] <- analogResults$sigmaDisappearance
  writeRaster(resultDataStructureRaster,file=paste0(resultsFolder.i,"/Raster/sigmaDisappearance.tif"), format="GTiff", overwrite=TRUE)
  
  resultDataStructureRaster <- dataStructureRaster
  resultDataStructureRaster[analogResults$cell] <- analogResults$analogDisappearance
  writeRaster(resultDataStructureRaster,file=paste0(resultsFolder.i,"/Raster/analogDisappearance.tif"), format="GTiff", overwrite=TRUE)
  
  resultDataStructureRaster <- dataStructureRaster
  resultDataStructureRaster[analogResults$cell] <- analogResults$refugiaWithinFocalPoll
  writeRaster(resultDataStructureRaster,file=paste0(resultsFolder.i,"/Raster/refugiaWithinFocalPoll.tif"), format="GTiff", overwrite=TRUE)
  
  resultDataStructureRaster <- dataStructureRaster
  resultDataStructureRaster[analogResults$cell] <- analogResults$refugia
  writeRaster(resultDataStructureRaster,file=paste0(resultsFolder.i,"/Raster/refugia.tif"), format="GTiff", overwrite=TRUE)
  
  resultDataStructureRaster <- dataStructureRaster
  resultDataStructureRaster[analogResults$cell] <- analogResults$analogNoveltyDist
  writeRaster(resultDataStructureRaster,file=paste0(resultsFolder.i,"/Raster/analogNoveltyDist.tif"), format="GTiff", overwrite=TRUE)
  
  resultDataStructureRaster <- dataStructureRaster
  resultDataStructureRaster[analogResults$cell] <- analogResults$analogDisappearanceDist
  writeRaster(resultDataStructureRaster,file=paste0(resultsFolder.i,"/Raster/analogDisappearanceDist.tif"), format="GTiff", overwrite=TRUE)
  
  # ------------
  
  if( "climateConnectivity.RData" %in% list.files(resultsFolder.i) ) {
  
    connectivityResults <- loadRData(paste0(resultsFolder.i,"/climateDissimilarity.RData"))
    
    corridor <- dataStructureRaster
    corridor[connectivityResults$cell] <- connectivityResults$corridorQuantile95
    writeRaster(corridor,file=paste0(resultsFolder,"/Raster/corridorQuantile95.tif"), format="GTiff", overwrite=TRUE)
    
    corridor <- dataStructureRaster
    corridor[connectivityResults$cell] <- connectivityResults$corridorQuantile75
    writeRaster(corridor,file=paste0(resultsFolder,"/Raster/corridorQuantile75.tif"), format="GTiff", overwrite=TRUE)
    
    analogDistanceVert <- dataStructureRaster
    analogDistanceVert[connectivityResults$cell] <- connectivityResults$analogDistanceVert
    writeRaster(analogDistanceVert,file=paste0(resultsFolder,"/Raster/analogDistanceVert.tif"), format="GTiff", overwrite=TRUE)
    
  }
   
  # ------------
  
  if( producePDF ) {
    
    if(! dir.exists(paste0(resultsFolder.i,"/Figures/"))) { dir.create(paste0(resultsFolder.i,"/Figures/"), recursive = T) }
    
    landmass <- loadRData(paste0(dataFolder,"/Spatial/globalLandmass/landmassGlobal.RData"))
    eez <- loadRData(paste0(dataFolder,"/Spatial/EEZ/EEZGlobal.RData"))
    mpa <- loadRData(paste0(dataFolder,"/Spatial/WDPA/WDPAGlobal.RData"))
    
    mainMap <- ggplot() + 
      geom_polygon(data = landmass, aes(long,lat,group=group), colour = "#E6E6E6" , fill="#E6E6E6", size=0.1 ) +
      xlab("Longitude") + 
      ylab("Latitude") + 
      theme_minimal() +
      theme(
        text = element_text(family = "Helvetica", color = "#22211d"),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        panel.grid.major = element_line(color = "black", size = 0.1),
        panel.grid.minor = element_blank(),
        plot.background = element_rect(fill = "#FFFFFF", color = NA), # F3F3F3
        panel.background = element_rect(fill = "#FFFFFF", color = NA), # F3F3F3
        legend.background = element_rect(fill = "#FFFFFF", color = NA),
        panel.border = element_blank() )
    
    # ---------
    
    resultsName <- "Climatic dissimilarity (sigma-dissimilarity)"
    dataRaster <- raster(paste0(resultsFolder,"Raster/regionalClimaticDissimilarity.tif"))
    dataRaster <- as.data.frame(dataRaster,xy=TRUE,na.rm=T)
    colnames(dataRaster) <- c("Lon","Lat","Val")
    myColors <- c("#6FBBE8","#A1ECD8","#F6F9AB","#FCB46D","#B21414")
    
    plot1 <- mainMap +
      geom_tile(data = dataRaster, aes(x=Lon,y=Lat,fill=Val)) +
      scale_fill_gradientn(colours=myColors, na.value='transparent')
    
    pdf(paste0(resultsFolder,"/Figures/",resultsName,".pdf"), width = 10 )
    print(plot1)
    dev.off()
    
    # ---------
    
    resultsName <- "Climatic dissimilarity to best analog (sigma-dissimilarity)"
    dataRaster <- raster(paste0(resultsFolder,"Raster/analogClimaticDissimilarity.tif"))
    dataRaster <- as.data.frame(dataRaster,xy=TRUE,na.rm=T)
    colnames(dataRaster) <- c("Lon","Lat","Val")
    myColors <- c("#6FBBE8","#A1ECD8","#F6F9AB","#FCB46D","#B21414")
    
    plot1 <- mainMap +
      geom_tile(data = dataRaster, aes(x=Lon,y=Lat,fill=Val)) +
      scale_fill_gradientn(colours=myColors, na.value='transparent')
    
    pdf(paste0(resultsFolder,"/Figures/",resultsName,".pdf"), width = 10 )
    print(plot1)
    dev.off()
    
    # ---------
    
    resultsName <- "Time of emergence"
    dataRaster <- raster(paste0(resultsFolder,"Raster/tEmergence.tif"))
    dataRaster <- as.data.frame(dataRaster,xy=TRUE,na.rm=T)
    colnames(dataRaster) <- c("Lon","Lat","Val")
    myColors <- c("#6FBBE8","#A1ECD8","#F6F9AB","#FCB46D","#B21414")
    
    plot1 <- mainMap +
      geom_tile(data = dataRaster, aes(x=Lon,y=Lat,fill=Val)) +
      scale_fill_gradientn(colours=myColors, na.value='transparent')
    
    pdf(paste0(resultsFolder,"/Figures/",resultsName,".pdf"), width = 10 )
    print(plot1)
    dev.off()
    
    # ---------
    
    resultsName <- "Refugia"
    dataRaster <- raster(paste0(resultsFolder,"Raster/refugia.tif"))
    dataRaster <- as.data.frame(dataRaster,xy=TRUE,na.rm=T)
    colnames(dataRaster) <- c("Lon","Lat","Val")
    myColors <- c("#A41F1F")
    
    plot1 <- mainMap +
      geom_tile(data = dataRaster, aes(x=Lon,y=Lat,fill=Val)) +
      scale_fill_gradientn(colours=myColors, na.value='transparent',limits=c(1,1))
    
    pdf(paste0(resultsFolder,"/Figures/",resultsName,".pdf"), width = 10 )
    print(plot1)
    dev.off()
    
    # ---------
    
    resultsName <- "Analog distance"
    dataRaster <- raster(paste0(resultsFolder,"Raster/analogDistance.tif"))
    dataRaster <- as.data.frame(dataRaster,xy=TRUE,na.rm=T)
    colnames(dataRaster) <- c("Lon","Lat","Val")
    myColors <- c("#6FBBE8","#A1ECD8","#F6F9AB","#FCB46D","#B21414")
    
    plot1 <- mainMap +
      geom_tile(data = dataRaster, aes(x=Lon,y=Lat,fill=Val)) +
      scale_fill_gradientn(colours=myColors, na.value='transparent')
    
    pdf(paste0(resultsFolder,"/Figures/",resultsName,".pdf"), width = 10 )
    print(plot1)
    dev.off()
    
    # ---------
    
    resultsName <- "Corridor"
    dataRaster <- raster(paste0(resultsFolder,"Raster/corridor.tif"))
    dataRaster <- as.data.frame(dataRaster,xy=TRUE,na.rm=T)
    colnames(dataRaster) <- c("Lon","Lat","Val")
    myColors <- c("#6FBBE8","#A1ECD8","#F6F9AB","#FCB46D","#B21414")
    
    plot1 <- mainMap +
      geom_tile(data = dataRaster, aes(x=Lon,y=Lat,fill=Val)) +
      scale_fill_gradientn(colours=myColors, na.value='transparent')
    
    pdf(paste0(resultsFolder,"/Figures/",resultsName,".pdf"), width = 10 )
    print(plot1)
    dev.off()
    
    # ---------
    
    resultsName <- "Corridor 95"
    dataRaster <- raster(paste0(resultsFolder,"Raster/corridorQuantile95.tif"))
    dataRaster <- as.data.frame(dataRaster,xy=TRUE,na.rm=T)
    colnames(dataRaster) <- c("Lon","Lat","Val")
    myColors <- c("#A41F1F")
    
    plot1 <- mainMap +
      geom_tile(data = dataRaster, aes(x=Lon,y=Lat,fill=Val)) +
      scale_fill_gradientn(colours=myColors, na.value='transparent',limits=c(1,1))
    
    pdf(paste0(resultsFolder,"/Figures/",resultsName,".pdf"), width = 10 )
    print(plot1)
    dev.off()
    
  }

}
