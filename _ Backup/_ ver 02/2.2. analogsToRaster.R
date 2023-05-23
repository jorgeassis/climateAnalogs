# ----------------------------------------------------------------------
# ----------------------------------------------------------------------
#
#
# ----------------------------------------------------------------------

setwd("/Volumes/Jellyfish/Dropbox/Manuscripts/Global connectivity corridors and refugia of climatic analogs for marine biodiversity/Code")

closeAllConnections()
rm(list=(ls()[ls()!="v"]))
gc(reset=TRUE)
source("/Volumes/Jellyfish/Dropbox/theMarineDataScientist/gitRepositories/climaticAnalogs/mainFunctions.R")

# ----------------------------------------------------------------------
# ----------------------------------------------------------------------

scenario <- "ssp119" # ssp119 ssp585
ROIName <- "Global"

dataFolder <- "../Data/"
resultsSubFolder <- paste0("",scenario,"/")
resultsFolder <- paste0("../Results/",resultsSubFolder)

# -----------

load(file=paste0(dataFolder,"/Spatial/globalLandmass/landmassGlobal.RData"))
load(file=paste0(dataFolder,"/Spatial/EEZ/EEZGlobal.RData"))
load(file=paste0(dataFolder,"/Spatial/WDPA/WDPAGlobal.RData"))

load(paste0(resultsFolder,"/dataStructure.RData"))
load(paste0(resultsFolder,"/novelClimates",ROIName,".RData"))

load(paste0(resultsFolder,"/novelClimatePredictorContrib.RData"))
load(paste0(resultsFolder,"/disappearClimatePredictorContrib.RData"))

# ------
# ------

if(! dir.exists(paste0(resultsFolder,"/Figures/"))) { dir.create(paste0(resultsFolder,"/Figures/"), recursive = T) }
if(! dir.exists(paste0(resultsFolder,"/Raster/"))) { dir.create(paste0(resultsFolder,"/Raster/"), recursive = T) }

# ------
# ------

novelClimateSigma <- rasterFromXYZ(dataStructureResult[,c("x","y","novelClimateSigma")])
writeRaster(novelClimateSigma,file=paste0(resultsFolder,"/Raster/novelClimateSigma.tif"), format="GTiff", overwrite=TRUE)

disappearClimateSigma <- rasterFromXYZ(dataStructureResult[,c("x","y","disappearClimateSigma")])
writeRaster(disappearClimateSigma,file=paste0(resultsFolder,"/Raster/disappearClimateSigma.tif"), format="GTiff", overwrite=TRUE)

refugia <- rasterFromXYZ(dataStructureResult[,c("x","y","refugia")])
writeRaster(refugia,file=paste0(resultsFolder,"/Raster/refugia.tif"), format="GTiff", overwrite=TRUE)

analogDistance <- rasterFromXYZ(dataStructureResult[,c("x","y","analogDistance")])
writeRaster(analogDistance,file=paste0(resultsFolder,"/Raster/analogDistance.tif"), format="GTiff", overwrite=TRUE)

# ----------------------------------------------------------------------
# ----------------------------------------------------------------------
# Map climatic analogs per cell

landmass <- ne_countries(scale = 10, returnclass = "sp")

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

resultsName <- "Novel climates (sigma-dissimilarity)"
dataRaster <- rasterFromXYZ(dataStructureResult[,c("x","y","novelClimateSigma")])
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

resultsName <- "Disappearing climates (sigma-dissimilarity)"
dataRaster <- rasterFromXYZ(dataStructureResult[,c("x","y","disappearClimateSigma")])
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
dataRaster <- rasterFromXYZ(dataStructureResult[,c("x","y","refugia")])
dataRaster <- as.data.frame(dataRaster,xy=TRUE,na.rm=T)
colnames(dataRaster) <- c("Lon","Lat","Val")
myColors <- c("#FCB46D")

plot1 <- mainMap +
  geom_tile(data = dataRaster, aes(x=Lon,y=Lat,fill=Val)) +
  scale_fill_gradientn(colours=myColors, na.value='transparent',limits=c(1,1)) + theme(legend.position="none")

pdf(paste0(resultsFolder,"/Figures/",resultsName,".pdf"), width = 10 )
print(plot1)
dev.off()

# ---------

resultsName <- "Analog distance"
dataRaster <- rasterFromXYZ(dataStructureResult[,c("x","y","analogDistance")])
dataRaster <- as.data.frame(dataRaster,xy=TRUE,na.rm=T)
colnames(dataRaster) <- c("Lon","Lat","Val")
myColors <- c("#6FBBE8","#A1ECD8","#F6F9AB","#FCB46D","#B21414")

plot1 <- mainMap +
  geom_tile(data = dataRaster, aes(x=Lon,y=Lat,fill=Val)) +
  scale_fill_gradientn(colours=myColors, na.value='transparent')

pdf(paste0(resultsFolder,"/Figures/",resultsName,".pdf"), width = 10 )
print(plot1)
dev.off()

# --------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------