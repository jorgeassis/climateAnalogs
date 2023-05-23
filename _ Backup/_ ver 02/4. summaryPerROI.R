# ----------------------------------------------------------------------
# ----------------------------------------------------------------------
#
#
# ----------------------------------------------------------------------

closeAllConnections()
rm(list=(ls()[ls()!="v"]))
gc(reset=TRUE)
source("/Volumes/Jellyfish/Dropbox/theMarineDataScientist/gitRepositories/climaticAnalogs/mainFunctions.R")

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

# ----------------------------------------------------------------------
# ----------------------------------------------------------------------

scenario <- "ssp119" # ssp119 ssp585
ROIName <- "Global"

dataFolder <- "../Data/"
resultsSubFolder <- paste0("",scenario,"/")
resultsFolder <- paste0("../Results/",resultsSubFolder)

landmass <- loadRData(paste0(dataFolder,"/Spatial/globalLandmass/landmassGlobal.RData"))
eez <- loadRData(paste0(dataFolder,"/Spatial/EEZ/EEZGlobal.RData"))
mpa <- loadRData(paste0(dataFolder,"/Spatial/WDPA/WDPAGlobal.RData"))

predictors <- colnames(loadRData(paste0(resultsFolder,"/novelClimatePredictorContrib.RData")))[-1]
predictorsReduced <- c("Oxygen","Max. Temp.","pH","Productivity")

load(paste0(resultsFolder,"/dataStructure.RData"))

# ----------------------------------------------------------------------
# ----------------------------------------------------------------------

polygonDFNames <- "analysesPerMPA" # analysesPerEEZ analysesPerMPA
polygonDF <- mpa # eez mpa

# --------------------
# --------------------

if( class(polygonDF)[1] == "sf") { polygonDF <- as(polygonDF, "Spatial") }

polygonDF$novelClimateSigmaMean <- NA
polygonDF$novelClimateSigmaMax <- NA
polygonDF$disappearClimateSigmaMean <- NA
polygonDF$disappearClimateSigmaMax <- NA
polygonDF$highNovelClimateProportion <- NA
polygonDF$analogDistanceMean <- NA
polygonDF$refugiaProportion <- NA

analogResults <- loadRData(paste0(resultsFolder,"/novelClimates",ROIName,".RData"))

# --------------------

novelClimateSigma <- rasterFromXYZ(analogResults[,c("x","y","novelClimateSigma")])
disappearClimateSigma <- rasterFromXYZ(analogResults[,c("x","y","disappearClimateSigma")])
analogDistance <- rasterFromXYZ(analogResults[,c("x","y","analogDistance")])
refugia <- rasterFromXYZ(analogResults[,c("x","y","refugia")])

for( i in 1:nrow(polygonDF) ) {
  
  cat("\014")  
  cat("# --------------------")
  cat("\n")  
  cat("# Progress:",i,"in",length(polygonDF) )
  cat("\n")  
  cat("# --------------------")
  
  polygonDF.i <- polygonDF[i,]
  bufferPolygon <- 0
  
  while( sum(!is.na(unlist(raster::extract(novelClimateSigma,polygonDF.i,small=TRUE)))) == 0 ) { 
    bufferPolygon <- bufferPolygon + 0.1
    polygonDF.i <- gBuffer(polygonDF.i,width=bufferPolygon) 
  }
  
  novelClimateSigma.i <- unlist(raster::extract(novelClimateSigma,polygonDF.i,small=TRUE)); novelClimateSigma.i <- novelClimateSigma.i[!is.na(novelClimateSigma.i)]
  disappearClimateSigma.i <- unlist(raster::extract(disappearClimateSigma,polygonDF.i,small=TRUE)); disappearClimateSigma.i <- disappearClimateSigma.i[!is.na(disappearClimateSigma.i)]
  analogDistance.i <- unlist(raster::extract(analogDistance,polygonDF.i,small=TRUE)); analogDistance.i <- analogDistance.i[!is.na(analogDistance.i)]
  refugia.i <- unlist(raster::extract(refugia,polygonDF.i,small=TRUE)); refugia.i <- refugia.i[!is.na(refugia.i)]
  
  polygonDF[i,"novelClimateSigmaMean"] <- round(mean(novelClimateSigma.i), digits = 2)
  polygonDF[i,"novelClimateSigmaMax"] <- round(max(novelClimateSigma.i), digits = 2)
  polygonDF[i,"disappearClimateSigmaMean"] <- round(mean(disappearClimateSigma.i), digits = 2)
  polygonDF[i,"disappearClimateSigmaMax"] <- round(max(disappearClimateSigma.i), digits = 2)
  polygonDF[i,"highNovelClimateProportion"] <- round(sum(novelClimateSigma.i >= 4) / length(novelClimateSigma.i), digits = 2)
  polygonDF[i,"analogDistanceMean"] <- round(mean(analogDistance.i), digits = 2)
  polygonDF[i,"refugiaProportion"] <- round(sum(refugia.i) / length(refugia.i), digits = 2)

}

write.csv(as.data.frame(polygonDF),file=paste0(resultsFolder,"/",polygonDFNames,".csv"), row.names = FALSE)

# ----------------------------------------------------------------------
# ----------------------------------------------------------------------
# Make maps

mainMap <- ggplot() + 
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
    plot.background = element_rect(fill = "#FFFFFF", color = NA),
    panel.background = element_rect(fill = "#FFFFFF", color = NA),
    legend.background = element_rect(fill = "#FFFFFF", color = NA),
    panel.border = element_blank() )

polygonDF$id <- 1:nrow(polygonDF)
polygonDF.f <- fortify(polygonDF)
polygonDF.f <- merge(polygonDF.f, polygonDF@data,by = "id")

names(polygonDF.f)
columnName <- "refugiaProportion" # novelClimateSigmaMean highNovelClimateProportion refugiaProportion
legendName <- "Proportion of refugia" # "Proportion of high novel climate"  "Proportion of refugia" "Mean novel climate"

fig1 <- mainMap +
  geom_polygon(data=polygonDF.f, aes(x = long, y = lat, group = group, fill = as.numeric(refugiaProportion)), size = 0.1) +
  labs(fill=legendName) +
  scale_fill_gradient2(low = "#8CD1DF",mid = "#EBE662", high = "#B30C0C",
                       midpoint = min(polygonDF.f[,columnName], na.rm=T) + ((max(polygonDF.f[,columnName], na.rm=T) - min(polygonDF.f[,columnName], na.rm=T)) / 2), na.value = NA,
                       aesthetics = "fill") +
  geom_polygon(data = landmass, aes(long,lat,group=group), colour = "#E6E6E6" , fill="#E6E6E6", size=0.1 )
  

pdf(paste0(resultsFolder,"/",polygonDFNames," ",legendName,".pdf"), width = 10 )
print(fig1)
dev.off()
