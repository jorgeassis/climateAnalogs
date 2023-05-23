# ----------------------------------------------------------------------
# ----------------------------------------------------------------------
#
#
# ----------------------------------------------------------------------

closeAllConnections()
rm(list=(ls()[ls()!="v"]))
gc(reset=TRUE)
source("mainFunctions.R")

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

scenario <- "ssp119" #  ssp585

dataFolder <- "../Data/"
resultsSubFolder <- paste0("finalRun/",scenario,"/")
resultsFolder <- paste0("../Results/",resultsSubFolder)

landmass <- loadRData(paste0(dataFolder,"/Spatial/globalLandmass/landmassGlobal.RData"))
eez <- loadRData(paste0(dataFolder,"/Spatial/EEZ/EEZGlobal.RData"))
mpa <- loadRData(paste0(dataFolder,"/Spatial/WDPA/WDPAGlobalSimpEEZAgg.RData"))

predictors <- colnames(loadRData(paste0(resultsFolder,"/dissimilarityPredictors.RData")))[-1]

load(paste0(resultsFolder,"/dataStructure.RData"))
load(paste0(resultsFolder,"/dataStructureRaster.RData"))

# ----------------------------------------------------------------------
# ----------------------------------------------------------------------

polygonDFNames <- "analysesPerEEZ"
polygonDF <- eez

# --------------------
# --------------------

polygonDF <- polygonDF[,c("EEZ","IHO_SEA","SOVEREIGN1","ISO_SOV1","AREA_KM2")]
names(polygonDF) <- c("EEZ","Region","Country","CountryCode","AreaKm2")

polygonDF$meanLocalClimaticDissimilarity <- NA
polygonDF$meanAnalogClimaticDissimilarity <- NA
polygonDF$noveltyAreaProportion <- NA
polygonDF$meanAnalogDistance <- NA
polygonDF$meanAnalogVerticalDistance <- NA
polygonDF$corridorAreaProportion <- NA
polygonDF$meanTimeEmergence <- NA
polygonDF$refugiaAreaProportion <- NA

analogResults <- loadRData(paste0(resultsFolder,"/analogResults.RData"))

structureCoordinates <- analogResults[,c("x","y")]
coordinates(structureCoordinates) <- ~x+y
crs(structureCoordinates) <- crs(polygonDF)

for( i in 1:nrow(polygonDF) ) {
  
  cat("\014")  
  cat("# --------------------")
  cat("\n")  
  cat("# Progress:",i,"in",length(polygonDF) )
  cat("\n")  
  cat("# --------------------")
  
  polygonDF.i <- polygonDF[i,]
  
  analogsOverPolygon <- which(!is.na(over(structureCoordinates,polygonDF.i)[,1]))
  
  if( length(analogsOverPolygon) == 0 ) { 
    
    analogsOverPolygon <- gDistance(structureCoordinates, polygonDF.i,byid=TRUE)
    analogsOverPolygon <- which.min(analogsOverPolygon)
    
  }
  
  analogResults.i <- analogResults[analogsOverPolygon,]

  nCells <- nrow(analogResults.i)
  
  polygonDF[i,"meanLocalClimaticDissimilarity"] <- round(mean(analogResults.i$localClimaticDissimilarity,na.rm=T), digits = 2)
  polygonDF[i,"meanAnalogClimaticDissimilarity"] <- round(mean(analogResults.i$analogClimaticDissimilarity,na.rm=T), digits = 2)
  polygonDF[i,"noveltyAreaProportion"] <- round(sum( analogResults.i$localClimaticDissimilarity > 4 ) / nCells, digits = 2)
  polygonDF[i,"meanAnalogDistance"] <- round(mean(analogResults.i$analogDistance,na.rm=T), digits = 2)
  polygonDF[i,"meanAnalogVerticalDistance"] <- round(mean(analogResults.i$analogDistanceVert,na.rm=T), digits = 2)
  polygonDF[i,"corridorAreaProportion"] <- round(sum(analogResults.i$corridorQuantile95) / nCells, digits = 2)
  polygonDF[i,"meanTimeEmergence"] <- round(mean(analogResults.i$tEmergence,na.rm=T), digits = 2)
  polygonDF[i,"refugiaAreaProportion"] <- round(sum(analogResults.i$absoluteRefugia) / nCells, digits = 2)

}

which(is.na(polygonDF$meanLocalClimaticDissimilarity))
write.csv(as.data.frame(polygonDF),file=paste0(resultsFolder,"/",polygonDFNames,".csv"), row.names = FALSE)

# ----------------------------

library("RColorBrewer")
display.brewer.pal(n = 4, name = 'YlOrRd')
pallete <- brewer.pal(n = 4, name = "YlOrRd")

assignGroup <- function(vector1,vector2,nclass) {
  min <- min(c(vector1,vector2),na.rm=T)
  max <- max(c(vector1,vector2),na.rm=T)
  list(as.numeric(cut(vector1, breaks=c(seq(min,max,length.out=nclass+1)[-length(seq(min,max,length.out=nclass+1))],max), labels=1:nclass, include.lowest = TRUE)),
       as.numeric(cut(vector2, breaks=c(seq(min,max,length.out=nclass+1)[-length(seq(min,max,length.out=nclass+1))],max), labels=1:nclass, include.lowest = TRUE)))
}

plotPolygonDF <- data.frame(name=c(polygonDF$EEZ,polygonDF$EEZ),
                            index=c(1:nrow(polygonDF),1:nrow(polygonDF)),
                            position=c(rep(1,nrow(polygonDF)),rep(2,nrow(polygonDF))),
                            color=c(pallete[assignGroup(polygonDF$meanLocalClimaticDissimilarity,polygonDF$meanLocalClimaticDissimilarity,4)[[1]]],
                                    pallete[assignGroup(polygonDF$meanLocalClimaticDissimilarity,polygonDF$meanLocalClimaticDissimilarity,4)[[2]]]))

ggplot() + geom_point( data=plotPolygonDF, aes(x = position, y = index) ,colour=plotPolygonDF$color)

# ----------------------------------------------------------------------
# ----------------------------------------------------------------------
# Make maps

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
    plot.background = element_rect(fill = "#FFFFFF", color = NA),
    panel.background = element_rect(fill = "#FFFFFF", color = NA),
    legend.background = element_rect(fill = "#FFFFFF", color = NA),
    panel.border = element_blank() )

# in order to plot polygons, first fortify the data
polygonDF@data$id <- rownames(polygonDF@data)
polygonDFDF <- fortify(polygonDF, region = "id")
polygonDFDF <- merge(polygonDFDF, polygonDF@data,by = "id")

names(polygonDF)
columnName <- "corridorAreaProportionSSP585"
legendName <- "Proportion of novel climate"

fig1 <- mainMap +
  geom_polygon(data=polygonDFDF, aes(x = long, y = lat, group = group,fill = get(columnName)), size = 0.1) +
  labs(fill=legendName) +
  scale_fill_gradient2(low = "#8CD1DF",mid = "#EBE662", high = "#B30C0C",
                       midpoint = min(polygonDFDF[,columnName], na.rm=T) + ((max(polygonDFDF[,columnName], na.rm=T) - min(polygonDFDF[,columnName], na.rm=T)) / 2), na.value = NA,
                       aesthetics = "fill")

pdf(paste0(resultsFolder,"/",polygonDFNames," ",legendName,".pdf"), width = 10 )
print(fig1)
dev.off()
