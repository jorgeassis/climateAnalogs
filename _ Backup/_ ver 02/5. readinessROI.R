# ----------------------------------------------------------------------
# ----------------------------------------------------------------------
#
#
# ----------------------------------------------------------------------

closeAllConnections()
rm(list=(ls()[ls()!="v"]))
gc(reset=TRUE)
source("/Volumes/Jellyfish/Dropbox/theMarineDataScientist/gitRepositories/climaticAnalogs/mainFunctions.R")
source("mainFunctions.R")

# ---------

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

scenario <- "ssp585" # ssp119 ssp585
ROIName <- "Global"

dataFolder <- "../Data/"
resultsSubFolder <- paste0("",scenario,"/")
resultsFolder <- paste0("../Results/",resultsSubFolder)
polygonAnalysesNames <- "readinessMPAinsideEEZ"

analogResults <- loadRData(paste0(resultsFolder,"/novelClimates",ROIName,".RData"))
load(paste0(resultsFolder,"/dataStructure.RData"))

nCores <- 32
iteractions <- 999

# --------------

landmass <- loadRData(paste0(dataFolder,"/Spatial/globalLandmass/landmassGlobal.RData"))
eez <- loadRData(paste0(dataFolder,"/Spatial/EEZ/EEZGlobal.RData"))
mpa <- loadRData(paste0(dataFolder,"/Spatial/WDPA/WDPAGlobal.RData"))

# --------------

polygonAnalyses <- mpa
polygonAnalysesWithin <- eez

# -----------------------------------------------
# Get climate readiness

novelClimateSigma <- rasterFromXYZ(analogResults[,c("x","y","novelClimateSigma")])

tik <- Sys.time()
Cluster <- makeCluster( nCores )
registerDoParallel( Cluster )

parallelProcess <- foreach(i = 1:nrow(polygonAnalyses), .combine=rbind, .verbose=FALSE, .packages=c("sf","rgeos","raster","sp")) %dopar% {
  
  sf::sf_use_s2(FALSE)
  
  polygonAnalyses.i <- polygonAnalyses[i,]
  polygonAnalyses.i.Name <- polygonAnalyses.i$Name
  polygonAnalyses.i.WithinName <- polygonAnalyses.i$EEZ
  
  polygonAnalysesWithin.i <- polygonAnalysesWithin[polygonAnalysesWithin$EEZ == polygonAnalyses.i.WithinName,]
  
  if(nrow(polygonAnalysesWithin.i) == 0) { stop("Wrong EEZ name")}
  
  novelClimateSigma.i <- crop(novelClimateSigma,polygonAnalysesWithin.i)
  novelClimateSigma.i <- mask(novelClimateSigma.i,polygonAnalysesWithin.i)
  
  polygonAnalyses.i.val <- unlist(raster::extract(novelClimateSigma.i,polygonAnalyses.i,small=TRUE))
  polygonAnalysesWithin.i.val <- unlist(raster::extract(novelClimateSigma.i,polygonAnalysesWithin.i,small=TRUE))
  
  polygonAnalyses.i.val <- polygonAnalyses.i.val[!is.na(polygonAnalyses.i.val)]
  polygonAnalysesWithin.i.val <- polygonAnalysesWithin.i.val[!is.na(polygonAnalysesWithin.i.val)]
  
  bufferPolygon <- 0
  while( length(polygonAnalyses.i.val) == 0 ) { 
    bufferPolygon <- bufferPolygon + 0.1; if(bufferPolygon > 5) { break }
    polygonAnalyses.i.val <- unlist(raster::extract(novelClimateSigma.i,st_buffer(polygonAnalyses.i,dist=bufferPolygon) ,small=TRUE))
    polygonAnalyses.i.val <- polygonAnalyses.i.val[!is.na(polygonAnalyses.i.val)]
  }
  
  if( length(polygonAnalyses.i.val) != 0 & length(polygonAnalysesWithin.i.val) != 0 ) { 
    
    ready <- numeric(iteractions)
    sampledRegion <- xyFromCell(novelClimateSigma.i,Which(!is.na(novelClimateSigma.i), cells=TRUE))
    
    for(int in 1:iteractions){
      sampledRegion.c <- sampledRegion[sample(1:nrow(sampledRegion),1),]
      sampledRegion.c <- sort(spDistsN1( as.matrix(sampledRegion) , matrix(sampledRegion.c , nrow=1 )), index.return=T)$ix[1:length(polygonAnalyses.i.val)]
      ready[int] <- mean(unlist(raster::extract(novelClimateSigma.i,matrix(sampledRegion[sampledRegion.c,], ncol=2),small=TRUE)), na.rm=T) < mean(polygonAnalyses.i.val, na.rm=T)
    }
    
    resultsDF <- data.frame( polygonName=polygonAnalyses.i.Name, polygonArea=as.numeric(st_area(polygonAnalyses.i)), withinName=polygonAnalyses.i.WithinName, polygonAverageSigma=mean(polygonAnalyses.i.val,na.rm=T), polygonReadynessPvalue = sum(ready > 0, na.rm=T) / sum( !is.na(ready) , na.rm=T)) 
    
  }
  
  if( length(polygonAnalyses.i.val) == 0 | length(polygonAnalysesWithin.i.val) == 0 ) { 
    
    resultsDF <- data.frame( polygonName=polygonAnalyses.i.Name, polygonArea=as.numeric(st_area(polygonAnalyses.i)), withinName=polygonAnalyses.i.WithinName, polygonAverageSigma=NA, polygonReadynessPvalue = 0 ) 
    
  }
  
  return(resultsDF)
  
}

stopCluster(Cluster); rm(Cluster)
closeAllConnections()
Sys.time() - tik

# --------------

polygonAnalyses$polygonArea <- parallelProcess$polygonArea
polygonAnalyses$withinName <- parallelProcess$withinName
polygonAnalyses$polygonAverageSigma <- parallelProcess$polygonAverageSigma
polygonAnalyses$polygonReadynessPvalue <- parallelProcess$polygonReadynessPvalue

save( polygonAnalyses, file=paste0(resultsFolder,"/",polygonAnalysesNames,".RData") )

polygonAnalysesDF <- as.data.frame(polygonAnalyses)[,-which(names(polygonAnalyses) == "geometry")]
for( c in 1:ncol(polygonAnalysesDF)) { polygonAnalysesDF[,c] <- gsub(",",";",polygonAnalysesDF[,c]) }
write.csv(polygonAnalysesDF,file=paste0(resultsFolder,"/",polygonAnalysesNames,".csv"), row.names = FALSE)

# -----------------------
# Histograms about no novelty, moderate, extreme


# -----------------------
# Aggregate data per Within areas
# Determine the proportion of area with better climates than in polygons 

novelClimateSigma <- rasterFromXYZ(analogResults[,c("x","y","novelClimateSigma")])

polygonAnalysesWithin$area <- NA
polygonAnalysesWithin$averageSigma <- NA
polygonAnalysesWithin$sdSigma <- NA
polygonAnalysesWithin$polygonNumber <- NA
polygonAnalysesWithin$polygonAverageSigma <- NA
polygonAnalysesWithin$polygonAreaTotal <- NA
polygonAnalysesWithin$polygonAreaProportion <- NA
polygonAnalysesWithin$polygonNumberModerate <- NA
polygonAnalysesWithin$polygonNumberModerateProportion <- NA
polygonAnalysesWithin$polygonNumberExtreme <- NA
polygonAnalysesWithin$polygonNumberExtremeProportion <- NA
polygonAnalysesWithin$polygonNumberReady <- NA
polygonAnalysesWithin$polygonNumberReadyProportion <- NA
polygonAnalysesWithin$polygonAreaReady <- NA
polygonAnalysesWithin$polygonAreaReadyProportion <- NA
polygonAnalysesWithin$betterClimateThanPolygon <- NA

for(i in 1:nrow(polygonAnalysesWithin)) {
  
  cat(i, "out of" , nrow(polygonAnalysesWithin) , "\n")
  polygonAnalysesWithin.i <- polygonAnalysesWithin[i,]
  polygonAnalyses.i <- polygonAnalyses[polygonAnalyses$withinName == polygonAnalysesWithin[i,]$EEZ,]
  
  if(nrow(polygonAnalyses.i) == 0) { next }
  
  # -----------------------
  
  analogResults <- unlist(raster::extract(novelClimateSigma, polygonAnalysesWithin[i,] ,small=TRUE))
  analogResults <- analogResults[!is.na(analogResults)]
  
  bufferPolygon <- 0
  while( length(analogResults) == 0 ) {
    sf::sf_use_s2(FALSE)
    bufferPolygon <- bufferPolygon + 0.1
    analogResults <- unlist(raster::extract(novelClimateSigma, st_buffer(polygonAnalysesWithin[i,],dist=bufferPolygon) ,small=TRUE))
    analogResults <- analogResults[!is.na(analogResults)]
  }
  
  analogResultsPoly <- unlist(raster::extract(novelClimateSigma,polygonAnalyses.i,small=TRUE))
  analogResultsPoly <- analogResultsPoly[!is.na(analogResultsPoly)]
  analogResultsPoly <- ifelse( length(analogResultsPoly) > length(analogResults) , rep(analogResultsPoly,length(analogResults)) , analogResultsPoly )
  if( length(analogResults) == 1 & length(analogResultsPoly) == 1) { analogResults <- rep( analogResults, 2); analogResultsPoly <- rep(analogResultsPoly,2) }
  
  betterClimatesProportion <- sum(apply(sapply(analogResultsPoly, function(x) analogResults < x),1,sum,na.rm=T) != 0) / length(analogResults)
  
  # -----------------------
  
  polygonAnalysesWithin[i,"area"] <- as.numeric(st_area(polygonAnalysesWithin[i,])) # m2
  polygonAnalysesWithin[i,"averageSigma"] <- mean(analogResults, na.rm=T)
  polygonAnalysesWithin[i,"sdSigma"] <- sd(analogResults, na.rm=T)
  polygonAnalysesWithin[i,"polygonAverageSigma"] <- mean(polygonAnalyses.i$polygonAverageSigma)
  polygonAnalysesWithin[i,"polygonNumber"] <- as.numeric(nrow(polygonAnalyses.i))
  polygonAnalysesWithin[i,"polygonAreaTotal"] <- sum(polygonAnalyses.i$polygonArea)
  
  polygonAnalysesWithin[i,"polygonAreaProportion"] <- sum(polygonAnalyses.i$polygonArea) / polygonAnalysesWithin[i,]$area
  
  polygonAnalysesWithin[i,"polygonNumberModerate"] <- sum(polygonAnalyses.i$polygonAverageSigma >= 2 & polygonAnalyses.i$polygonAverageSigma < 4)
  polygonAnalysesWithin[i,"polygonNumberModerateProportion"] <- sum(polygonAnalyses.i$polygonAverageSigma >=2 & polygonAnalyses.i$polygonAverageSigma < 4) / nrow(polygonAnalyses.i)
  polygonAnalysesWithin[i,"polygonNumberExtreme"] <- sum(polygonAnalyses.i$polygonAverageSigma >= 4)
  polygonAnalysesWithin[i,"polygonNumberExtremeProportion"] <- sum(polygonAnalyses.i$polygonAverageSigma >= 4) / nrow(polygonAnalyses.i)
  
  polygonAnalysesWithin[i,"polygonNumberReady"] <- sum(polygonAnalyses.i$polygonReadynessPvalue > 0.95)
  polygonAnalysesWithin[i,"polygonNumberReadyProportion"] <- sum(polygonAnalyses.i$polygonReadynessPvalue > 0.95) /  nrow(polygonAnalyses.i)
  polygonAnalysesWithin[i,"polygonAreaReady"] <- sum(polygonAnalyses.i[polygonAnalyses.i$polygonReadynessPvalue > 0.95,]$polygonArea)
  polygonAnalysesWithin[i,"polygonAreaReadyProportion"] <- sum(polygonAnalyses.i[polygonAnalyses.i$polygonReadynessPvalue > 0.95,]$polygonArea) / sum(polygonAnalyses.i[,]$polygonArea)
  polygonAnalysesWithin[i,"betterClimateThanPolygon"] <- betterClimatesProportion
  
}

save( polygonAnalysesWithin, file=paste0(resultsFolder,"/",polygonAnalysesNames,"Agg.RData") )

polygonAnalysesWithinDF <- as.data.frame(polygonAnalysesWithin)[,-which(names(polygonAnalysesWithin) == "geometry")]
for( c in 1:ncol(polygonAnalysesWithinDF)) { polygonAnalysesWithinDF[,c] <- gsub(",",";",polygonAnalysesWithinDF[,c]) }
write.csv(polygonAnalysesWithinDF,file=paste0(resultsFolder,"/",polygonAnalysesNames,"Agg.csv"), row.names = FALSE)

# ----------------------------------------------------------------------
# ----------------------------------------------------------------------
# Make maps

# polygon scale

polygonAnalyses <- as(polygonAnalyses,"Spatial")
polygonAnalyses$id <- polygonAnalyses$Name
polygonAnalyses$ready <- ifelse(polygonAnalyses$polygonReadynessPvalue > 0.95 ,1,0)
polygonAnalyses$moderate <- ifelse( polygonAnalyses$polygonAverageSigma >= 2 & polygonAnalyses$polygonAverageSigma < 4 ,1,0)
polygonAnalyses$extreme <- ifelse(polygonAnalyses$polygonAverageSigma >= 4,1,0)

polygonAnalyses.f <- fortify(polygonAnalyses)
polygonAnalyses.f <- merge(polygonAnalyses.f, polygonAnalyses@data,by = "id")
head(polygonAnalyses.f)

names(polygonAnalyses.f)
columnName <- "polygonAverageSigma"
legendName <- "Novel climates" 

# "Novel climates" : polygonAverageSigma
# "Climate ready MPA" : ready
# "MPA exposed to moderate novel climate" : moderate
# "MPA exposed to extreme novel climate" : extreme

fig1 <- mainMap +
  geom_polygon(data=polygonAnalyses.f, aes(x = long, y = lat, group = group), size = 0.1, fill= c("#8CD1DF","#B30C0C")[polygonAnalyses.f[,columnName] + 1]) +
  labs(fill=legendName) +
  geom_polygon(data = landmass, aes(long,lat,group=group), colour = "#E6E6E6" , fill="#E6E6E6", size=0.1 )

pdf(paste0(resultsFolder,"/",polygonAnalysesNames," ",legendName,".pdf"), width = 10 )
print(fig1)
dev.off()

# ----------------
# ----------------
# Within scale

polygonAnalysesWithin <- as(polygonAnalysesWithin,"Spatial")
polygonAnalysesWithin$id <- polygonAnalysesWithin$EEZ
polygonAnalysesWithin.f <- fortify(polygonAnalysesWithin)
polygonAnalysesWithin.f <- merge(polygonAnalysesWithin.f, polygonAnalysesWithin@data,by = "id")
head(polygonAnalysesWithin.f)

names(polygonAnalysesWithin.f)
columnName <- "polygonNumberReadyProportion"
legendName <- "Proportion of climate ready MPA (number)" 

# "Proportion of climate ready MPA (number)" : polygonNumberReadyProportion
# "Proportion of climate ready MPA (area)" : polygonAreaReadyProportion
# "Proportion of better climates than in MPA"  : betterClimateThanProportion

fig1 <- mainMap +
  geom_polygon(data=polygonAnalysesWithin.f, aes(x = long, y = lat, group = group, fill = as.numeric(polygonNumberReadyProportion)), size = 0.1) +
  labs(fill=legendName) +
  scale_fill_gradient2(low = "#8CD1DF",mid = "#EBE662", high = "#B30C0C",
                       midpoint = min(polygonAnalysesWithin.f[,columnName], na.rm=T) + ((max(polygonAnalysesWithin.f[,columnName], na.rm=T) - min(polygonAnalysesWithin.f[,columnName], na.rm=T)) / 2), na.value = NA,
                       aesthetics = "fill") +
  geom_polygon(data = landmass, aes(long,lat,group=group), colour = "#E6E6E6" , fill="#E6E6E6", size=0.1 )


pdf(paste0(resultsFolder,"/",polygonAnalysesNames," ",legendName,".pdf"), width = 10 )
print(fig1)
dev.off()

# ----------------------------------------------------------------------
# ----------------------------------------------------------------------