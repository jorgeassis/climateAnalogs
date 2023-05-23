# ----------------------------------------------------------------------
# ----------------------------------------------------------------------
#
#
# ----------------------------------------------------------------------

closeAllConnections()
rm(list=(ls()[ls()!="v"]))
gc(reset=TRUE)
source("mainFunctions.R")
  
library(stringdist)
library(parallel)

# ----------------------------------------------------------------------
# ----------------------------------------------------------------------

dataFolder <- "../Data/"
dumpFolder <- "../Data/" 
regionExtent <- extent(-180,180,-90,90) # extent(-30,65,-45,43) 

# -----------------------------------------
# -----------------------------------------
# Landmass

world <- ne_countries(scale = "large", returnclass = "sp")
landmass <- world
save(landmass,file=paste0(dumpFolder,"/Spatial/globalLandmass/landmassGlobal.RData"))

landmass <- crop(world,regionExtent)

landmassNames <- landmass$name
landmass <- world[world$name %in% landmassNames,]
save(landmass,file=paste0(dumpFolder,"/Spatial/globalLandmass/landmassRegion.RData"))
plot(landmass)

# -----------------------------------------
# -----------------------------------------
# EEZ boundaries shapefile from marineregions.org

EEZ <- shapefile(paste0(dataFolder,"/Spatial/EEZ/Intersect_EEZ_IHO_v4_2020.shp"))
EEZDF <- as.data.frame(EEZ)
EEZ <- gUnaryUnion(EEZ, id = EEZ@data$EEZ, checkValidity=2L)

mathingVector <- unlist(sapply(1:length(EEZ), function(x) which(EEZDF$EEZ ==  slot(slot(EEZ, "polygons")[[x]],"ID") )[1]))
length(mathingVector) == length(EEZ)

EEZDF <- EEZDF[mathingVector,]
rownames(EEZDF) <- unlist(sapply(1:length(EEZ), function(x) slot(slot(EEZ, "polygons")[[x]],"ID") ))
EEZ <- SpatialPolygonsDataFrame(EEZ, EEZDF )

EEZ <- st_as_sf(EEZ)
EEZ <- st_simplify(EEZ, preserveTopology=TRUE, dTolerance = 0.01)

EEZ <- EEZ[,c("fid","EEZ","IHO_SEA","SOVEREIGN1","ISO_SOV1","geometry")]
names(EEZ) <- c("fid","EEZ","seaBasin","country","countryCode","geometry")

st_write(st_as_sf(EEZ), paste0(dumpFolder,"/Spatial/EEZ/EEZGlobal.shp"), delete_layer = TRUE) # overwrites
save(EEZ,file=paste0(dumpFolder,"/Spatial/EEZ/EEZGlobal.RData"))

# ---------------------

load(file=paste0(dataFolder,"/Spatial/EEZ/EEZGlobal.RData"))

MPA1 <- shapefile(paste0(dataFolder,"/Spatial/WDPA/Shp0/WDPA_Mar2022_Public_shp-polygons.shp"))
MPA1 <- MPA1[MPA1$MARINE == "2",]

MPA2 <- shapefile(paste0(dataFolder,"/Spatial/WDPA/Shp1/WDPA_Mar2022_Public_shp-polygons.shp"))
MPA2 <- MPA2[MPA2$MARINE == "2",]

MPA3 <- shapefile(paste0(dataFolder,"/Spatial/WDPA/Shp2/WDPA_Mar2022_Public_shp-polygons.shp"))
MPA3 <- MPA3[MPA3$MARINE == "2",]
MPA <- bind(MPA1, MPA2, MPA3)

MPADF <- as.data.frame(MPA)
MPADF$country <- sapply(MPADF$ISO3, function(x) { EEZ[which(EEZ$countryCode == x)[1] , ]$country })

MPA <- gUnaryUnion(MPA, id = MPA@data$NAME, checkValidity=2L)

mathingVector <- unlist(sapply(1:length(MPA), function(x) which(MPADF$NAME ==  slot(slot(MPA, "polygons")[[x]],"ID") )[1]))
length(mathingVector) == length(MPA)
MPADF <- MPADF[mathingVector,]
rownames(MPADF) <- unlist(sapply(1:length(MPA), function(x) slot(slot(MPA, "polygons")[[x]],"ID") ))
MPA <- SpatialPolygonsDataFrame(MPA, MPADF )
MPA$EEZ <- NA

MPA <- MPA[,c("NAME","EEZ","country","WDPAID","DESIG_ENG","IUCN_CAT","NO_TAKE","NO_TK_AREA")]
names(MPA) <- c("Name","EEZ","country","WDPAId","Designation","Category","noTake","noTakeArea")

sf::sf_use_s2(FALSE)
MPA.sf <- gBuffer(MPA, byid = T, width = 0.000000000001)
MPA.sf <- st_as_sf(MPA.sf)
MPA.sf <- st_simplify(MPA.sf, preserveTopology=TRUE, dTolerance = 0.001)

assignEEZName <- st_intersects(MPA.sf,EEZ)
assignEEZName <- sapply(1:length(assignEEZName),function(x) { ifelse( length(assignEEZName[x][[1]]) == 1 , assignEEZName[x][[1]] , NA ) } )

for(i in which(is.na(assignEEZName)) ) {
  
  cat(i ,"out of",nrow(MPA.sf),"\n")
  EEZName <- which( st_intersects(EEZ,MPA.sf[i,], sparse = FALSE)[,1] )
  
  if( length(EEZName) == 0) {
      bufferPolygon <- 0
      while( length(EEZName) == 0 ) { 
        bufferPolygon <- bufferPolygon + 0.1
        EEZName <- which( st_intersects(EEZ,st_buffer(MPA.sf[i,],dist=bufferPolygon), sparse = FALSE)[,1] )
  }  }
  
  if( length(EEZName) > 1 & ! is.na(MPA.sf[i,]$country) ) { EEZName <- which.min(stringdist(EEZName,MPA.sf[i,]$country)) }
  if( length(EEZName) > 1 & is.na(MPA.sf[i,]$country) ) { EEZName <- which.min(stringdist(EEZName,MPA.sf[i,]$Name)) }
  if( length(EEZName) != 1) { stop(i) }
  
  assignEEZName[i] <- EEZName

}

MPA.sf$EEZ <- EEZ[assignEEZName,]$EEZ

plot(MPA.sf[630,1])
MPA.sf[630,]

st_write(MPA.sf, paste0(dumpFolder,"/Spatial/WDPA/WDPAGlobal.shp"), delete_layer = TRUE)
save(MPA.sf,file=paste0(dumpFolder,"/Spatial/WDPA/WDPAGlobal.RData"))

# ----------------------------------------------------------------------
# ----------------------------------------------------------------------

load(file=paste0(dataFolder,"/Spatial/globalLandmass/landmassRegion.RData"))
load(file=paste0(dataFolder,"/Spatial/EEZ/EEZ.RData"))
load(file=paste0(dataFolder,"/Spatial/WDPA/WDPA.RData"))

# Make maps

mainMap <- ggplot() + 
  geom_polygon(data = landmass, aes(long,lat,group=group), colour = "#E6E6E6" , fill="#E6E6E6", size=0.1 ) +
  xlab("Longitude") + 
  ylab("Latitude") + 
  coord_map() + 
  theme_minimal() +
  theme(
    text = element_text(family = "Helvetica", color = "#22211d"),
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks = element_blank(),
    #axis.title.x = element_blank(),
    #axis.title.y = element_blank(),
    panel.grid.major = element_line(color = "black", size = 0.1),
    panel.grid.minor = element_blank(),
    plot.background = element_rect(fill = "#FFFFFF", color = NA), # F3F3F3
    panel.background = element_rect(fill = "#FFFFFF", color = NA), # F3F3F3
    legend.background = element_rect(fill = "#FFFFFF", color = NA),
    panel.border = element_blank() )

fig1 <- mainMap +
  geom_polygon(data=regionEZZ, aes(x = long, y = lat, group = group,fill = "#8CD1DF"), size = 0.3) +
  labs(fill='Economic Exclusive Zones boundaries') + theme(legend.position="none")

fig2 <- mainMap +
  geom_polygon(data=regionMPA, aes(x = long, y = lat, group = group,fill = "#8CD1DF"), size = 0.1) +
  labs(fill='Marine Protected Areas') + theme(legend.position="none")

dumpFolder <- "/Volumes/Jellyfish/Dropbox/Manuscripts/Global connectivity corridors and refugia of climatic analogs for marine biodiversity/Results/"

pdf(paste0(dumpFolder,"/Landmass.pdf"), width = 10 )
print(mainMap)
dev.off()

pdf(paste0(dumpFolder,"/EEZ.pdf"), width = 10 )
print(fig1)
dev.off()

pdf(paste0(dumpFolder,"/MPA.pdf"), width = 10 )
print(fig2)
dev.off()
