## ---------------------------------------------------------------------
## ---------------------------------------------------------------------
##      
##  biodiversityDS [ biodiversitydatascience.com ]
##
## ---------------------------------------------------------------------
## ---------------------------------------------------------------------

closeAllConnections()
rm(list=(ls()[ls()!="v"]))
gc(reset=TRUE)

source("mainFunctions.R")
library('sf')

## -----------------

myTaxaDB.1 <- read.csv("../../Data/listIUCNSpecies.csv") 
myTaxaDB.2 <- read.csv("../../Data/listFishSpecies.csv")  
myTaxaDB.3 <- read.csv("../../Data/AquacultureSpeciesListAll.csv", sep=";")  
speciesToGetRecords <- unique(c(myTaxaDB.1$species,myTaxaDB.2$species,myTaxaDB.3$Species))

baseMapFile <- "MaskRes025.tif"
fileAquamaps <- "/Volumes/StingRay/Dropbox/Data/Biodiversity Data/Occurrence records/Biodiversity database Aquamaps.db"
filesIUCN <- list.files("/Volumes/StingRay/Dropbox/Data/Biodiversity Data/IUCN range maps", pattern = ".shp", recursive = TRUE, full.names = TRUE)
filesIUCN <- filesIUCN[!grepl("xml",filesIUCN)]

# -----------------------
# -----------------------

sqlite.driver <- dbDriver("SQLite")
db <- dbConnect(sqlite.driver, dbname = fileAquamaps)
taxaDB <- dbReadTable(db,"taxa")
taxaDB$scientificName <- taxa <- paste0(taxaDB$Genus," ",taxaDB$Species)

nativemaps <- dbReadTable(db,"nativemaps")
cells <- dbReadTable(db,"hcaf")
occurrence <- dbReadTable(db,"occ")

# -----------------------
# -----------------------

IUCNSpeciesDBAvailable <- data.frame()
for( f in filesIUCN  ) {
  filesIUCN.i <- read_sf(f)
  IUCNSpeciesDBAvailable <- rbind(IUCNSpeciesDBAvailable,data.frame(file=f,species=unique(filesIUCN.i$binomial)))
  rm(filesIUCN.i)
  gc(reset=TRUE)
}

# -----------------------
# -----------------------

summary <- data.frame()
rasterShape <- raster(baseMapFile)
rasterShape05 <- aggregate(rasterShape,2)
rgdal::setCPLConfigOption("GDAL_PAM_ENABLED", "FALSE")

for( sp in speciesToGetRecords[1:length(speciesToGetRecords)]) {
  
  cat("\014")  
  cat("\n")  
  cat("# --------------")
  cat("\n")  
  cat("Process",which(speciesToGetRecords == sp),"out of",length(speciesToGetRecords),"\n")
  cat("\n")  
  cat("# --------------")
  cat("\n")  

  # -----------------------------
  # IUCN 
  
  file.i <- IUCNSpeciesDBAvailable[IUCNSpeciesDBAvailable$species == sp,1]
  
  if(length(file.i) > 0) { 
    
    file.i <- read_sf(unique(file.i[1]))
    databaseIUCN.range <- file.i[file.i$binomial == sp,]
    databaseIUCN.range <- as_Spatial(databaseIUCN.range)
    iucnRangeMap <- raster::resample(rasterize(databaseIUCN.range,rasterShape05),rasterShape)
    iucnRangeMap[iucnRangeMap != 0] <- 1
    iucnRangeMap[iucnRangeMap != 1] <- 0
    iucnRangeMap[ is.na(iucnRangeMap)] <- 0
    iucnRangeMap <- mask(iucnRangeMap,rasterShape)
    # plot(iucnRangeMap)
    
    if( 1 %in% cellStats(iucnRangeMap,range) ) {
    
        rasterShape.i.pts <- iucnRangeMap
        rasterShape.i.pts <- xyFromCell(rasterShape.i.pts,Which(rasterShape.i.pts == 1 , cells=T))
        rasterShape.i.pts <- data.frame(decimalLongitude=rasterShape.i.pts[,1],decimalLatitude=rasterShape.i.pts[,2], sourceType="Native Maps")
        resultsDirectory.i <- paste0("../../Data/Rasters/",sp)
        if( ! dir.exists(resultsDirectory.i) ) { dir.create(resultsDirectory.i, recursive = T) }
        
        summary <- rbind(summary,data.frame(speciesName=sp, source="IUCN", occurrenceRecords=nrow(rasterShape.i.pts)))
        
        write.csv(rasterShape.i.pts,file=paste0(resultsDirectory.i,"/IUCNRecords.csv"), row.names = FALSE)
        raster::writeRaster(iucnRangeMap,filename=paste0(resultsDirectory.i,"/IUCNRecords.tif"),format="GTiff",overwrite=T)
    
    }
  }
  
  # -----------------------------
  # Aquamaps 

  databaseAquamaps.occ <- data.frame()
  databaseAquamaps.range <- data.frame()
  
  # Occurrences
  
  sp.id <- taxaDB[which(taxaDB$scientificName == sp),"SPECIESID"]
  
  if( length(sp.id) == 0) { next }
  
  sp.records <- occurrence[occurrence$SpeciesID == sp.id & occurrence$GoodCell == 1,"CsquareCode"]
  sp.records <- cells[cells$CsquareCode  %in% sp.records,c("CenterLong","CenterLat")]
  
  if(nrow(sp.records) > 0) {
    databaseAquamaps.occ <- data.frame( decimalLongitude=sp.records[,1],decimalLatitude=sp.records[,2] , speciesName=sp, sourceType="Occurrence Data")
  }
  
  # Native maps
  
  sp.records <- nativemaps[nativemaps$SpeciesID == sp.id & nativemaps$FAOAreaYN == 1 & nativemaps$probability >= 0.5 ,c("CsquareCode")]
  sp.records <- cells[cells$CsquareCode %in% sp.records,c("CenterLong","CenterLat")]
  
  if(nrow(sp.records) > 0) {
    databaseAquamaps.range <- data.frame( decimalLongitude=sp.records[,1],decimalLatitude=sp.records[,2] , speciesName=sp, sourceType="Native Maps")
  }

  databaseAquamaps <- rbind(databaseAquamaps.occ,databaseAquamaps.range)
  
  # ------------

  summary <- rbind(summary,data.frame(speciesName=sp, source="Aquamaps", occurrenceRecords=nrow(databaseAquamaps)))
  
  # ------------
  
  aquamapsRangeMap <- raster::resample(rasterize(databaseAquamaps[,1:2],rasterShape05, field=1),rasterShape)
  aquamapsRangeMap[aquamapsRangeMap != 1] <- 0
  aquamapsRangeMap[is.na(aquamapsRangeMap)] <- 0
  aquamapsRangeMap <- mask(aquamapsRangeMap,rasterShape)
  # plot(aquamapsRangeMap)
  
  resultsDirectory.i <- paste0("../../Data/Rasters/",sp)
  if( ! dir.exists(resultsDirectory.i) ) { dir.create(resultsDirectory.i, recursive = T) }
  
  write.csv(databaseAquamaps,file=paste0(resultsDirectory.i,"/aquamapsRecords.csv"), row.names = FALSE)
  raster::writeRaster(aquamapsRangeMap,filename=paste0(resultsDirectory.i,"/aquamapsRecords.tif"),format="GTiff",overwrite=T)
  
}

# ------------

write.csv(summary,file="../../Results/summaryRecords.csv", row.names = FALSE)

# --------------------------------------------------
# --------------------------------------------------