# ----------------------------------------------------------------------
# ----------------------------------------------------------------------
#
#
# ----------------------------------------------------------------------

closeAllConnections()
rm(list=(ls()[ls()!="v"]))
gc(reset=TRUE)
source("/Volumes/Jellyfish/Dropbox/theMarineDataScientist/gitRepositories/climaticAnalogs/mainFunctions.R")
source("~/Dropbox/theMarineDataScientist/gitRepositories/climaticAnalogs/mainFunctions.R")

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
# Importance of predictors

# Boxplot per scenario

relativeImportancePredictors <- loadRData(paste0(resultsFolder,"/novelClimatePredictorContrib.RData"))
relativeImportancePredictors <- loadRData(paste0(resultsFolder,"/disappearClimatePredictorContrib.RData"))

data <- data.frame()
dataMeans <- numeric(0)
for( pred in predictors) { 
  data <- rbind(data,data.frame(predictor=pred,value=relativeImportancePredictors[,pred]))
  dataMeans <- c(dataMeans,mean(relativeImportancePredictors[,pred], na.rm=T))
}
data <- data[!is.na(data$value),]

dataMeansDF <- data.frame(predictor=predictors,contribution=dataMeans)
dataMeans <- dataMeans[sort(predictors, index.return=T)$ix]

plot1 <- ggplot(data=data) + geom_boxplot(aes(x=predictor, y=value, fill=predictor), outlier.size = 0.25, outlier.color = "gray",notch = FALSE, fatten = NULL) +
  theme( legend.position="none", plot.title = element_text(size=11)) +
  ggtitle("Relative contribution of variables") +
  xlab("")

for(i in 1:length(predictors)) { 
  plot1 <- plot1 + geom_line(data = data.frame(x=c(-0.375+i,0.375+i), y=c(dataMeans[i],dataMeans[i])), aes(x=x, y=y), linetype = "dotted") 
}

pdf(file=paste0(resultsFolder,"/Figures/",scenario,"relativeImportanceVariables.pdf"), width=14, height=5)
plot1
dev.off()

# Donuts plot per scenario

data <- dataMeansDF
data$predictor <- predictorsReduced 

data$fraction <- data$contribution / sum(data$contribution)
data$ymax <- cumsum(data$fraction)
data$ymin <- c(0, head(data$ymax, n=-1))
data$labelPosition <- (data$ymax + data$ymin) / 2
data$label <- paste0(data$predictor, " (", round(data$contribution, digits = 1),"%)")

plot1 <- ggplot(data, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=predictor)) +
  geom_rect() + geom_text( x=1.5, aes(y=labelPosition, label=label), size=3) + 
  scale_fill_brewer(palette=3) + scale_color_brewer(palette=3) +
  coord_polar(theta="y") + xlim(c(-1, 4)) + theme_void() + theme(legend.position = "none")

pdf(file=paste0(resultsFolder,"/Figures/",scenario,"relativeImportanceVariablesDonut.pdf"), width=14, height=5)
plot1
dev.off()

# ---------------
# Mapping importance of predictors

relativeImportancePredictors <- loadRData(paste0(resultsFolder,"/novelClimatePredictorContrib.RData"))
relativeImportancePredictors <- loadRData(paste0(resultsFolder,"/disappearClimatePredictorContrib.RData"))

for( pred in predictors) {
  
  dataRaster <- rasterFromXYZ( data.frame(dataStructure[match(relativeImportancePredictors$cell,dataStructure$cell),c("x","y")],relativeImportancePredictors[,pred]) )
  writeRaster(dataRaster,file=paste0(resultsFolder,"/Raster/relativeImportance",pred,".tif"), format="GTiff", overwrite=TRUE)
  
  resultsName <- pred
  dataRaster <- as.data.frame(dataRaster,xy=TRUE,na.rm=T)
  colnames(dataRaster) <- c("Lon","Lat","Val")
  dataRaster <- dataRaster[dataRaster$Val > 0,]
  myColors <- c("#6FBBE8","#A1ECD8","#F6F9AB","#FCB46D","#B21414")
  
  plot1 <- mainMap +
    geom_tile(data = dataRaster, aes(x=Lon,y=Lat,fill=Val)) +
    scale_fill_gradientn(colours=myColors, na.value='transparent')
  
  pdf(file=paste0(resultsFolder,"/Figures/",scenario,"relativeImportance",pred,".pdf"), width=14, height=5)
  print(plot1)
  dev.off()
  
}

# ----------------------------------------------------------------------
# ----------------------------------------------------------------------
# Histogram of area per sigma class [classes of 0.1]

analogResultsSSP119 <- loadRData(paste0(gsub("ssp585","ssp119",resultsFolder),"/novelClimates",ROIName,".RData"))
analogResultsSSP585 <- loadRData(paste0(gsub("ssp119","ssp585",resultsFolder),"/novelClimates",ROIName,".RData"))

data <- rbind(data.frame(scenario="SSP119",novelClimateSigma=analogResultsSSP119$novelClimateSigma),
              data.frame(scenario="SSP585",novelClimateSigma=analogResultsSSP585$novelClimateSigma))

plot1 <- ggplot(data,aes(x=novelClimateSigma)) + 
  geom_histogram(data=subset(data,scenario == 'SSP119'),fill = "blue", alpha = 0.5,aes(y=..count../sum(..count..))) +
  geom_histogram(data=subset(data,scenario == 'SSP585'),fill = "red", alpha = 0.5,aes(y=..count../sum(..count..))) +
  theme_minimal() +
  xlab("Climatic dissimilarity") + 
  ylab("Proportion of the ocean") +
  geom_vline(xintercept=c(2), linetype="dashed", size=0.2) +
  annotate(geom="text", x=2.1, y=0.35, label="Moderate", color="black", hjust = 0) +
  geom_vline(xintercept=c(4), linetype="dashed", size=0.2) +
  annotate(geom="text", x=4.1, y=0.35, label="Extreme", color="black", hjust = 0)

pdf(file=paste0(resultsFolder,"/Figures/areaPerSigmaNovel.pdf"), width=14, height=5)
plot1
dev.off()

# --------

analogResultsSSP119 <- loadRData(paste0(gsub("ssp585","ssp119",resultsFolder),"/novelClimates",ROIName,".RData"))
analogResultsSSP585 <- loadRData(paste0(gsub("ssp119","ssp585",resultsFolder),"/novelClimates",ROIName,".RData"))

data <- rbind(data.frame(scenario="SSP119",disappearClimateSigma=analogResultsSSP119$disappearClimateSigma),
              data.frame(scenario="SSP585",disappearClimateSigma=analogResultsSSP585$disappearClimateSigma))

plot1 <- ggplot(data,aes(x=disappearClimateSigma)) + 
  geom_histogram(data=subset(data,scenario == 'SSP119'),fill = "blue", alpha = 0.5,aes(y=..count../sum(..count..))) +
  geom_histogram(data=subset(data,scenario == 'SSP585'),fill = "red", alpha = 0.5,aes(y=..count../sum(..count..))) +
  theme_minimal() +
  xlab("Disappearing climate") + 
  ylab("Proportion of the ocean") +
  geom_vline(xintercept=c(2), linetype="dashed", size=0.2) +
  annotate(geom="text", x=2.1, y=0.35, label="Moderate", color="black", hjust = 0) +
  geom_vline(xintercept=c(4), linetype="dashed", size=0.2) +
  annotate(geom="text", x=4.1, y=0.35, label="Extreme", color="black", hjust = 0)

pdf(file=paste0(resultsFolder,"/Figures/areaPerSigmaDisappearing.pdf"), width=14, height=5)
plot1
dev.off()

sum(analogResultsSSP119$novelClimateSigma < 2, na.rm=T) / length(analogResultsSSP119$novelClimateSigma)
sum(analogResultsSSP119$novelClimateSigma > 2 & analogResultsSSP119$novelClimateSigma < 4, na.rm=T) / length(analogResultsSSP119$novelClimateSigma)
sum(analogResultsSSP119$novelClimateSigma > 4, na.rm=T) / length(analogResultsSSP119$novelClimateSigma)

sum(analogResultsSSP585$novelClimateSigma < 2, na.rm=T) / length(analogResultsSSP585$novelClimateSigma)
sum(analogResultsSSP585$novelClimateSigma > 2 & analogResultsSSP585$novelClimateSigma < 4, na.rm=T) / length(analogResultsSSP585$novelClimateSigma)
sum(analogResultsSSP585$novelClimateSigma > 4, na.rm=T) / length(analogResultsSSP585$novelClimateSigma)

sum(analogResultsSSP119$disappearClimateSigma < 2, na.rm=T) / length(analogResultsSSP119$disappearClimateSigma)
sum(analogResultsSSP119$disappearClimateSigma > 2 & analogResultsSSP119$disappearClimateSigma < 4, na.rm=T) / length(analogResultsSSP119$disappearClimateSigma)
sum(analogResultsSSP119$disappearClimateSigma > 4, na.rm=T) / length(analogResultsSSP119$disappearClimateSigma)

sum(analogResultsSSP585$disappearClimateSigma < 2, na.rm=T) / length(analogResultsSSP585$disappearClimateSigma)
sum(analogResultsSSP585$disappearClimateSigma > 2 & analogResultsSSP585$disappearClimateSigma < 4, na.rm=T) / length(analogResultsSSP585$disappearClimateSigma)
sum(analogResultsSSP585$disappearClimateSigma > 4, na.rm=T) / length(analogResultsSSP585$disappearClimateSigma)

# ----------------------------------------------------------------------
# ----------------------------------------------------------------------
# Distance to climate analogs per scenario [remove end-on analogs]

analogResultsSSP119 <- loadRData(paste0(gsub("ssp585","ssp119",resultsFolder),"/novelClimates",ROIName,".RData"))
analogResultsSSP585 <- loadRData(paste0(gsub("ssp119","ssp585",resultsFolder),"/novelClimates",ROIName,".RData"))

data <- rbind(data.frame(scenario="SSP119",analogDistance=analogResultsSSP119$analogDistance),
              data.frame(scenario="SSP585",analogDistance=analogResultsSSP585$analogDistance))

plot1 <- ggplot(data,aes(x=analogDistance)) + 
  geom_histogram(data=subset(data,scenario == 'SSP119'),fill = "blue", alpha = 0.5,aes(y=..count../sum(..count..))) +
  geom_histogram(data=subset(data,scenario == 'SSP585'),fill = "red", alpha = 0.5,aes(y=..count../sum(..count..))) +
  theme_minimal() +
  xlab("Distance to analogue (km)") + 
  ylab("Proportion of the ocean")

pdf(file=paste0(resultsFolder,"/Figures/analogDistance.pdf"), width=14, height=5)
plot1
dev.off()

# ----------------------------------------------------------------------
# ----------------------------------------------------------------------
# Raster data per latitudinal bin

analogResultsSSP119 <- loadRData(paste0(gsub("ssp585","ssp119",resultsFolder),"/novelClimates",ROIName,".RData"))
analogResultsSSP585 <- loadRData(paste0(gsub("ssp119","ssp585",resultsFolder),"/novelClimates",ROIName,".RData"))

analogResultsSSP119$area <- raster::area(rasterFromXYZ(analogResultsSSP119[,c("x","y","novelClimateSigma")]))[analogResultsSSP119$cell]
analogResultsSSP585$area <- raster::area(rasterFromXYZ(analogResultsSSP585[,c("x","y","novelClimateSigma")]))[analogResultsSSP585$cell]

colnames(analogResultsSSP119)
resultsName <- "analogDistance" # refugia novelClimate novelClimateArea highNovelClimateArea disappearClimate analogDistance

perLatitude <- data.frame()
bins <- seq(-90,90, by=1)
bins <- data.frame(from=bins[-length(bins)],to=bins[-1])

for(i in 1:nrow(bins)) {
  
  if(resultsName == "refugia") {
    
    resultsNamePlot <- "Total area of refugia (km2)"
    resultSSP119 <- analogResultsSSP119[which(analogResultsSSP119$y >= bins[i,1] & analogResultsSSP119$y < bins[i,2] ),"refugia"]
    resultSSP585 <- analogResultsSSP585[which(analogResultsSSP585$y >= bins[i,1] & analogResultsSSP585$y < bins[i,2] ),"refugia"]
    
    resultSSP119Area <- analogResultsSSP119[which(analogResultsSSP119$y >= bins[i,1] & analogResultsSSP119$y < bins[i,2] ),"area"]
    resultSSP585Area <- analogResultsSSP585[which(analogResultsSSP585$y >= bins[i,1] & analogResultsSSP585$y < bins[i,2] ),"area"]
    
    perLatitude <- rbind(perLatitude,
                         data.frame(bin=bins[i,1],
                                    meanSSP119=sum(resultSSP119 * resultSSP119Area),
                                    meanSSP585=sum(resultSSP585 * resultSSP119Area)))
  }

  if(resultsName == "analogDistance") {
    
    resultsNamePlot <- "Mean distance to analogue (km)"
    resultSSP119 <- analogResultsSSP119[which(analogResultsSSP119$y >= bins[i,1] & analogResultsSSP119$y < bins[i,2] ),"analogDistance"]
    resultSSP585 <- analogResultsSSP585[which(analogResultsSSP585$y >= bins[i,1] & analogResultsSSP585$y < bins[i,2] ),"analogDistance"]
    
    perLatitude <- rbind(perLatitude,
                         data.frame(bin=bins[i,1],
                                    meanSSP119=mean(resultSSP119 , na.rm=T),
                                    meanSSP585=mean(resultSSP585 , na.rm=T)))
  }
  
  if(resultsName == "disappearClimate") {
    
    resultsNamePlot <- "Mean disappear climate (sigma)"
    resultSSP119 <- analogResultsSSP119[which(analogResultsSSP119$y >= bins[i,1] & analogResultsSSP119$y < bins[i,2] ),"disappearClimateSigma"]
    resultSSP585 <- analogResultsSSP585[which(analogResultsSSP585$y >= bins[i,1] & analogResultsSSP585$y < bins[i,2] ),"disappearClimateSigma"]
    
    perLatitude <- rbind(perLatitude,
                         data.frame(bin=bins[i,1],
                                    meanSSP119=mean(resultSSP119 , na.rm=T),
                                    meanSSP585=mean(resultSSP585 , na.rm=T)))
  }
  
  if(resultsName == "novelClimate") {
    
    resultsNamePlot <- "Mean novel climate (sigma)"
    resultSSP119 <- analogResultsSSP119[which(analogResultsSSP119$y >= bins[i,1] & analogResultsSSP119$y < bins[i,2] ),"novelClimateSigma"]
    resultSSP585 <- analogResultsSSP585[which(analogResultsSSP585$y >= bins[i,1] & analogResultsSSP585$y < bins[i,2] ),"novelClimateSigma"]

    perLatitude <- rbind(perLatitude,
                         data.frame(bin=bins[i,1],
                                    meanSSP119=mean(resultSSP119 , na.rm=T),
                                    meanSSP585=mean(resultSSP585 , na.rm=T)))
  }
  
  if(resultsName == "novelClimateArea") {
    
    resultsNamePlot <- "Total area of moderate new climates (km2)"
    resultSSP119 <- analogResultsSSP119[which(analogResultsSSP119$y >= bins[i,1] & analogResultsSSP119$y < bins[i,2] ),"novelClimateSigma"]
    resultSSP585 <- analogResultsSSP585[which(analogResultsSSP585$y >= bins[i,1] & analogResultsSSP585$y < bins[i,2] ),"novelClimateSigma"]
    
    resultSSP119 <- resultSSP119 >= 2 & resultSSP119 < 4
    resultSSP585 <- resultSSP585 >= 2 & resultSSP585 < 4
    
    resultSSP119Area <- analogResultsSSP119[which(analogResultsSSP119$y >= bins[i,1] & analogResultsSSP119$y < bins[i,2] ),"area"]
    resultSSP585Area <- analogResultsSSP585[which(analogResultsSSP585$y >= bins[i,1] & analogResultsSSP585$y < bins[i,2] ),"area"]
    
    perLatitude <- rbind(perLatitude,
                         data.frame(bin=bins[i,1],
                                    meanSSP119=sum(resultSSP119 * resultSSP119Area),
                                    meanSSP585=sum(resultSSP585 * resultSSP119Area)))
    
  }
  
  if(resultsName == "highNovelClimateArea") {
    
    resultsNamePlot <- "Total area of extreme novel climate (km2)"
    resultSSP119 <- analogResultsSSP119[which(analogResultsSSP119$y >= bins[i,1] & analogResultsSSP119$y < bins[i,2] ),"novelClimateSigma"]
    resultSSP585 <- analogResultsSSP585[which(analogResultsSSP585$y >= bins[i,1] & analogResultsSSP585$y < bins[i,2] ),"novelClimateSigma"]
    
    resultSSP119 <- resultSSP119 >= 4
    resultSSP585 <- resultSSP585 >= 4
    
    resultSSP119Area <- analogResultsSSP119[which(analogResultsSSP119$y >= bins[i,1] & analogResultsSSP119$y < bins[i,2] ),"area"]
    resultSSP585Area <- analogResultsSSP585[which(analogResultsSSP585$y >= bins[i,1] & analogResultsSSP585$y < bins[i,2] ),"area"]
    
    perLatitude <- rbind(perLatitude,
                         data.frame(bin=bins[i,1],
                                    meanSSP119=sum(resultSSP119 * resultSSP119Area),
                                    meanSSP585=sum(resultSSP585 * resultSSP119Area)))
    
  }
  
}

fig3 <- ggplot() +
  geom_line( data=perLatitude,aes(x=bin, y=meanSSP119) , color="blue", size=0.4) +
  geom_line( data=perLatitude,aes(x=bin, y=meanSSP585) , color="red", size=0.4) +
  theme_minimal() +
  xlab("Latitude") + 
  ylab(resultsNamePlot) + geom_hline(yintercept=0, linetype="dotted", color = "black", size=0.7) + coord_flip()
fig3

pdf(paste0(resultsFolder,"/Figures/",resultsName," Latitude.pdf"), width = 10 )
fig3
dev.off()

# ----------------------------------------------------------------------
# ----------------------------------------------------------------------
# Cumulative % of ocean with new climate [sigma > 4] over time

# Needs time construction !

analogResultsSSP119 <- loadRData(paste0(gsub("ssp585","ssp119",resultsFolder),"/climaticAnalogsSigma.RData"))
analogResultsSSP585 <- loadRData(paste0(gsub("ssp119","ssp585",resultsFolder),"/climaticAnalogsSigma.RData"))

cumNewClimateAreaSSP119 <- numeric(ncol(analogResultsSSP119))
cumNewClimateAreaSSP585 <- numeric(ncol(climaticAnalogsSigmaSSP585))

for(t in 1:ncol(analogResultsSSP119)) {
  cumNewClimateAreaSSP119[t] <- sum(analogResultsSSP119[,t] > 4) / nrow(analogResultsSSP119)
  cumNewClimateAreaSSP585[t] <- sum(climaticAnalogsSigmaSSP585[,t] > 4) / nrow(climaticAnalogsSigmaSSP585)
}

dataToPlot <- data.frame(Time=2021:2100,NewclimateSSP119=cumNewClimateAreaSSP119,NewclimateSSP585=cumNewClimateAreaSSP585)

fig1 <- ggplot() +
  geom_line( data=dataToPlot,aes(x=Time, y=NewclimateSSP119) , color="blue", size=0.2) +
  geom_line( data=dataToPlot,aes(x=Time, y=NewclimateSSP585) , color="red", size=0.2) +
  theme_minimal() +
  xlab("Time") + 
  ylab("Novel climates (proportion of the ocean)")

pdf(paste0(resultsFolder,"/Figures/Proportion of the ocean with very high novel climate.pdf"), width = 10 )
fig1
dev.off()

analogResultsSSP119 <- loadRData(paste0(gsub("ssp585","ssp119",resultsFolder),"/novelClimates",ROIName,".RData"))
analogResultsSSP585 <- loadRData(paste0(gsub("ssp119","ssp585",resultsFolder),"/novelClimates",ROIName,".RData"))

dataToPlot <- data.frame(Time=2021:2100,
                         DistToAnalogSSP119=apply(distToAnalogSSP119,2,mean),
                         DistToAnalogSSP585=apply(distToAnalogSSP585,2,mean) )

fig2 <- ggplot() +
  geom_line( data=dataToPlot,aes(x=Time, y=DistToAnalogSSP119) , color="blue", size=0.2) +
  geom_line( data=dataToPlot,aes(x=Time, y=DistToAnalogSSP585) , color="red", size=0.2) +
  theme_minimal() +
  xlab("Time") + 
  ylab("Average distance to best climatic analog (km)")

pdf(paste0(resultsFolder,"/Figures/Distance to analog.pdf"), width = 10 )
fig2
dev.off()
