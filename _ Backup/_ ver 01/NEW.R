closeAllConnections()
rm(list=(ls()[ls()!="v"]))
gc(reset=TRUE)
source("Dependencies/Functions.R")
#-------------------------------------------------------

# Gradien based climate velocity ----
## Load data ----
OST <- stack(list.files(path = "A:/CCMAR/Temperature/Yearly/SSP85/MEAN", pattern ='.tif', full.names = TRUE))

## Change raster resolution ----
OST <- aggregate(OST, fact = 2) 
res(OST)

## Temperature trend ----
vt <- tempTrend(OST,
                th = nlayers(OST)) # Number of years

## Spatial gradient ----
vg <- spatGrad(OST,
               th = 0.0001,   #Temperature threshold to be consider as change
               projected = FALSE)

## Temperature velocity ----
gv <- gVoCC(vt, 
            vg) 

vel <- gv[[1]] # Intensity
ang <- gv[[2]] # Direction

# Temperature Trajectories ----

## Creating lat long matrix ----
lonlat <- data.frame(xyFromCell(vel, 
                                1:ncell(vel)))
lonlat$vel <- raster::extract(vel, 
                              lonlat)
lonlat$ang <- raster::extract(ang, 
                              lonlat[,1:2])
lonlat$mn <- raster::extract(mn, 
                             lonlat[,1:2])
lonlat <- na.omit(lonlat)

## Mean temperature for the whole period ----
mn <- calc(OST, 
           mean, 
           na.rm = TRUE)
## Trajectories in parallel ----
cores <- detectCores()
ncores<- cores[1] ##ALL CORES
cuts <- cut(1:nrow(lonlat), 
            ncores)
cl <- makeCluster(ncores)
registerDoParallel(cl)
traj <- foreach(x = levels(cuts), 
                .combine = rbind, 
                .packages = c('raster','sp','rgeos','geosphere','rgdal','VoCC'), 
                .multicombine = TRUE) %dopar% {
                  voccTraj(lonlat[cuts == x,], 
                           vel, 
                           ang, 
                           mn, 
                           tyr = nlayers(OST), 
                           trajID = as.numeric(rownames(lonlat[cuts == x,])), 
                           correct = TRUE)}



stopCluster(cl)

## Trajectories clasification ----
traj_cl <- trajClas(traj, 
                    vel, 
                    ang, 
                    mn, 
                    trajSt = 4, 
                    tyr =  nlayers(OST), 
                    nmL = 20, 
                    smL = 100,
                    Nend = 45, 
                    Nst = 15, 
                    NFT = 70, 
                    DateLine = FALSE)

TrajClas <- traj_cl$TrajClas

## Plot trajectories ----
ggplot() +
  geom_tile(data = gplot_data(TrajClas), 
            aes(x = x, 
                y = y, 
                fill = factor(value))) +
  scale_fill_manual(values = c("#e1e1e1", "#85ce99", "#aa5f40",
                               "#ee393f", "#1e50a3", "#89cc9a",
                               "#c05ea6", "#01bbe7", "#f1eb79"), 
                    breaks = 1:9,
                    labels = c("Non-moving", "Slow-moving", "Internal sinks", 
                               "Boundary sinks", "Sources", "Relative sinks", 
                               "Corridors", "Divergence", "Convergence"),
                    name = "Trajectory class")

# Currents ----

## Load data ----
SC_x <- calc(stack(list.files(path = "A:/CCMAR/Currents/Yearly/SSP85/MAX", 
                              pattern ='X', 
                              full.names = TRUE)), 
             mean)
SC_y <- calc(stack(list.files(path = "A:/CCMAR/Currents/Yearly/SSP85/MAX", 
                              pattern ='Y', 
                              full.names = TRUE)), 
             mean)

SC <-  raster::stack(SC_x, 
                     SC_y)

## Change raster resolution ----
OST <- aggregate(SC, 
                 fact = 2) 

## Two currents directions ----
u <- SC[[1]]
v <- SC[[2]]

## Creating currents direction raster ----
Cur_Dir <- atan2(u, 
                 v)
Cur_Dir <- rad2deg(Cur_Dir)
Cur_Dir[Cur_Dir < 0] <- 360 + Cur_Dir[Cur_Dir < 0]

## Temperature velovity direction ----
Vel_Dir <- ang

## Both rasters at same resolution and projection ----
Cur_Dir <- projectRaster(Cur_Dir,
                         Vel_Dir,
                         method = 'bilinear') 
## Direction aggregment between temperature direction and current direction ----
Dir_agr <- Cur_Dir - Vel_Dir
Dir_agr <- cos(Dir_agr * pi / 180)

# Overlapping trajectories and aggregment direction ----

# Corridors <-  TrajClas == 7
# Slow_moving <- TrajClas == 2
# Non_moving <- TrajClas == 1

all <- TrajClas %in% c(1,7)
final_all <- all*Dir_agr
final_all_pos <- reclassify(final_all,
                            c( -Inf, 0 , 0))

## Plotting ----
ggplot() +
  geom_tile(data = gplot_data(final_all), 
            aes(x = x, 
                y = y, 
                fill = value))+
  scale_fill_gradient2(
    low = muted("#60e219"),
    mid = "gray",
    high = muted("#e31f03"),
    midpoint = 0)+
  labs(x = "Long", 
       y = "Lat",
       fill = "Directional\nAgreement")