### Novelty_Oceans Functions
### KE Lotterhos
### Oct 2018 - mod. June 2019
### Northeastern University
### Mod. Áki Nov. 2018


#Create function that removes previous user installed packages to avoid masking
#clean_pkgs<-function(){
#  lapply(paste('package:',names(sessionInfo()$otherPkgs),sep=""),detach,character.only=TRUE)
#}
#clean_pkgs() #Remove all non-essential previously called packages

packages_needed <- c("raster", "FNN", "RColorBrewer", "colorRamps", "adehabitatLT",
                     "data.table", "tidyverse", "fields", "ggplot2", "hexbin",
                     "rgdal", "tmap", "gstat", "sp", "maptools", "sf", "fasterize",
                     "fansi", "raster", "tmap", "gstat", "ContourFunctions","ash",
                     "gridExtra", "shape", "reshape2"
                     )

for (i in 1:length(packages_needed)){
  if(!(packages_needed[i] %in% installed.packages())){install.packages(packages_needed[i])}
}

#  require(ggplot2)}

for (i in 1:length(packages_needed)){
  library( packages_needed[i], character.only = TRUE)
}

#--------------------------------
#### 40-year climate normals ####
#--------------------------------
# summer  
# (if Lat > 0, months 6,7,8)
# (if Lat < 0, months 12,1,2)
# winter
# (if Lat < 0, months 6,7,8)
# (if Lat > 0, months 12,1,2)
# means: SST_summer, SST_winter, Arag_summer, Arag_winter, 
# Calc_summer, Calc_winter, pH_summer, pH_winter

calculate_normals <- function(dat1){
  # input the data frame for the span of years you want the normals for
  
  month1 <- c(6,7,8)
  month2 <- c(12,1,2)
  
  ## Summer calculations
  x_sum <- dat1 %>% filter((Lat <= 0 & Month %in% month2)|
                             (Lat >= 0 & Month %in% month1))
  SST_sum <- tapply(x_sum$SST,INDEX = x_sum$No,mean, rm.na=TRUE)
  #length(SST_sum)
  Arag_sum <- tapply(x_sum$Arag,INDEX = x_sum$No,mean, rm.na=TRUE)
  #length(Arag_sum)
  pH_sum <- tapply(x_sum$pH,INDEX = x_sum$No,mean, rm.na=TRUE)
  
  
  # check that the column names are identical:
  if(!(identical(names(SST_sum), names(Arag_sum)) |
     #  identical(names(Arag_sum), names(Calc_sum)) |
       identical(names(pH_sum), names(SST_sum)))
    ){break}
  
  # create summer data frame for each station
  smr <- data.frame(No=as.integer(names(SST_sum)), SST_sum,Arag_sum, pH_sum  )
  head(smr)
  
  ## Winter calculations    
  x_win <- dat1 %>% filter((Lat > 0 & Month %in% month2)|
                             (Lat <0 & Month %in% month1))
  SST_win <- tapply(x_win$SST,INDEX = x_win$No,mean, rm.na=TRUE)
  Arag_win <- tapply(x_win$Arag,INDEX = x_win$No,mean, rm.na=TRUE)
  pH_win <- tapply(x_win$pH,INDEX = x_win$No,mean, rm.na=TRUE)

  wnt <- data.frame(No=as.integer(names(SST_win)), SST_win,Arag_win, pH_win  )
  head(wnt)

  # merge summer and winter data frames
  normals <- full_join(smr, wnt, by="No")
  #head(normals)
  cond <- which(!complete.cases(normals))
  normals[cond,]
  return(normals[order(normals$No),])
}



calculate_normals_annual <- function(dat1){
  # input the data frame for the span of years you want the normals for
  
  ## Summer calculations
  SST <- tapply(dat1$SST,INDEX = dat1$No,mean, rm.na=TRUE)
  #length(SST_sum)
  Arag <- tapply(dat1$Arag,INDEX = dat1$No,mean, rm.na=TRUE)
  #length(Arag_sum)
  pH <- tapply(dat1$pH,INDEX = dat1$No,mean, rm.na=TRUE)
  
  
  # check that the column names are identical:
  if(!(identical(names(SST), names(Arag)) |
       #  identical(names(Arag_sum), names(Calc_sum)) |
       identical(names(pH), names(SST)))
  ){break}
  
  # create summer data frame for each station
  normals <- data.frame(No=as.integer(names(SST)), SST,Arag, pH  )
  return(normals[order(normals$No),])
}

#--------------------------------  
### Calculate Sigma Dissimilarity ####
#--------------------------------
#initiate the data frame to store the projected sigma dissimilarity of 
# best analogs for each grid cell. 

loop_sigma_D <- function(A, B, C, append="", makePlot=FALSE){
  # A is baseline, the base that a station from B is compared to
  # If A is past and B is future, then degree of novelty
  # If A is future and B is past, then degree of disappearance
  # B will be subset to a particular station for analysis
  # C is the data frame used to calculate the ICV, will also be subset to a particular station for analysis
  if(!identical(A$No, B$No)){print("Error"); break}
  
  NN.sigma <- data.frame(No=A$No, NN.sigma=NA, NN.station=NA, NN.Mdist=NA, numPCs=NA)
  
  for(j in 1:nrow(NN.sigma)){ 
    NN.sigma[j,2:5] <- calc_sigma_D(A, B, C, NN.sigma$No[j], append, makePlot)
    if(j%%100==0){print(c(j, NN.sigma$No[j], nrow(NN.sigma)))}
  }
  names(NN.sigma)[2:5] <- paste0(names(NN.sigma)[2:5],append)
  return(NN.sigma)
}

  
calc_sigma_D <- function(A, B, C, whichStation, append="", makePlot=FALSE){
  # A is baseline, the base that a station from B is compared to
    # If A is past and B is future, then degree of novelty
    # If A is future and B is past, then degree of disappearance
  # B will be subset to a particular station for analysis
  # C is the data frame used to calculate the ICV, will also be subset to a particular station for analysis
  if(!identical(A$No, B$No)){print("Error"); break}
  if(!identical(dim(A), dim(B))){print("Error A and B different dimensions")}
  if(ncol(A) != ncol(C)){print("Error A and C different number of columns")}
  
  C.id <- C$No
  proxy <- B$No
  length(proxy)
  proxy2 <- sort(unique(proxy))
  if(!identical(proxy, proxy2)){break}
  
  # Principal component truncation rule
  trunc.SDs <- 0.1 #truncation 

  j <- whichStation
    # 29736, 29789
      
    # run the novelty calculation once for each ICV proxy. 
    # Takes about 1.5 sec/iteration on a typical laptop. 
    
    ## Select data relevant to ICV proxy j
    Bj <- B[which(proxy==j),]   # set of data from station j
    # select locations for which ICV proxy j is the closest ICV proxy. 
    Cj <- C[which(C.id==j),]    # reference period ICV at ICV proxy j
    
    ## Step 1: express climate data as standardized anomalies of reference period 
    #  ICV at ICV proxy j. 
    Cj.sd <- apply(Cj,2,sd, na.rm=T)  #standard deviation of interannual variability in each climate variable, ignoring missing years
    #standard deviation of variability in each climate 
    # variable, ignoring missing years
    A.prime <- sweep(A[,-1],MARGIN=2,Cj.sd[-1],`/`) #standardize the reference ICV
    # a <- matrix(c(1,2,3,4,5,6), nrow=2)
    # sweep(a, MARGIN =2, STATS=c(2,3,4)) # subtracts STATs from each column
    # sweep(a, MARGIN =2, STATS=c(2,3,4), FUN=`/`) # divides each column by STATS
    Bj.prime <- sweep(Bj[,-1],MARGIN=2,Cj.sd[-1],`/`) #standardize the analog pool    
    Cj.prime <- sweep(Cj[,-1],MARGIN=2,Cj.sd[-1],`/`) #standardize the projected future conditions of grid cells represented by ICV proxy j
    
    colnames(Cj.prime) <- colnames(A.prime)
    ## Step 2: Extract the principal components (PCs) of the reference period ICV 
    # and project all data onto these PCs
    PCA <- prcomp(Cj.prime[!is.na(apply(Cj.prime,1,mean)),])   
    # Principal components analysis. The !is.na(apply(...)) term is there 
    # simply to select all years with complete observations in all variables. 
    PCA$rotation
    
    #plot(PCA$rotation[,1], PCA$rotation[,2], xlim=c(-0.6, 0.1))
    #text(PCA$rotation[1:4,1], PCA$rotation[1:4,2], rownames(PCA$rotation)[1:4])
    # SST right of PC space, 
    # Arag and Calc and pH in upper left of PC space
    
    #plot(PCA$rotation[,1], PCA$rotation[,3], xlim=c(-0.6, 0.1))
    #text(PCA$rotation[1:4,1], PCA$rotation[1:4,3], rownames(PCA$rotation)[1:4])
    # separates pH from Calc and Arag
    
    #plot(PCA$rotation[,1], PCA$rotation[,4], xlim=c(-0.6, 0.1))
    #text(PCA$rotation[1:4,1], PCA$rotation[1:4,4], rownames(PCA$rotation)[1:4])
    # separates Calc from Arag
    
    #plot(PCA$sdev)
    #round(PCA$sdev, 2)
    
    (PCs <- max(c(which(unlist(summary(PCA)[1])>trunc.SDs),1)))
    # find the number of PCs to retain using the PC truncation 
    # rule of eigenvector stdev > the truncation threshold
    X <- as.data.frame(predict(PCA,A.prime))   # X is new baseline
    # project the analog pool onto the PCs
    head(X)
    
    Yj <- as.data.frame(predict(PCA,Bj.prime)) 
    # project the queried location onto the PCs
    
    Zj <- as.data.frame(predict(PCA,Cj.prime)) 
    # project the reference ICV onto the PCs
    
    
    #plot(X[,1], X[,2], xlim=c(-60, 20)) # analog
    #points(Zj[,1], Zj[,2], pch=19, col=rgb(1,0,0, 0.5)) # reference ICV
    #points(Yj[,1], Yj[,2], pch=19, col=rgb(0,1,0)) # future conditions
    
    
    ## Step 3a: express PC scores as standardized anomalies of reference interannual variability 
    Zj.sd <- apply(Zj,2,sd, na.rm=T)     
    #standard deviation of 1951-1990 interannual variability in each principal component, ignoring missing years
    #Zj.sd
    X.prime <- sweep(X,MARGIN=2,Zj.sd,`/`) # standardized baseline
    #standardize the analog pool  
    #head(X.prime)
    Yj.prime <- sweep(Yj,MARGIN=2,Zj.sd,`/`) 
    Zj.prime <- sweep(Zj, MARGIN=2, Zj.sd,`/`)
    
    #plot(X.prime[,1], X.prime[,2], xlim=c(-60, 20)) # analog
    #points(Zj.prime[,1], Zj.prime[,2], pch=19, col=rgb(1,0,0, 0.5)) # reference ICV
    #points(Yj.prime[,1], Yj.prime[,2], pch=19, col=rgb(0,1,0)) # future conditions
    #standardize the projected conditions   
    #Yj.prime
    ## Step 3b: find the sigma dissimilarity of each projected condition with 
    # its best analog (Euclidean nearest neighbour) in the observed analog pool.
    #X.prime <- X.prime[complete.cases(X.prime),] # standardized baseline
    nnd <- get.knnx(data=X.prime[,1:PCs],
                    query=Yj.prime[,1:PCs],
                    k=1,algorithm="brute")
    NN.dist <- as.vector(nnd[[2]]) 
    # Euclidean nearest neighbour distance in the z-standardized PCs of 
    # interannual climatic variability, i.e. the Mahalanobian nearest neighbour. 
    NN.chi <- pchi(NN.dist,PCs, rel.tol=.Machine$double.eps^0.8) 
    # percentile of the nearest neighbour 
    # increase the precision with 'rel.tol' #(default is ^0.5) to help with rounding error issues
    # distance on the chi distribution with degrees of freedom 
    # equaling the dimensionality of the distance measurement (PCs)
    if(NN.chi>=(1-1e-16)){NN.chi=1-1e-16}
      # sometimes with rounding error the NN.chi is printed as 1 but so close to 1 that
      # small changes in the decimal place equal large changes in sigma.
      # Also, if NN.chi equals 1 exactly than NN.sigma is infinite
      # this slight transformation gives the largest possible value of NN.sigma
    NN.sigma <- qchi(NN.chi,1) 
    # values of the chi percentiles on a standard half-normal distribution (chi distribution with one degree of freedom)
    if(NN.dist>9){NN.sigma <- qchi(1-1e-16,1) }
      # there are still issues with rounding error in the tail of the chi squared distribution
      # see my notebook post in the manuscript draft about this
    NN.station <- A$No[nnd$nn.index]
    #
    # Evaluation and plotting code
    #NN.sigma
    #X.prime[nnd$nn.index, 1:PCs]
    #Yj.prime[,1:PCs]
    #rbind(Bj,A[nnd$nn.index,])
    
    
    if(whichStation %in% c(18906,18952, 29736, 29789, 8193, 29319)){makePlot=TRUE}
      # 18906 is equator, high novelty and high disappearance
      # 18952 is n. pacific, low novelty but high disappearance
    if(makePlot==TRUE){
      
    png(paste0(append, "_",whichStation,  "_sigma=", round(NN.sigma,2),"_Md=",round(NN.dist,2),"_numPCs=",PCs,"_blue-focal_red-NN.png"), 
        height=8, width=8, units="in", res=600)
    par(mar=c(4,4,1,1), mfrow=c(2,2))
    
    scale=20
    n = nrow(PCA$rotation)
    
    # plot first and second PC
    x <- c(X.prime[,1], Yj.prime[1,1], Zj.prime[,1])
    y <- c(X.prime[,2], Yj.prime[1,2], Zj.prime[,2])
    plot(x, y,
         xlab="PC1", ylab="PC2", bty="l", las=1, col=adjustcolor("grey", 0.5),
         xlim=c(min(c(x, PCA$rotation[1:n,1]*PCA$sdev[1]*scale)), 
                max(c(x, PCA$rotation[1:n,1]*PCA$sdev[1]*scale))),
         ylim=c(min(c(y,(PCA$rotation[1:n,2]*PCA$sdev[2]*scale))),
                max(c(y,(PCA$rotation[1:n,2]*PCA$sdev[2]*scale))))
         )
    points(Zj.prime[,1], Zj.prime[,2], col=adjustcolor("magenta", 0.3), cex=1)
    points(Yj.prime[1,1], Yj.prime[1,2], col="blue", pch=19, cex=1)
    points(X.prime[nnd$nn.index,1], X.prime[nnd$nn.index,2], col=adjustcolor("black",0.8), bg=adjustcolor("green",0.7), pch=23, cex=1)
    Arrows(x0=rep(0, 6), y0=rep(0, 6), 
           x1=(PCA$rotation[1:n,1]*PCA$sdev[1]*scale), 
           y1=(PCA$rotation[1:n,2]*PCA$sdev[2]*scale),
           code=2, col=adjustcolor("black", 0.2), arr.type="curved")
    text(x=(PCA$rotation[1:n,1]*PCA$sdev[1]*scale), 
         y=(PCA$rotation[1:n,2]*PCA$sdev[2]*scale), 
         labels=rownames(PCA$rotation)[1:6], col=adjustcolor("black",0.8), adj=0.5)
    xa <- par("usr")[1] + (par("usr")[2]-par("usr")[1])*0.1
    ya <- par("usr")[4] - (par("usr")[4]-par("usr")[3])*0.1
    text(xa, ya, "A",  cex=2) 
    legend(xa, ya*0.95, c("Focal station", "NN"), col=c("blue", "black"), pch=c(19,23), 
           pt.bg=c("black", "green"), bty="n")
    
    # plot first and third PC
    x <- c(X.prime[,1], Yj.prime[1,1], Zj.prime[,1])
    y <- c(X.prime[,2], Yj.prime[1,3], Zj.prime[,3])
    plot(x, y,
         xlab="PC1", ylab="PC3", bty="l", las=1, col=adjustcolor("grey", 0.5),
         xlim=c(min(c(x, PCA$rotation[1:n,1]*PCA$sdev[1]*scale)), 
                max(c(x, PCA$rotation[1:n,1]*PCA$sdev[1]*scale))),
         ylim=c(min(c(y,(PCA$rotation[1:n,3]*PCA$sdev[3]*scale))),
                max(c(y,(PCA$rotation[1:n,3]*PCA$sdev[3]*scale))))
    )
    points(Zj.prime[,1], Zj.prime[,3], col=adjustcolor("magenta", 0.3), cex=1)
    points(Yj.prime[1,1], Yj.prime[1,3], col="blue", pch=19, cex=1)
    points(X.prime[nnd$nn.index,1], X.prime[nnd$nn.index,3], 
           col=adjustcolor("black",0.8), bg=adjustcolor("green",0.7), pch=23, cex=1)
    Arrows(x0=rep(0, 6), y0=rep(0, 6), 
           x1=(PCA$rotation[1:n,1]*PCA$sdev[1]*scale), 
           y1=(PCA$rotation[1:n,3]*PCA$sdev[3]*scale),
           code=2, col=adjustcolor("black", 0.2), arr.type="curved")
    text(x=(PCA$rotation[1:n,1]*PCA$sdev[1]*scale), 
         y=(PCA$rotation[1:n,3]*PCA$sdev[3]*scale), 
         labels=rownames(PCA$rotation)[1:n], col=adjustcolor("black",0.8), adj=0.5)
    xa <- par("usr")[1] + (par("usr")[2]-par("usr")[1])*0.1
    ya <- par("usr")[4] - (par("usr")[4]-par("usr")[3])*0.1
    text(xa, ya, "B",  cex=2) 
    
    # plot(c(A$SST_sum, Bj$SST_sum),c(A$Arag_sum,Bj$Arag_sum), xlab="Summer SST", ylab="Summer Aragonite S.S.",
    #      bty="l", las=1, col=adjustcolor("grey", 0.5))
    # points(Cj$SST, Cj$Arag, col=adjustcolor("magenta", 0.3))
    # points(Bj$SST_sum, Bj$Arag_sum, col="blue", pch=19, cex=2)
    # points(A$SST_sum[nnd$nn.index], A$Arag_sum[nnd$nn.index], col="black", bg="red", pch=23, cex=2)
    # xa <- par("usr")[1] + (par("usr")[2]-par("usr")[1])*0.1
    # ya <- par("usr")[4] - (par("usr")[4]-par("usr")[3])*0.1
    # text(xa, ya, "D",  cex=2) 
    # 
    # plot(c(A$pH_sum,Bj$pH_sum),c(A$Arag_sum,Bj$Arag_sum), xlab="Summer pH", ylab="Summer Aragonite S.S.",
    #      bty="l", las=1, col=adjustcolor("grey", 0.5))
    # points(Cj$pH, Cj$Arag, col=adjustcolor("magenta", 0.3))
    # points(Bj$pH_sum, Bj$Arag_sum, col="blue", pch=19, cex=2)
    # points(A$pH_sum[nnd$nn.index], A$Arag_sum[nnd$nn.index], col="black", bg="red", pch=23, cex=2)
    # xa <- par("usr")[1] + (par("usr")[2]-par("usr")[1])*0.1
    # ya <- par("usr")[4] - (par("usr")[4]-par("usr")[3])*0.1
    # text(xa, ya, "D",  cex=2) 

    plot(c(A$SST, Bj$SST, Cj$SST),c(A$Arag,Bj$Arag, Cj$Arag), xlab="SST", ylab="Aragonite S.S.",
          bty="l", las=1, col=adjustcolor("grey", 0.5))
     points(Cj$SST, Cj$Arag, col=adjustcolor("magenta"))
     points(Bj$SST, Bj$Arag, col="blue", pch=19, cex=1)
     points(A$SST[nnd$nn.index], A$Arag[nnd$nn.index], 
            col=adjustcolor("black",0.8), bg=adjustcolor("green",0.7), pch=23, cex=1)
     xa <- par("usr")[1] + (par("usr")[2]-par("usr")[1])*0.1
     ya <- par("usr")[4] - (par("usr")[4]-par("usr")[3])*0.1
     text(xa, ya, "C",  cex=2)
    #
     plot(c(A$pH,Bj$pH, Cj$pH),c(A$Arag,Bj$Arag, Cj$Arag), xlab="pH", ylab="Aragonite S.S.",
          bty="l", las=1, col=adjustcolor("grey", 0.5))
     points(Cj$pH, Cj$Arag, col=adjustcolor("magenta"))
     points(Bj$pH, Bj$Arag, col="blue", pch=19, cex=1)
     points(A$pH[nnd$nn.index], A$Arag[nnd$nn.index], 
            col=adjustcolor("black",0.8), bg=adjustcolor("green",0.7), pch=23, cex=1)
     xa <- par("usr")[1] + (par("usr")[2]-par("usr")[1])*0.1
     ya <- par("usr")[4] - (par("usr")[4]-par("usr")[3])*0.1
     text(xa, ya, "D",  cex=2)
    dev.off()
    
    }
    
    #names(NN.station) <- paste0("NN.station",names(NN.station))
    
    return(data.frame(NN.sigma, NN.station, NN.Mdist=NN.dist, num_PCs=PCs))
}


#--------------------------------  
### function for plotting ####
#--------------------------------
world <- map_data("world2")
Plot_nonInt<-function(lat, long, var, refMap, legend_name){
  sampData <- data.frame(lat, long, var)
  
  ggplot()+
    geom_polygon(data = refMap, aes(x=long, y = lat, group=group))+
    stat_summary_2d(data=sampData, aes(x=long, y = lat, z= var), bins=80, alpha = 0.8)+
    #theme(legend.position=c(50,100))+
    #scale_fill_brewer(palette="Dark2") +
    scale_fill_gradientn(name=legend_name, 
                         colours=two.colors(40,start = "blue", 
                                            end="red", middle="orange")) +
    coord_fixed() 
}


## Function to create a random grid empty grid 
# n = number of cells, increase n to increase resolution
makeGrid<-function(EB2){
  repeat{
    gr <- as.data.frame(spsample(EB2, 'regular', n  = 50000))
    names(gr) <- c('X', 'Y')
    coordinates(gr) <- c('X', 'Y')
    gridded(gr) <- TRUE  
    fullgrid(gr) <- TRUE  # Create SpatialGrid object
    try(proj4string(gr) <- proj4string(EB2))
    if(is.na(proj4string(gr))==FALSE) return(gr)
  }
}

#Function to convert raster objects as tibbles, written by Sébastien Rochette
gplot_data <- function(x, maxpixels = 100000)  {
  x <- raster::sampleRegular(x, maxpixels, asRaster = TRUE)
  coords <- raster::xyFromCell(x, seq_len(raster::ncell(x)))
  ## Extract values
  dat <- utils::stack(as.data.frame(raster::getValues(x))) 
  names(dat) <- c('value', 'variable')
  
  dat <- dplyr::as.tbl(data.frame(coords, dat))
  
  if (!is.null(levels(x))) {
    dat <- dplyr::left_join(dat, levels(x)[[1]], 
                            by = c("value" = "ID"))
  }
  dat
}

### Function for interpolated plotting ####
Plot_interp<-function(data,column, B, title){ #Data is expected to be data.table format, column should be in the format 'data$column', title is a character string for the figure legend.
  FD12<-data[,c("long","lat")]
  FD12$sigma<-column
  FD12<-FD12[!is.na(FD12$sigma),]
  FD12[long>360,long:=long-360] #Correct format
  #FD12$long<-FD12$long-180 # Convert to WGS 1984 bounds
  FD12[long>180,long:=long-360]
  FD12_SP <- SpatialPoints(FD12) # Spatial points df
  proj4string(FD12_SP) <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84")
  
  boolFALSE<-F
  MakeRast<-function(data){
    while(boolFALSE==F){
      tryCatch({
        gr <- as.data.frame(spsample(data, 'regular', n  = 1000000))#n  = 10000000))
        names(gr) <- c('X', 'Y')
        coordinates(gr) <- c('X', 'Y')
        gridded(gr) <- TRUE 
        fullgrid(gr) <- TRUE
        proj4string(gr) <- proj4string(data)
        data.idw <- idw(sigma~1, data, newdata = gr, idp = 10)
        return(raster(data.idw))
        boolFalse<-T
      },
      error=function(e){
      },finally={})
    }
  }
  
  r<-MakeRast(FD12_SP)
  
  data('World', package = 'tmap')
  #get_projection(World)
  
  tm_shape(r) + 
    tm_raster(breaks = B, palette = 'plasma', # n = 10 may be better
               legend.show=FALSE) + 
    #tm_legend(legend.outside = TRUE) +
    tm_layout(main.title=title, main.title.size=1,
           #   outer.margins = c(0,0,0,0),
          #    inner.margins = c(0,0,0,0),
          #    between.margin = c(0,0,0,0)
          ) + 
    tm_shape(World, projection="longlat") +
    tm_fill()
}


### Function for interpolated plotting only legend ####
Plot_interp_legend <-function(data,column, B, title){ #Data is expected to be data.table format, column should be in the format 'data$column', title is a character string for the figure legend.
  FD12<-data[,c("long","lat")]
  FD12$sigma<-column
  FD12<-FD12[!is.na(FD12$sigma),]
  FD12[long>360,long:=long-360] #Correct format
  #FD12$long<-FD12$long-180 # Convert to WGS 1984 bounds
  FD12[long>180,long:=long-360]
  FD12_SP <- SpatialPoints(FD12) # Spatial points df
  proj4string(FD12_SP) <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84")
  
  boolFALSE<-F
  MakeRast<-function(data){
    while(boolFALSE==F){
      tryCatch({
        gr <- as.data.frame(spsample(data, 'regular', n  = 10000))#DO NOT INCREASE n FOR LEGEND))
        names(gr) <- c('X', 'Y')
        coordinates(gr) <- c('X', 'Y')
        gridded(gr) <- TRUE 
        fullgrid(gr) <- TRUE
        proj4string(gr) <- proj4string(data)
        data.idw <- idw(sigma~1, data, newdata = gr, idp = 10)
        return(raster(data.idw))
        boolFalse<-T
      },
      error=function(e){
      },finally={})
    }
  }
  
  r<-MakeRast(FD12_SP)
  
  data('World', package = 'tmap')
  #get_projection(World)
  
  tm_shape(r) + 
    tm_raster(breaks = B, palette = 'plasma', # n = 10 may be better
    title=title) + 
    tm_legend(legend.outside = FALSE, legend.position=c("left","center"),
              legend.title.size=1, legend.text.size=1) +
    tm_layout(title=title, legend.only=TRUE, main.title=title, main.title.size=1) + 
    tm_shape(World, projection="longlat") +
    tm_fill()
}
