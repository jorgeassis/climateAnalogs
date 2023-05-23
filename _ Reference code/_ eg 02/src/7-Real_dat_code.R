## 2019-07-09
## Áki Jarl - Ocean Climate Novelty script to view variance at existing oceanographic meaurement opperations
##
source("src/0-Novelty_Oceans_Functions.R")

require(ggplot2)
require(bit64)
require(cowplot)

#Read in real data
hots <- fread("Data/real_data/HOT_bottle_data_mod.txt", sep = ",")
hots <- hots[hots$`date(mmddyy)` != -9] # Remove data which doesn't contain a date
hots[hots == -9] <- NA # Convert missing value columns from the numeric -9 to an NA 
bats <- fread("Data/real_data/BATS_bottle_data_mod.txt", sep = "\t")
bats[bats == -999] <- NA # Convert missing value columns from the numeric -9 to an NA 
bbh <- fread("Data/real_data/BBH.csv", sep = ",")
cml <- fread("Data/real_data/CML.csv", sep = ",")
cml[cml == -999] <- NA
cml[cml == "NaN"] <- NA

nh1 <- fread("Data/real_data/NH_70W_43N_2011.csv", sep = ",")
nh1 <- nh1[-c(1:4),c(4,25)]
nh1<-nh1[-1,]
colnames(nh1)<-c("Date","ph")
nh2 <- fread("Data/real_data/NH_70W_43N_2012.csv", sep = ",")
nh2 <- nh2[,c(4,25)]
colnames(nh2)<-c("Date","ph")
nh3 <- fread("Data/real_data/NH_70W_43N_2013.csv", sep = ",")
nh3 <- nh3[-c(1:4),c(4,25)]
nh3<-nh3[-1,]
colnames(nh3)<-c("Date","ph")
nh4 <- fread("Data/real_data/NH_70W_43N_2014.csv", sep = ",")
nh4 <- nh4[,c(4,25)]
colnames(nh4)<-c("Date","ph")
nh5 <- fread("Data/real_data/NH_70W_43N_2015.csv", sep = ",")
nh5 <- nh5[-c(1:4),c(4,25)]
nh5<-nh5[-1,]
colnames(nh5)<-c("Date","ph")
nh6 <- fread("Data/real_data/NH_70W_43N_2016.csv", sep = ",")
nh6 <- nh6[-c(1:4),c(4,25)]
nh6<-nh6[-1,]
colnames(nh6)<-c("Date","ph")

nh_all<-rbind(nh1,nh2)
nh_all<-rbind(nh_all,nh3)
nh_all<-rbind(nh_all,nh4)
nh_all<-rbind(nh_all,nh5)
nh_all<-rbind(nh_all,nh6)

nh_all[nh_all == -999] <- NA

#Read in projection data
dat_2100_4.5_new <- fread("Data/model/Temp_Arag_2070_2100_RCP45.txt", sep = ",")
colnames(dat_2100_4.5_new)<-c("No", "longitude", "latitude", "year", "month", "SST", "Aragonite","pH")
dat_2100_8.5_new <- fread("Data/model/Temp_Arag_2070_2100_RCP85.txt", sep = ",")
colnames(dat_2100_8.5_new)<-c("No", "longitude", "latitude", "year", "month", "SST", "Aragonite","pH")
dat_1800_new <- fread("Data/model/Temp_Arag_1800_2000.txt", sep = ",")
colnames(dat_1800_new)<-c("No", "longitude", "latitude", "year", "month", "SST", "Aragonite","pH")

#Reformat the data in the 'date' column to only include the year 
years<-NULL
for(i in as.numeric(substr(hots$`date(mmddyy)`,nchar(hots$`date(mmddyy)`)-1,nchar(hots$`date(mmddyy)`)))){
  if(i>80){
    years<-c(years,as.numeric(paste(19,i,sep="")))
  }
  else if(nchar(i)==2){
    years<-c(years,as.numeric(paste(20,i,sep="")))
  }
  else{
    years<-c(years,as.numeric(paste(200,i,sep="")))
  }
}

hots$Years<-years

mon<-NULL
for(i in 1:length(hots$`date(mmddyy)`)){
  if(nchar(hots$`date(mmddyy)`[i])==6){
    mon<-c(mon,substr(hots$`date(mmddyy)`[i],1,2))
  }else{
    mon<-c(mon,paste(0,substr(hots$`date(mmddyy)`[i],0,1),sep=""))
  }
}

hots$Month<-mon

bats$Years<-as.numeric(substr(bats$yyyymmdd,1,4))

bbh$Years<-as.numeric(substr(bbh$COLLECTION_DATE,nchar(bbh$COLLECTION_DATE)-3,nchar(bbh$COLLECTION_DATE)))

cml$Years<-as.numeric(substr(cml$time_UTC,1,4))
cml$YearsMonths<-substr(cml$time_UTC,1,7)

nh_all$Years<-as.numeric(paste("20",substr(nh_all$Date,nchar(nh_all$Date)-1,nchar(nh_all$Date)),sep=""))

#Spot check number of years represented
summary(hots$Years)
summary(as.factor(hots$Years))
min(summary(as.factor(hots$Years)))
min(summary(as.factor(hots$Years)))
mean(summary(as.factor(hots$Years)))
hist(hots$Years)
summary(bats$Years)
summary(as.factor(bats$Years))
min(summary(as.factor(bats$Years)))
mean(summary(as.factor(bats$Years)))
hist(bats$Years)
summary(as.factor(bbh$Years))
min(summary(as.factor(bbh$Years)))
mean(summary(as.factor(bbh$Years)))
hist(bbh$Years)
summary(as.factor(cml$Years[!is.na(cml$temperature_C)]))

min(summary(as.factor(cml$Years)))
mean(summary(as.factor(cml$Years)))
hist(cml$Years)
summary(as.factor(nh_all$Years))
min(summary(as.factor(nh_all$Years)))
mean(summary(as.factor(nh_all$Years)))

mean(mean(summary(as.factor(hots$Years))),mean(summary(as.factor(bbh$Years))),mean(summary(as.factor(cml$Years))),mean(summary(as.factor(nh_all$Years))))

#Get surface data from Hawaii station, calculate s.d. per year
hots_sd<-hots[hots$`press(dbar)`<10] %>%
  group_by(Years) %>%
  summarise(sd(`temp(ITS-90)`,na.rm=T),sd(ph,na.rm=T))
colnames(hots_sd)<-c("Years","Temp_sd","pH_sd")

#Get surface and tropical data from Bermuda station, calculate s.d. per year
bats_sd<-bats[bats$Depth < 10 & bats$latN<30]  %>%
  group_by(Years,latN) %>%
  summarise(sd(Temp,na.rm=T))
colnames(bats_sd)<-c("Years","Latitude","Temp_sd")

#Calculate s.d. per year for Maine station
bbh_sd<-bbh  %>%
  group_by(Years) %>%
  summarise(sd(`Sea Surface Temp Ave C`,na.rm=T))
colnames(bbh_sd)<-c("Years","Temp_sd")

#Get surface and temperate data from New Hampshire station, calculate s.d. per year
cml_sd<-cml %>%
  group_by(Years) %>%
  summarise(sd(temperature_C,na.rm=T),sd(pH,na.rm=T))
colnames(cml_sd)<-c("Years","Temp_sd","pH_sd")

#Calculate s.d. per year for second New Hampshire station
nh_sd<-nh_all  %>%
  group_by(Years) %>%
  summarise(sd(ph,na.rm=T))
colnames(nh_sd)<-c("Years","pH_sd")

#Explicitly set Latitude for Hawaii and Maine station
hots_sd$Latitude<-22.75
bbh_sd$Latitude<-43.84
cml_sd$Latitude<-43.07
nh_sd$Latitude<-43.02

#Specify data source
hots_sd$Source<-as.factor("HOTS - Hawaii")
bats_sd$Source<-as.factor("Bermuda")
bbh_sd$Source<-as.factor("BBH - Maine")
cml_sd$Source<-as.factor("CML/NOAA - New Hampshire")
nh_sd$Source<-as.factor("New Hampshire")

#Code below specific to SST analysis
#Get tropical subset of data from Reconstructed dataset
dat_1800_Tropic_sd<-dat_1800_new[dat_1800_new$latitude>20 & dat_1800_new$latitude<25 & dat_1800_new$longitude>160 & dat_1800_new$longitude<230]  %>%
  group_by(year,latitude) %>%
  summarise(sd(SST),sd(pH))
colnames(dat_1800_Tropic_sd)<-c("Years", "Latitude", "Temp_sd","pH_sd")
#Specify data source
dat_1800_Tropic_sd$Source<-as.factor("Model")

#Alternate narrower window for tropics which is the same range as the temperate
# dat_1800_Tropic_sd<-dat_1800_new[dat_1800_new$latitude>20 & dat_1800_new$latitude<25 & dat_1800_new$longitude>180 & dat_1800_new$longitude<210]  %>%
#   group_by(year,latitude) %>%
#   summarise(sd(SST),sd(pH))
# colnames(dat_1800_Tropic_sd)<-c("Years", "Latitude", "Temp_sd","pH_sd")
# #Specify data source
# dat_1800_Tropic_sd$Source<-as.factor("Model")

#Get temprate subset of data from Reconstructed dataset
dat_1800_Temprate_sd<-dat_1800_new[dat_1800_new$latitude>40 & dat_1800_new$latitude<45 & dat_1800_new$longitude>290 & dat_1800_new$longitude<320]  %>%
  group_by(year,latitude) %>%
  summarise(sd(SST),sd(pH))
colnames(dat_1800_Temprate_sd)<-c("Years", "Latitude", "Temp_sd","pH_sd")
#Specify data source
dat_1800_Temprate_sd$Source<-as.factor("Model")

#Get tropical subset of data from Projected 4.5 dataset
# dat_4.5_Tropic_sd<-dat_2100_4.5_new[dat_2100_4.5_new$latitude>20 & dat_2100_4.5_new$latitude<25 & dat_2100_4.5_new$longitude>180 & dat_2100_4.5_new$longitude<210]  %>%
#   group_by(year,latitude) %>%
#   summarise(sd(SST),sd(pH))
# colnames(dat_4.5_Tropic_sd)<-c("Years", "Latitude", "Temp_sd","pH_sd")
# #Specify data source
# dat_4.5_Tropic_sd$Source<-as.factor("Projection_4.5")

#Get tropical subset of data from Projected 4.5 dataset
#Alternate narrower window for tropics which is the same range as the temperate
dat_4.5_Tropic_sd<-dat_2100_4.5_new[dat_2100_4.5_new$latitude>20 & dat_2100_4.5_new$latitude<25 & dat_2100_4.5_new$longitude>160 & dat_2100_4.5_new$longitude<230]  %>%
  group_by(year,latitude) %>%
  summarise(sd(SST),sd(pH))
colnames(dat_4.5_Tropic_sd)<-c("Years", "Latitude", "Temp_sd","pH_sd")
#Specify data source
dat_4.5_Tropic_sd$Source<-as.factor("Projection_4.5")

#Get temprate subset of data from Projected 4.5 dataset
dat_4.5_Temprate_sd<-dat_2100_4.5_new[dat_2100_4.5_new$latitude>40 & dat_2100_4.5_new$latitude<45 & dat_2100_4.5_new$longitude>290 & dat_2100_4.5_new$longitude<320]  %>%
  group_by(year,latitude) %>%
  summarise(sd(SST),sd(pH))
colnames(dat_4.5_Temprate_sd)<-c("Years", "Latitude","Temp_sd","pH_sd")
#Specify data source
dat_4.5_Temprate_sd$Source<-as.factor("Projection_4.5")

#Get tropical subset of data from Projected 8.5 dataset
dat_8.5_Tropic_sd<-dat_2100_8.5_new[dat_2100_8.5_new$latitude>20 & dat_2100_8.5_new$latitude<25 & dat_2100_8.5_new$longitude>160 & dat_2100_8.5_new$longitude<230]  %>%
  group_by(year,latitude) %>%
  summarise(sd(SST),sd(pH))
colnames(dat_8.5_Tropic_sd)<-c("Years", "Latitude","Temp_sd","pH_sd")
#Specify data source
dat_8.5_Tropic_sd$Source<-as.factor("Projection_8.5")

#Get temprate subset of data from Projected 8.5 dataset
dat_8.5_Temprate_sd<-dat_2100_8.5_new[dat_2100_8.5_new$latitude>40 & dat_2100_8.5_new$latitude<45 & dat_2100_8.5_new$longitude>290 & dat_2100_8.5_new$longitude<320]  %>%
  group_by(year,latitude) %>%
  summarise(sd(SST),sd(pH))
colnames(dat_8.5_Temprate_sd)<-c("Years", "Latitude", "Temp_sd","pH_sd")
#Specify data source
dat_8.5_Temprate_sd$Source<-as.factor("Projection_8.5")

#Code below specific to SST analysis with longitude included 
#Get tropical subset of data from Reconstructed dataset
dat_1800_Tropic_sd_LL<-dat_1800_new[dat_1800_new$latitude>20 & dat_1800_new$latitude<25 ] %>%#& dat_1800_new$longitude>160 & dat_1800_new$longitude<230]  %>%
  group_by(year,latitude, longitude) %>%
  summarise(sd(SST),sd(pH))
colnames(dat_1800_Tropic_sd_LL)<-c("Years", "Latitude", "Longitude", "Temp_sd","pH_sd")
#Specify data source
dat_1800_Tropic_sd_LL$Source<-as.factor("Model")

#Get temprate subset of data from Reconstructed dataset
dat_1800_Temprate_sd_LL<-dat_1800_new[dat_1800_new$latitude>40 & dat_1800_new$latitude<45 ] %>% #& dat_1800_new$longitude>290 & dat_1800_new$longitude<320]  %>%
  group_by(year,latitude,longitude) %>%
  summarise(sd(SST),sd(pH))
colnames(dat_1800_Temprate_sd_LL)<-c("Years", "Latitude", "Longitude", "Temp_sd","pH_sd")
#Specify data source
dat_1800_Temprate_sd_LL$Source<-as.factor("Model")

#Get tropical subset of data from Projected 4.5 dataset
dat_4.5_Tropic_sd_LL<-dat_2100_4.5_new[dat_2100_4.5_new$latitude>20 & dat_2100_4.5_new$latitude<25 ] %>% #& dat_2100_4.5_new$longitude>160 & dat_2100_4.5_new$longitude<230]  %>%
  group_by(year,latitude,longitude) %>%
  summarise(sd(SST),sd(pH))
colnames(dat_4.5_Tropic_sd_LL)<-c("Years", "Latitude", "Longitude", "Temp_sd","pH_sd")
#Specify data source
dat_4.5_Tropic_sd_LL$Source<-as.factor("Projection_4.5")

#Get temprate subset of data from Projected 4.5 dataset
dat_4.5_Temprate_sd_LL<-dat_2100_4.5_new[dat_2100_4.5_new$latitude>40 & dat_2100_4.5_new$latitude<45 ] %>% #& dat_2100_4.5_new$longitude>290 & dat_2100_4.5_new$longitude<320]  %>%
  group_by(year,latitude,longitude) %>%
  summarise(sd(SST),sd(pH))
colnames(dat_4.5_Temprate_sd_LL)<-c("Years", "Latitude" , "Longitude","Temp_sd","pH_sd")
#Specify data source
dat_4.5_Temprate_sd_LL$Source<-as.factor("Projection_4.5")

#Get tropical subset of data from Projected 8.5 dataset
dat_8.5_Tropic_sd_LL<-dat_2100_8.5_new[dat_2100_8.5_new$latitude>20 & dat_2100_8.5_new$latitude<25 ] %>% #& dat_2100_8.5_new$longitude>160 & dat_2100_8.5_new$longitude<230]  %>%
  group_by(year,latitude,longitude) %>%
  summarise(sd(SST),sd(pH))
colnames(dat_8.5_Tropic_sd_LL)<-c("Years", "Latitude", "Longitude","Temp_sd","pH_sd")
#Specify data source
dat_8.5_Tropic_sd_LL$Source<-as.factor("Projection_8.5")

#Get temprate subset of data from Projected 8.5 dataset
dat_8.5_Temprate_sd_LL<-dat_2100_8.5_new[dat_2100_8.5_new$latitude>40 & dat_2100_8.5_new$latitude<45 ] %>% #& dat_2100_8.5_new$longitude>290 & dat_2100_8.5_new$longitude<320]  %>%
  group_by(year,latitude,longitude) %>%
  summarise(sd(SST),sd(pH))
colnames(dat_8.5_Temprate_sd_LL)<-c("Years", "Latitude", "Longitude", "Temp_sd","pH_sd")
#Specify data source
dat_8.5_Temprate_sd_LL$Source<-as.factor("Projection_8.5")

#Combine all data into one for visualization
#Tdat<-bind_rows(cml_sd[,-3],hots_sd[,c(1:2,4:5)],bbh_sd,dat_4.5_Tropic_sd[,-4],dat_4.5_Temprate_sd[,-4],dat_8.5_Tropic_sd[,-4],dat_8.5_Temprate_sd[,-4],dat_1800_Tropic_sd[,-4],dat_1800_Temprate_sd[,-4])
#Dataset focused on 1970-2020
Tdat<-bind_rows(cml_sd[,-3],hots_sd[,c(1:2,4:5)],bbh_sd,dat_1800_Tropic_sd[dat_1800_Tropic_sd$Years>1830,-4],dat_1800_Temprate_sd[dat_1800_Temprate_sd$Years>1830,-4])
Tdat_LL<-bind_rows(cml_sd[,-3],hots_sd[,c(1:2,4:5)],bbh_sd,dat_1800_Tropic_sd_LL[dat_1800_Tropic_sd_LL$Years>1830,-5],dat_1800_Temprate_sd_LL[dat_1800_Temprate_sd_LL$Years>1830,-5])

#Visualize data
ggplot(data=Tdat,aes(x=Years,y=Temp_sd,color=Latitude, shape=Source))+
  geom_point(size=2)+
  xlab("Years")+
  ylab("Surface Sea Temperature (standard deviation)")+
  geom_hline(yintercept=mean(hots_sd$Temp_sd, na.rm = T),col="red")+
  geom_hline(yintercept=mean(dat_1800_Tropic_sd$Temp_sd, na.rm = T),col="red",linetype="dashed")+
  geom_hline(yintercept=mean(bbh_sd$Temp_sd, na.rm = T),col="cyan")+
  geom_hline(yintercept=mean(dat_1800_Temprate_sd$Temp_sd, na.rm = T),col="cyan",linetype="dashed")+
  #scale_shape_manual(values=c(19,15,3,4,17))+
  scale_colour_continuous(name="Latitude",type = "gradient",low="red",high="cyan")

#Why is 1950 an outlier for BBH? Checking the year with the greatest number of NA's
NAs<-NULL
for(i in min(bbh$Years):max(bbh$Years)){
  NAs<-c(NAs,sum(is.na(bbh[bbh$Years==i,4])))
}
yrs<-min(bbh$Years):max(bbh$Years)
NA_count<-data.frame(yrs,NAs)
NA_count[NA_count$NAs==max(NA_count$NAs),]
#1950 has the greatest number of NA's (285 out of 365, so 78% of the year is missing)

#Why is 2018an outlier for CML? Checking the year with the greatest number of NA's
NAs<-NULL
for(i in min(cml$Years):max(cml$Years)){
  NAs<-c(NAs,sum(is.na(cml[cml$Years==i,6])))
}
yrs<-min(cml$Years):max(cml$Years)
NA_count<-data.frame(yrs,NAs)
NA_count[NA_count$NAs==max(NA_count$NAs),]

counts<-data.frame(summary(as.factor(cml$Years)))
colnames(counts)<-"surveys"
counts$yrs<-rownames(counts)
cml_counts<-merge(counts, NA_count,by="yrs")
cml_counts$prop<-cml_counts$NAs/cml_counts$surveys
max(cml_counts$prop)
cml_counts[cml_counts$prop==min(cml_counts$prop),]

#It's' not the number (raw or proportional) of NA's, but it might be the season in which 2018 was sampled
summary(as.factor(cml$YearsMonths[cml$Years==2018]))
#Yep, only Winter-early Spring was sampled in 2018

omit<-c(1905:1969,1988,2018)
#Removing 1950 and 2018 from dataset and visualizing again
SST<-
  ggplot(data=Tdat[!Tdat$Years%in%omit,],aes(x=Years,y=Temp_sd,color=Latitude, shape=Source))+
  geom_point(size=4)+
  xlab("Years")+
  ylab("Surface Sea Temperature (standard deviation)")+
  xlim(1969,2020)+
  geom_hline(yintercept=mean(hots_sd$Temp_sd, na.rm = T),col="red")+
  geom_hline(yintercept=mean(dat_1800_Tropic_sd[dat_1800_Tropic_sd$Years>1830,-5]$Temp_sd, na.rm = T),col="red",linetype="dashed")+
  geom_hline(yintercept=mean(bbh_sd$Temp_sd, na.rm = T),col="blue")+
  geom_hline(yintercept=mean(dat_1800_Temprate_sd[dat_1800_Temprate_sd$Years>1830,-5]$Temp_sd, na.rm = T),col="blue",linetype="dashed")+
  scale_shape_manual(values=c(17,15,18,21))+
  scale_colour_continuous(name="Latitude",type = "gradient",low="red",high="blue")+
  theme_classic()+
  theme(text = element_text(size = 20))

SST_LL<-
  ggplot(data=Tdat_LL[!Tdat_LL$Years%in%omit,],aes(x=Years,y=Temp_sd,color=Latitude, shape=Source))+
  #ggplot(data=Tdat_LL[Tdat_LL$Latitude<40,],aes(x=Years,y=Temp_sd,color=Latitude, shape=Source))+
  geom_point(size=4)+
  xlab("Years")+
  ylab("Surface Sea Temperature (standard deviation)")+
  xlim(1969,2020)+
  geom_hline(yintercept=mean(hots_sd$Temp_sd, na.rm = T),col="red")+
  geom_hline(yintercept=mean(dat_1800_Tropic_sd_LL[dat_1800_Tropic_sd_LL$Years>1830,-5]$Temp_sd, na.rm = T),col="red",linetype="dashed")+
  geom_hline(yintercept=mean(bbh_sd$Temp_sd, na.rm = T),col="blue")+
  geom_hline(yintercept=mean(dat_1800_Temprate_sd_LL[dat_1800_Temprate_sd_LL$Years>1830,-5]$Temp_sd, na.rm = T),col="blue",linetype="dashed")+
  #scale_shape_manual(values=c(19,15,17,18))+
  scale_shape_manual(values=c(17,15,18,21))+
  scale_colour_continuous(name="Latitude",type = "gradient",low="red",high="blue")+
  theme_classic()+
  theme(text = element_text(size = 20))

#Code below is specific to pH data visualization
#Combine all data into one for visualization
#Tdat<-bind_rows(nh_sd,hots_sd,dat_4.5_Tropic_sd,dat_4.5_Temprate_sd,dat_8.5_Tropic_sd,dat_8.5_Temprate_sd,dat_1800_Tropic_sd,dat_1800_Temprate_sd)
TdatPh<-bind_rows(nh_sd,hots_sd,dat_1800_Tropic_sd[dat_1800_Tropic_sd$Years>1830,],dat_1800_Temprate_sd[dat_1800_Temprate_sd$Years>1830,])
TdatPh_LL<-bind_rows(nh_sd,hots_sd,dat_1800_Tropic_sd_LL[dat_1800_Tropic_sd_LL$Years>1830,],dat_1800_Temprate_sd_LL[dat_1800_Temprate_sd_LL$Years>1830,])

omitPh<-c(1905:1969)
pH<-
  ggplot(data=TdatPh[!TdatPh$Years%in%omitPh,],aes(x=Years,y=pH_sd,color=Latitude, shape=Source))+
  geom_point(size=4,show.legend=FALSE)+
  #geom_point(size=3)+
  xlab("Years")+
  ylab("pH (standard deviation)")+
  xlim(1969,2020)+
  geom_hline(yintercept=mean(hots_sd$pH_sd, na.rm = T),col="red")+
  geom_hline(yintercept=mean(nh_sd$pH_sd, na.rm = T),col="blue")+
  geom_hline(yintercept=mean(dat_1800_Tropic_sd[dat_1800_Tropic_sd$Years>1830,]$pH_sd, na.rm = T),col="red",linetype="dashed")+
  geom_hline(yintercept=mean(dat_1800_Temprate_sd[dat_1800_Temprate_sd$Years>1830,]$pH_sd, na.rm = T),col="blue",linetype="dashed")+
  #scale_shape_manual(values=c(17,3,15))+
  scale_shape_manual(values=c(17,15,21))+
  scale_colour_continuous(name="Latitude",type = "gradient",low="red",high="blue")+
  #guides(
  #  color = guide_colorbar(order = 1),
  #  fill = guide_legend(order = 0)
  #)+
  theme_classic()+
  theme(text = element_text(size = 20))

pH_LL<-ggplot(data=TdatPh_LL[!TdatPh_LL$Years%in%omitPh,],aes(x=Years,y=pH_sd,color=Latitude, shape=Source))+
  geom_point(size=4,show.legend=FALSE)+
  #geom_point(size=3)+
  xlab("Years")+
  ylab("pH (standard deviation)")+
  xlim(1969,2020)+
  geom_hline(yintercept=mean(hots_sd$pH_sd, na.rm = T),col="red")+
  geom_hline(yintercept=mean(nh_sd$pH_sd, na.rm = T),col="blue")+
  geom_hline(yintercept=mean(dat_1800_Tropic_sd_LL[dat_1800_Tropic_sd_LL$Years>1830,]$pH_sd, na.rm = T),col="red",linetype="dashed")+
  geom_hline(yintercept=mean(dat_1800_Temprate_sd_LL[dat_1800_Temprate_sd_LL$Years>1830,]$pH_sd, na.rm = T),col="blue",linetype="dashed")+
  #scale_shape_manual(values=c(17,3,15))+
  scale_shape_manual(values=c(17,15,21))+
  scale_colour_continuous(name="Latitude",type = "gradient",low="red",high="blue")+
  #guides(
  #  color = guide_colorbar(order = 1),
  #  fill = guide_legend(order = 0)
  #)+
  theme_classic()+
  theme(text = element_text(size = 20))

legend <- get_legend(SST)
SST<-SST+ theme(legend.position="none")

legend_LL <- get_legend(SST_LL)
SST_LL<-SST_LL+ theme(legend.position="none")

library(gridExtra)
grid.arrange(SST,pH, legend,ncol=3,nrow=1,widths=c(2,2, 1))

grid.arrange(SST_LL,pH_LL, legend_LL,ncol=3,nrow=1,widths=c(2,2, 1))

# tmp <- arrangeGrob(SST, SST_LL, pH, pH_LL, layout_matrix = matrix(c(1,2,3,4), nrow = 2))
# grid.arrange(tmp,legend,ncol=2,nrow=1,widths=c(4, 1))





