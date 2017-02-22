# Author: Kevin
# Date: 22th Feb, 2017
# Propose: test space-time modeling with gstat 
# Environment: R 3.3
####################################

library(gstat)
library(sp)
library(spacetime)
library(raster)
library(rgdal)
library(rgeos) 

#set the work directory
#setwd("...")

#load the data set
#data <- read.table("D:\\Code\\WorkSpace\\ThesisCode\\R_src\\ozon_tram1_14102011_14012012.csv", sep=",", header=T)

#demo copy ability
paste(data$generation_time[1])

#transform the UNIX time from numerical to character format
#demo
as.POSIXlt(as.numeric(substr(paste(data$generation_time[1]), start=1, stop=10)), origin="1970-01-01")
#work
data$TIME <- as.POSIXlt(as.numeric(substr(paste(data$generation_time), 1, 10)), origin="1970-01-01")


#transform geographical coordinates to degree
data$LAT <- as.numeric(substr(paste(data$latitude),1,2))+(as.numeric(substr(paste(data$latitude),3,10))/60)
data$LON <- as.numeric(substr(paste(data$longitude),1,1))+(as.numeric(substr(paste(data$longitude),2,10))/60) 
#excluede the NAs
data <- na.omit(data) 
#find the date start and stop point
min(data$TIME)
max(data$TIME)

#put aside a subset
sub <- data[data$TIME>=as.POSIXct('2011-12-12 00:00 CET')&data$TIME<=as.POSIXct('2011-12-14 23:00 CET'),]
nrow(sub) #size of the sub-dataset

#Create a SpatialPointsDataFrame
coordinates(sub)=~LON+LAT
projection(sub)=CRS("+init=epsg:4326")

#Transform into Mercator Projection
ozone.UTM <- spTransform(sub,CRS("+init=epsg:3395")) 
#now there exists a object ozone.UTM, which is 
#coordinates with metres

###########################################
#create spatial points object with locations of 
#sensors at any given time
ozoneSP <- SpatialPoints(ozone.UTM@coords,CRS("+init=epsg:3395"))

#check if we have duplicate positions
#to prevent kriging failing due to "singular matrix"
dupl <- zerodist(ozoneSP) 
#create the  temporal component 
ozoneDF <- data.frame(PPB=ozone.UTM$ozone_ppb[-dupl[,2]]) 
ozoneTM <- as.POSIXct(ozone.UTM$TIME[-dupl[,2]],tz="CET")
timeDF <- STIDF(ozoneSP,ozoneTM,data=ozoneDF)
stplot(timeDF) 

#computation of the variogram
var <- variogramST(PPB~1,data=timeDF,tunit="hours",assumeRegular=F,na.omit=T) 
#plot as series
plot(var,map=F) 
#plot as map
plot(var,map=T) 
#plot as 3d wireframe
plot(var,wireframe=T) 
