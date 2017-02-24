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
setwd("...")

#load the data set
data <- read.table("D:\\Code\\WorkSpace\\ThesisCode\\R_src\\ozon_tram1_14102011_14012012.csv", sep=",", header=T)

#demo copy ability
paste(data$generation_time[1])

#transform the UNIX time from numerical to character format
#demo
as.POSIXlt(as.numeric(substr(paste(data$generation_time[1]), start=1, stop=10)), origin="1970-01-01")
#work
data$TIME <- as.POSIXlt(as.numeric(substr(paste(data$generation_time), 1, 10)), origin="1970-01-01")

#extract only numbers needed and converting this
#format into degree
data$LAT <- as.numeric(substr(paste(data$latitude),1,2))+(as.numeric(substr(paste(data$latitude),3,10))/60)
data$LON <- as.numeric(substr(paste(data$longitude),1,1))+(as.numeric(substr(paste(data$longitude),2,10))/60) 
#excluede the NAs
data <- na.omit(data) 

#find the date start and stop point
#min(data$TIME)
#max(data$TIME)

#set a subset to reduce computations
sub <- data[data$TIME>=as.POSIXct('2011-12-12 00:00 CET')&data$TIME<=as.POSIXct('2011-12-13 23:00 CET'),]
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
#two arguments: 1)matrix with coordinates of each point
#2)set the projection in UTM
ozoneSP <- SpatialPoints(ozone.UTM@coords,CRS("+init=epsg:3395"))

#check if we have duplicate positions
#to prevent kriging failing due to "singular matrix"
dupl <- zerodist(ozoneSP) 

#remove spatial duplicates 
ozoneDF <- data.frame(PPB=ozone.UTM$ozone_ppb[-dupl[,2]]) 

#remove temporal duplicates
ozoneTM <- as.POSIXct(ozone.UTM$TIME[-dupl[,2]],tz="CET")

#combine the objects to cmpute the variogram
#and perform the spatiao-temporal interpolation
timeDF <- STIDF(ozoneSP,ozoneTM,data=ozoneDF)
stplot(timeDF) 

#computation of the variogram
#arguments: formula, data source, time units/time lags,
#takes quite a long time
#arguments: formula, data
var <- variogramST(PPB~1,data=timeDF,tunit="hours",assumeRegular=F,progress = interactive(),na.omit=T) 

#plot as series
plot(var,map=F) 
#plot as map
plot(var,map=T) 
#plot as 3d wireframe
plot(var,wireframe=T) 

#fit a model to the variogram
demo(stkrige) 

###############################################

#5 options regarding the variogram models
#  seprarable, product sum, metric, sum metric, 
#  simple sum metric

#set model parametric, automaticly fitting
sumMetric <- vgmST("sumMetric", space = vgm(psill=5,"Sph", range=500, nugget=0),time = vgm(psill=500,"Sph", range=500, nugget=0), joint = vgm(psill=1,"Sph", range=500, nugget=10), stAni=500) 
sumMetric_Vgm <- fit.StVariogram(var, sumMetric, method="L-BFGS-B",lower=pars.l,upper=pars.u,tunit="hours")
attr(sumMetric_Vgm, "MSE")

#set prodiction grids
roads <- shapefile("VEC25_str_l_Clip/VEC25_str_l.shp")
#extract road objects
Klass1 <- roads[roads$objectval=="1_Klass",] 
#change the projection type from CH93 to UTM
Klass1.UTM <- spTransform(Klass1,CRS("+init=epsg:3395"))
#crop the file obtain road in the studies area
Klass1.cropped <- crop(Klass1.UTM,ozone.UTM) 
#show the results
plot(Klass1.cropped)
plot(ozone.UTM,add=T,col="red")
# create a random grid of points along the road lines
sp.grid.UTM <- spsample(Klass1.cropped,n=1500,type="random") 
# create temporal compent
tm.grid <- seq(as.POSIXct('2011-12-12 06:00 CET'),as.POSIXct('2011-12-14 09:00 CET'),length.out=5)
grid.ST <- STF(sp.grid.UTM,tm.grid)

#perform the interpolation
pred <- krigeST(PPB~1, data=timeDF, modelList=sumMetric_Vgm, newdata=grid.ST) 
vignette("st", package = "gstat")
demo(stkrige)