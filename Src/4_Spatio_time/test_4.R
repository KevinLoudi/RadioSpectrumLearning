# Author: Kevin
# Date: 7th March, 2017
# Propose: test scripts for spatial-temporal
# Environment: R 3.3
####################################

##############test scripts for class 'SDSTF'#############
# create object of 'SDSTF'
library(R.matlab)
graphics.off()
path<-"D:/Code/WorkSpace/ThesisCode/Src/4_Spatio_time/SpatialDataset/Spdata_1710_1740.mat"
sp_data<-readMat(path)
path<-"D:/Code/WorkSpace/ThesisCode/Src/4_Spatio_time/SpatialDataset/Splocation_1710_1740.mat"
sp_loc<-readMat(path)
path<-"D:/Code/WorkSpace/ThesisCode/Src/4_Spatio_time/SpatialDataset/Sptime_1710_1740.mat"
sp_ti<-readMat(path)

sp_data<-sp_data$data.sp
sp_loc<-sp_loc$locateion.sp
sp_ti<-as.POSIXct(sp_ti$dateStamp)

grid=cbind(x=sp_loc[1,1:9],y=sp_loc[2,1:9])
row.names(grid) = paste("point", 1:nrow(grid), sep="")
grid=SpatialPoints(grid)
plot(grid@coords)
library(xts)
len=length(sp_ti)
grid_time=xts(1:len,as.POSIXct(sp_ti))
plot(grid_time)
m_g=c(1:9)
#m_data=rnorm(length(grid)*length(grid_time),mean=rep(m_g, len))
m_data<-as.vector(sp_data)
plot(m_data)
IDs = paste("ID",1:length(m_data)) # ID for each data point
mydata = data.frame(values = signif(m_data,9), ID=IDs) # 12 observers
grid_df = STFDF(grid, grid_time, mydata) #create a 'SDSTF' object
plot(grid_df)



grid_df
grid_df@time
grid_df@time
grid_df@endTime
grid_df@sp@coords
grid_df@sp@bbox

sp = cbind(x = c(0,0,1), y = c(0,1,1)) #combine two vector by columns
row.names(sp) = paste("point", 1:nrow(sp), sep="")
library(sp)
library(spacetime)

sp = SpatialPoints(sp)#SPATIAL CLASS #three location positions
plot(sp@coords)

library(xts)
time = xts(1:4, as.POSIXct("2010-08-05")+3600*(10:13)) #four time slots
plot(time)

m = c(10,20,30) # means for each of the 3 point locations #1*3 column

mydata = rnorm(length(sp)*length(time),mean=rep(m, 4)) #4*3 data point

IDs = paste("ID",1:length(mydata)) # ID for each data point
mydata = data.frame(values = signif(mydata,3), ID=IDs) # 12 observers

stfdf = STFDF(sp, time, mydata) #create a 'SDSTF' object
plot(stfdf)


stfdf
stfdf@data #data points
stfdf@time #time slots
stfdf@sp@coords #grids
stfdf@sp@bbox  #range

stsdf = as(stfdf, "STSDF")

stsdf[1:2,] #first two locations
stsdf[,1:2] #first two time-slots
stsdf[,,2]
stsdf[,,"values"]

stsdf[1,]
stsdf[,2]
# examples for [[, [[<-, $ and $<- 
stsdf[[1]]
stsdf[["values"]]
stsdf[["newVal"]] <- rnorm(12)
stsdf$ID
stsdf$ID = paste("OldIDs", 1:12, sep="")
stsdf$NewID = paste("NewIDs", 12:1, sep="")
stsdf
x = stsdf[stsdf,]
x = stsdf[stsdf[1:2,],]
all.equal(x, stsdf[1:2,])


####################ST class
# PLEASE read the vignette of package spacetime for a more
# clever way to do all this!
library(sp)
library(gstat)
library(rgdal)
library(maptools)
# load wind data, run test:
example(wind)

m = map2SpatialLines(
  map("worldHires", xlim = c(-11,-5.4), ylim = c(51,55.5), plot=F))
proj4string(m) = "+proj=longlat +datum=WGS84 +ellps=WGS84"
m = spTransform(m, CRS("+proj=utm +zone=29 +datum=WGS84 +ellps=WGS84"))

# model temporal autocorrelation
acf(wind[7])
tdiscr = 0:40
lines(tdiscr, exp(- tdiscr/1.5))

# set up data, last year
years = 61
months = 1
jday = c(1,6,11,16,21,26)
sel = wind[wind$year %in% years & 
             wind$month %in% months &
             wind$jday %in% jday,]

#stations = 4:15
stations = 4:15

sels = stack(sel[stations])
sels$t = rep(sel$jday, length(stations))
sels$x = coordinates(wind.loc)[match(sels$ind, wind.loc$Code),1]
sels$y = coordinates(wind.loc)[match(sels$ind, wind.loc$Code),2]
summary(sels)

coordinates(sels) = ~x+y
proj4string(sels) = "+proj=longlat +datum=WGS84 +ellps=WGS84"
sels = spTransform(sels, CRS("+proj=utm +zone=29 +datum=WGS84 +ellps=WGS84"))
grd = makegrid(m, n = 1000)
grd$t = rep(1, nrow(grd))
coordinates(grd) = ~x1+x2
gridded(grd)=TRUE
proj4string(grd) = proj4string(sels)

#sels = as(sels, "data.frame")

# setup grid
covfn = function(x, y = x) { 
  u = spDists(coordinates(x), coordinates(y))
  t = abs(outer(x$t,y$t,"-"))
  0.6 * exp(-u/750000) * exp(-t/1.5)
}
for (i in 1:120) {
  grd$t = rep(i/4, nrow(grd))
  n = paste("t", i/4, sep="")
  grd[[n]] = krige0(sqrt(values)~1, sels, grd, covfn)
}
grd$t = NULL
#grd$pr = out$pred
#library(lattice)
#levelplot(pr~x1+x2|t,grd,col.regions=bpy.colors())
spl = list(list("sp.points", sels,first=F, cex=.5),
           list("sp.lines", m, col='grey'))
spplot(grd, sp.layout = spl, col.regions=bpy.colors())