# Author: Kevin
# Date: 7th March, 2017
# Propose: test scripts for spatial-temporal
# Environment: R 3.3
####################################

##############test scripts for class 'SDSTF'#############
# create object of  'SDSTF'
sp = cbind(x = c(0,0,1), y = c(0,1,1)) #combine two vector by columns
row.names(sp) = paste("point", 1:nrow(sp), sep="")
library(sp)

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
plot(sp)

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