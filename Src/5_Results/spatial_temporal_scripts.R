# Author: Kevin
# Date: 20th March, 2017
# Propose: spatial-temporal related scripts
# Environment: R 3.3
####################################

# prepare library and status
library(R.matlab)
library(xts)
library(spacetime)
require(sp) 
require(gstat) #load/install + load installr

show_fig<-FALSE #do not show unimportant figures
show_details<-TRUE #do show details in console

  #global functions
    #normalize elements of matrix
    Normalize<-function(x)
    {
      x<-(x-min(x,na.rm=TRUE))/(max(x,na.rm=TRUE)-min(x,na.rm=TRUE))
    }
#===================================================

#load data set
origin_path<-"D:/Code/WorkSpace/ThesisCode/Src/4_Spatio_time/SpatialDataset/"

path<-paste(origin_path,"Spdata_1730_1740.mat",sep = "") #data matrix
sp_data<-readMat(path,open='r')
path<-paste(origin_path,"Splocation_1730_1740.mat",sep = "") #location set
sp_loc<-readMat(path,open='r')
path<-paste(origin_path,"Sptime_1730_1740.mat",sep = "") #time stamp
sp_ti<-readMat(path,open='r')
path<-paste(origin_path,"SensorIds3.mat",sep = "") #sensor/device names
sp_de<-readMat(path,open='r')

#===================================================

#reduce data: sp_data,sp_loc,sp_ti,sp_de
start_obs=250
stop_obs=2200  #cut off by index
sp_data<-sp_data$data.sp[start_obs:stop_obs,] #data list
sp_loc<-sp_loc$locateion.sp #location list
sp_ti<-as.POSIXct(sp_ti$dateStamp[start_obs:stop_obs]) #transform to standard time stampe
sp_de<-sp_de$SensorIds #sensor list
de_num<-length(sp_de) #sensor num
#===================================================

#create spatial-temporal object

  #normalize loaction 
x=sp_loc[1,1:de_num]
x_nor=(x-min(x))/(max(x)-min(x))
y=sp_loc[2,1:de_num]
y_nor=(y-min(y))/(max(y)-min(y))

  #create spatial grid
grid=cbind(x_nor,y=y_nor)
row.names(grid) = paste("point", 1:nrow(grid), sep="")
grid=SpatialPoints(grid)

  #create extensible time stamp 'grid_time'
len=length(sp_ti)
grid_time=xts(1:len,as.POSIXct(sp_ti))
plot(grid_time)
m_g=c(1:de_num)

  #organize data in frame
m_data<-as.vector(sp_data)
IDs = paste("ID",1:length(m_data)) # ID for each data point
mydata = data.frame(values = signif(m_data,de_num), ID=IDs) # 12 observers

  #create Spatial-temporal object
library(sp)
grid_df=STFDF(grid, grid_time, mydata) #create a 'SDSTF' object
grid_df=as(grid_df,"STSDF")
#===================================================

#create spatial object
  #select the fisrt time data
s_data<-data.frame(value=mydata$values[1:de_num],ID=mydata$ID[1:de_num])
s_data$value[30]<-0.30;
#s_data$value<-(s_data$value-min(s_data$value,na.rm=TRUE))/(max(s_data$value,na.rm=TRUE)-min(s_data$value,na.rm=TRUE))
#s_data$value<-Normalize(s_data$value)
s_obj<-SpatialPointsDataFrame(grid,s_data)
bubble(s_obj,na.rm=FALSE,main = "SCR的空间分布（2015-12-16 17:35）")
#===================================================

#spatial analysis
  #summarize dataset
grid_df@data$values<-Normalize(grid_df@data$values)
var(grid_df@data$values,na.rm=TRUE)
length(grid_df@sp)

  #calculate the empirical variogram
  save_path="D:/Code/RData/"
  empVgm <- variogramST(values~1, grid_df, tlags=0:10,na.omit=TRUE)
  empVgm$dist  <- empVgm$dist*1000
  empVgm$avgDist  <- empVgm$avgDist*1000
  empVgm$spacelag <- empVgm$spacelag*1000
  
  linStAni <- linStAni/1000
  save(empVgm,file = paste(save_path,"empVgm.RData",sep=""))
  #plot 3-D varigram map
  plot(empVgm, wireframe=T, scales=list(arrows=F))
  #empVgm
  plot(empVgm)

  #fit the 0-lag empirical variogram to purely theoretical spatial models
  spEmpVgm <- empVgm[empVgm$timelag == 0,] #get the lag0 in time
  class(spEmpVgm) <- c("gstatVariogram", "data.frame")
  spEmpVgm <- spEmpVgm[-1,1:3]
  spEmpVgm$dir.hor <- 0
  spEmpVgm$dir.ver <- 0
  theVgm<-vgm(1, "Mat", 1, kappa=.3)
  spVgmMod <- fit.variogram(spEmpVgm, vgm(1, "Mat", 1, kappa=.3),fit.sills = TRUE,
                            fit.ranges = TRUE,fit.method = 7,debug.level = 1,
                            warn.if.neg = FALSE,fit.kappa = TRUE)
  print(spVgmMod,digit=3)
  plot(spEmpVgm, spVgmMod)
  
  linStAni <- estiStAni(empVgm, c(50000,200000))
  plot(gamma~dist, empVgm[empVgm$timelag == 0,], ylim=c(0,100), xlim=c(0,800000))
  points(empVgm[empVgm$spacelag == 0,]$timelag*linStAni, empVgm[empVgm$spacelag == 0,]$gamma, col="red")
#===================================================
  
#spatial-temporal analysis
  separableModel <- vgmST("separable", 
                          space=vgm(0.9,"Exp", 200, 0.1),
                          time =vgm(0.9,"Sph", 3.5, 0.1),
                          sill=120)
  
  fitSepModel <- fit.StVariogram(empVgm, separableModel, fit.method = 7, 
                                 stAni = linStAni, method = "L-BFGS-B", 
                                 control = list(parscale=c(100,1,10,1,100)),
                                 lower = c(10,0,.1,0,0.1), 
                                 upper = c(2000,1,12,1,200))
  
  attr(fitSepModel, "optim.output")$value
  
  plot(empVgm, fitSepModel, wireframe=T, all=T, scales=list(arrows=F), zlim=c(0,135))
  
  krigeST(data~1,data = grid_df, newdata = grid_pred, fitSepModel, nmax = 50)
  
  
#===============================Part 3: Spatial kriging===================================
  
  library(sp)
  library(gstat)
  library(R.matlab)
  
  #load data
  origin_path<-"D:/Code/WorkSpace/ThesisCode/Src/5_Results/Datasets/"
  path<-paste(origin_path,"Data.mat",sep = "") #data set
  sp_data<-readMat(path,open='r')
  path<-paste(origin_path,"Position(2015).mat",sep = "") #location set
  sp_loc<-readMat(path,open='r')
  sp_ti="2015-12-15 19:00:00"
  sp_de<-sp_loc$ID
  

  #reduce data: sp_data,sp_loc,sp_ti,sp_de
  start_obs=1
  stop_obs=1  #cut off by index
  sp_data<-sp_data$shown.occ[start_obs:stop_obs,] #data list
  sp_loc<-[sp_loc$Lon,sp_loc$Lat] #location list
  sp_ti<-as.POSIXct(sp_ti[start_obs:stop_obs]) #transform to standard time stampe
  sp_de<-sp_de$SensorIds #sensor list
  de_num<-length(sp_de) #sensor num

 
  