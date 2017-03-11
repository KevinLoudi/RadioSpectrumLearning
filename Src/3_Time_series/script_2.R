# Author: Kevin
# Date: 9th March, 2017
# Propose: matlab script_1.m's counterpart
# Environment: R 3.3
####################################

#%%%%%%%%%%%%%%%%%%%%%%%%%Part 1: Index sequence analysis
#load data
library(R.matlab)
path<-"D:/Code/WorkSpace/ThesisCode/Src/1_Generate_Dataset/ChannelStauseDataset/Duty_cycle_1710_1740.mat"
dc<-readMat(path)
path<-"D:/Code/WorkSpace/ThesisCode/Src/1_Generate_Dataset/ChannelStauseDataset/Timeindex_1710_1740.mat"
ti_ix<-readMat(path)
u <- runif(2000,0,1)

#translate into time series
ix=ts(as.POSIXct(ti_ix$tot.time))
dc=ts(dc$dc)


Box.test(dc,lag = 20,type="Ljung-Box")
Box.test(u,lag = 20,type="Ljung-Box")

library("fUnitRoots")
adfTest(dc,lags=10,type = "nc",title=NULL,description = NULL)
adfTest(u,lags=10,type = "nc",title=NULL,description = NULL)

#differ to reach stationary
dc.dif<-diff(dc,differences=1)
plot.ts(dc.dif)
adfTest(dc.dif,lags=10,type = "nc",title=NULL,description = NULL)

#lag p exceeds the significance bounds p=2
acf(dc.dif,lag.max=10,main='Autocorrect function')
#slowly decreasing in magnitude with increasing lag q=10
pacf(dc.dif,lag.max=30,main='Partial autocorrect function')

#AIC and BIC (1,1,2)
require(forecast)
model_par<-forecast::auto.arima(dc.dif,ic="bic")
model_par

#estimate model parameters
dc.arima<-arima(dc,order=c(1,1,2))
dc.arima
dc.forecast<-forecast.Arima(dc.arima,h=20,level=c(99.5))
plot(dc.forecast)

#fit error analysis
plot_forecast_error(dc.arima$residuals)
plot_forecast_error(dc.forecast$residuals)

##############test scripts for class 'SDSTF'#############
# prepare library
library(R.matlab)
library(xts)
library(spacetime)
require(sp) 
require(gstat) #load/install + load installr
show_fig<-FALSE #do not show unimportant figures
show_details<-TRUE #do show details in console

#load data set
origin_path<-"D:/Code/WorkSpace/ThesisCode/Src/4_Spatio_time/SpatialDataset/"
start_freq<-1068 #kHz
stop_freq<-1068
path<-paste(origin_path,"Spdata_1068_1068.mat",sep = "") #data matrix
sp_data<-readMat(path,open='r')
path<-paste(origin_path,"Splocation_1068_1068.mat",sep = "") #location set
sp_loc<-readMat(path,open='r')
path<-paste(origin_path,"Sptime_1068_1068.mat",sep = "") #time stamp
sp_ti<-readMat(path,open='r')
path<-paste(origin_path,"SensorIds.mat",sep = "") #sensor/device names
sp_de<-readMat(path,open='r')

#reserve partial data
start_obs=500
stop_obs=700  #cut off by index
sp_data<-sp_data$data.sp[start_obs:stop_obs,] #data list
sp_loc<-sp_loc$locateion.sp #location list
sp_ti<-as.POSIXct(sp_ti$dateStamp[start_obs:stop_obs]) #transform to standard time stampe
sp_de<-sp_de$SensorIds #sensor list
de_num<-length(sp_de) #sensor num


#create spatial grid as coordinates
grid=cbind(x=sp_loc[1,1:de_num],y=sp_loc[2,1:de_num])
row.names(grid) = paste("point", 1:nrow(grid), sep="")
#create spatialpoints object 'grid'
grid=SpatialPoints(grid)
if (show_fig) {plot(grid@coords)}
  
#create extensible time stamp 'grid_time'
len=length(sp_ti)
grid_time=xts(1:len,as.POSIXct(sp_ti))
plot(grid_time)
m_g=c(1:de_num)

#assign ID for data point
m_data<-as.vector(sp_data)
plot(m_data)
IDs = paste("ID",1:length(m_data)) # ID for each data point

#organize data in frame
mydata = data.frame(values = signif(m_data,de_num), ID=IDs) # 12 observers

#create Spatial-temporal object
grid_df = STFDF(grid, grid_time, mydata) #create a 'SDSTF' object
grid_df=as(grid_df,"STSDF")

##statistics for data set
grid_dd<-grid_df@data
if (show_fig) {plot(grid_dd)}

#summary observations of each station
barplot(sort(table(grid_df@index[,1])),
        main="reported days per station",
        ylab="number of days", xaxt="n")

var(grid_df@data$values)

#numbers of stations
length(grid_df@sp)

#calculate the empirical variogram: the first important step
#formulation // exclude NaN in calculation
empVgm <- variogramST(values~1, grid_df, tlags=0:10,na.omit=TRUE)
if (show_fig) {plot(empVgm)}
#plot 3-D varigram map
plot(empVgm, wireframe=T, scales=list(arrows=F))
if (show_details) {empVgm}

# fit of theoretical purely spatial models 
spEmpVgm <- empVgm[empVgm$timelag == 0,] #get the lag0 in time
class(spEmpVgm) <- c("gstatVariogram", "data.frame")
spEmpVgm <- spEmpVgm[-1,1:3]
spEmpVgm$dir.hor <- 0
spEmpVgm$dir.ver <- 0

#combined two model: nugget model + exponential model
spVgmMod <- fit.variogram(spEmpVgm, vgm(800,"Exp",300000,20),fit.sills = TRUE,
                          fit.ranges = TRUE,fit.method = 7,debug.level = 1,
                          warn.if.neg = FALSE,fit.kappa = FALSE)
if (show_details) {print(spVgmMod,digit=3)} 

plot(spEmpVgm, spVgmMod)

#estimation the ansiotrpy
linStAni <- estiStAni(empVgm, c(50000,200000))

#demo variogram-distance relationship
plot(gamma~dist, empVgm[empVgm$timelag == 0,], ylim=c(0,100), xlim=c(0,800000))
points(empVgm[empVgm$spacelag == 0,]$timelag*linStAni, empVgm[empVgm$spacelag == 0,]$gamma, col="red")


#separableModel for ST kriging
separableModel <- vgmST("separable", 
                        space=vgm(0.9,"Exp", 200, 0.1),
                        time =vgm(0.9,"Sph", 3.5, 0.1),
                        sill=120)
print(separableModel)

#fit the ST model
fitSepModel <- fit.StVariogram(empVgm, separableModel, fit.method = 7, 
                               stAni = linStAni, method = "L-BFGS-B", 
                               control = list(parscale=c(100,1,10,1,100)),
                               lower = c(10,0,.1,0,0.1), 
                               upper = c(2000,1,12,1,200))

attr(fitSepModel, "optim.output")$value

plot(empVgm, fitSepModel, wireframe=T, all=T, scales=list(arrows=F), zlim=c(0,135))
