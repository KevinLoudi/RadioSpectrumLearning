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
# create object of 'STSTF'
library(R.matlab)
graphics.off()
path<-"D:/Code/WorkSpace/ThesisCode/Src/4_Spatio_time/SpatialDataset/Spdata_1710_1740.mat"
sp_data<-readMat(path)
path<-"D:/Code/WorkSpace/ThesisCode/Src/4_Spatio_time/SpatialDataset/Splocation_1710_1740.mat"
sp_loc<-readMat(path)
path<-"D:/Code/WorkSpace/ThesisCode/Src/4_Spatio_time/SpatialDataset/Sptime_1710_1740.mat"
sp_ti<-readMat(path)


start_obs=500
stop_obs=700
sp_data<-sp_data$data.sp[start_obs:stop_obs,]
sp_loc<-sp_loc$locateion.sp
sp_ti<-as.POSIXct(sp_ti$dateStamp[start_obs:stop_obs])

grid=cbind(x=sp_loc[1,1:9],y=sp_loc[2,1:9])
row.names(grid) = paste("point", 1:nrow(grid), sep="")
grid=SpatialPoints(grid)
plot(grid@coords)
library(xts)
len=length(sp_ti)
grid_time=xts(1:len,as.POSIXct(sp_ti))
plot(grid_time)
m_g=c(1:9)

m_data<-as.vector(sp_data)
plot(m_data)
IDs = paste("ID",1:length(m_data)) # ID for each data point
mydata = data.frame(values = signif(m_data,9), ID=IDs) # 12 observers
grid_df = STFDF(grid, grid_time, mydata) #create a 'SDSTF' object

grid_df=as(grid_df,"STSDF")

##statistics for data set
grid_dd<-grid_df@data
#plot(grid_dd)
#summary observations of each station
barplot(sort(table(grid_df@index[,1])),
        main="reported days per station",
        ylab="number of days", xaxt="n")

acf(grid_df[sample(9,1),,drop=T]@data)
var(grid_df@data$values)

#numbers of stations
length(grid_df@sp)

#calculate the empirical variogram: the first important step
#formulation // exclude NaN in calculation
empVgm <- variogramST(values~1, grid_df, tlags=0:10,na.omit=TRUE)
plot(empVgm)
empVgm

# fit of theoretical purely spatial models #
spEmpVgm <- empVgm[empVgm$timelag == 0,]
class(spEmpVgm) <- c("gstatVariogram", "data.frame")
spEmpVgm <- spEmpVgm[-1,1:3]
spEmpVgm$dir.hor <- 0
spEmpVgm$dir.ver <- 0

#combined two model: nugget model + exponential model
spVgmMod <- fit.variogram(spEmpVgm, vgm(800,"Exp",300000,20),fit.sills = TRUE,
                          fit.ranges = TRUE,fit.method = 7,debug.level = 1,
                          warn.if.neg = FALSE,fit.kappa = FALSE)
print(spVgmMod,digit=3) 

plot(spEmpVgm, spVgmMod)

#estimation the ansiotrpy
linStAni <- estiStAni(empVgm, c(50000,200000))


