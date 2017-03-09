# Author: Kevin
# Date: 9th March, 2017
# Propose: matlab script_1.m's counterpart
# Environment: R 3.3
####################################

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

#AIC and BIC 
require(forecast)
model_par<-forecast::auto.arima(dc.dif,ic="bic")
model_par