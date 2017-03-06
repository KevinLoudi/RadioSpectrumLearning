# Author: Kevin
# Date: 5th March, 2017
# Propose: arima sim and test
# Environment: R 3.3
####################################

#load data
dc <- read.csv(file="D:/Code/WorkSpace/ThesisCode/Src/3_Time_series/TimeSeriesDataset/dc_2_1740_1750.csv", header=FALSE)
#generate random noise
u <- runif(200,0,1)

data<-ts(dc$V1,frequency = 3)
plot.ts(data)
data.log<-log(data)
plot.ts(data.log)

library("TTR")
#smooth with moving average of 8
#provide a trend 
smo<-SMA(dc$V1,n=10)
plot.ts(smo)

#sensoal decompose
dc_compent <- decompose(data)
plot(dc_compent)
dc_sensonal_adj<-data-dc_compent$seasonal
plot(dc_sensonal_adj)

#simple exponential smoothing
dc_forecasts<-HoltWinters(dc$V1,beta=FALSE, gamma=FALSE)
#value of "alpha" is not that close to zero, base on
#  more time-faraway results
dc_forecasts
#forcast results
dc_forecasts$fitted
plot(dc_forecasts)
#accuracy of the forecasts: sum of the squared error
dc_forecasts$SSE

#make forecast for further time points
library("forecast")
#pass in predictive model fitted, h=how many further points
dc_forecasts_further<-forecast.HoltWinters(dc_forecasts,h=10)
#show 80% and 95% prediction interval
plot.forecast(dc_forecasts_further)
plot.ts(dc_forecasts_further$residuals)
Box.test(dc_forecasts_further$residuals,lag=20,type = "Ljung-Box")


