# Author: Kevin
# Date: 5th March, 2017
# Propose: arima modeling
# Environment: R 3.3
####################################

#load data
dc <- read.csv(file="D:/Code/WorkSpace/ThesisCode/Src/3_Time_series/TimeSeriesDataset/dc_2_1740_1750.csv", header=FALSE)
u <- runif(200,0,1)
data<-ts(dc$V1)

#difference the time series
# until the mean value appear stationary
data.dif<-diff(data,differences=3)
plot.ts(data.dif)
#make sure the differenced series is stationary and can pass adf test
adfTest(dc$V1,lags=10,type = "nc",title=NULL,description = NULL)

#select model order for ARIMA with the help of the correlogram and partial
# correlogram
#lag p exceeds the significance bounds p=2
acf(data.dif,lag.max=20)
#slowly decreasing in magnitude with increasing lag q=10
pacf(data.dif,lag.max=30)








