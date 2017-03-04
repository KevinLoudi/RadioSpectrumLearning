# Author: Kevin
# Date: 4th March, 2017
# Propose: time series analysi in dc and scr
# Environment: R 3.3
####################################

dc <- read.csv(file="D:/Code/WorkSpace/ThesisCode/Src/3_Time_series/TimeSeriesDataset/dc_2_1740_1750.csv", header=FALSE)
u <- runif(200,0,1)
dc.dif=diff(dc$V1,lag=1,differences = 1)

par(mar = rep(2, 4))
acf(dc,main="acf of duty cycle")
acf(dc.dif,main="acf of diff of duty cycle")
acf(u,main="acf of random noise")
len=length(dc)
plot(dc$V1,type="s",main="duty cycle")
plot(u,type="s",main="random noise")

#to see if the series auto-correlated
Box.test(dc,lag=10,type=c("Ljung-Box"))
Box.test(dc.dif,lag=10,type=c("Ljung-Box"))
Box.test(u,lag=10,type=c("Ljung-Box")) #high than 5%, reject the correlation presume

#decompose of time series
dc.decom<-decompose(dc,type="multiplicative")
trend<-dc.decom$trend
sensonal<-dc.decom$sensonal
random<-dc.decom$random

dc.dif<-diff(dc$V1)
par(mfrow=c(2,1),cex.axis=1.5,cex.lab=1.5)
plot(dc$V1,type="o",xlab="time slots")
plot(dc.dif,type="o",xlab="time slots")

par(mfrow=c(1,1))
par(mfrow=c(2,1),cex.axis=1.5,cex.lab=1.5)
acf(dc.dif)
pacf(dc.dif)
par(mfrow=c(1,1))

dc.fit=sarima(dc,0,1,1,detailss=T)
