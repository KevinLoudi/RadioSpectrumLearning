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
u <- runif(200,0,1)

#translate into time series
ix=as.POSIXct(ti_ix$tot.time[dc$dc!=0])
dc=ts(dc$dc[dc$dc!=0])

i<-1:162 #every hour a index
i<-i*12  #max index: 1944

ixx<-ix[i]
dcc<-dc[i]
plot.ts(dcc,type='h')

Box.test(dcc,lag = 20,type="Ljung-Box")
Box.test(u,lag = 20,type="Ljung-Box")

library("fUnitRoots")
adfTest(dcc,lags=10,type = "nc",title=NULL,description = NULL)
adfTest(u,lags=10,type = "nc",title=NULL,description = NULL)

#differ to reach stationary
dcc.dif<-diff(dcc,differences=1)
plot.ts(dcc.dif)
adfTest(dcc.dif,lags=10,type = "nc",title=NULL,description = NULL)

#lag p exceeds the significance bounds p=2
graphics.off()
par(mar=c(3,3,3,3))
attach(mtcars)
par(mfrow=c(1,2))
acf(dcc,lag.max=20,main="ACF",xlab="Lags",ylab="corr.")
#slowly decreasing in magnitude with increasing lag q=10
pacf(dcc,lag.max=20,main="PACF",xlab="Lags",ylab="corr.")


#AIC and BIC (1,1,2)
require(forecast)
model_par<-forecast::auto.arima(dcc,ic="bic")
model_par
fit <- Arima(dcc, order=c(3,0,1), seasonal=FALSE)
plot_forecast_error(fit$residuals)
plot(fit$x, col='red',main="ARIMA(3,0,1)拟合结果[实线]与SCR样本[空心点]",
     xlab="监测周期数",ylab="SCR指标")
lines(fitted(fit), col='blue')
fit$sigma2
SCR=dcc
dcc.for<-sarima.for(SCR,24,3,0,1)
plot(dcc)
lines(dcc.for$pred,col='red')


fit <- Arima(dcc, order=c(3,0,1), seasonal=FALSE)
plot_forecast_error(fit$residuals)
plot(fit$x, col='red',main="ARIMA(3,0,1) fitted result and SCR sample",
     xlab="Time",ylab="SCR")
lines(fitted(fit), col='blue')
fit$sigma2


#
require(astsa)
graphics.off()
dcc.sarima<-sarima(dcc, 3,0,1,1,0,1,24)
dcc.sarima
#using standard model
mod2<-Arima(dcc,order=c(3, 0, 1),
            seasonal=list(order=c(1, 0, 1), period=24))
plot(mod2$x, col='red',main="SARIMA(3,0,1)x(1,0,1)[24] fitted result and SCR sample",
     xlab="Time",ylab="SCR")
lines(fitted(mod2), col='blue')
SCR=dcc
sarima.for(SCR, 24, 3,0,1,1,0,1,24)
predict(mod2, n.ahead=24)

#fit error analysis
plot_forecast_error(dc.arima$residuals)
plot_forecast_error(dc.forecast$residuals)

#SARIMA
#require(astsa)
#dc.sarima<-sarima(dc,1,1,2,0,0,0,177,details = TRUE,Model=TRUE)

fit <- Arima(dc, order=c(2,1,2), seasonal=c(1,1,1))
tsdisplay(residuals(fit))
auto.arima(dc,stepwise=FALSE, approximation=FALSE)

#check if ARCH effect exist
require(FinTS)
require(rugarch)
ArchTest(dc.arima$residuals)
arch.spec=ugarchspec(variance.model=list(garchOrder=c(7,0)),
                     mean.model=list(armaOrder=c(0,0)))


############1 Autoregressive models
#In a multiple regression model, we forecast the 
#variable of interest using a linear combination of 
#predictors. In an autoregression model, we forecast 
#the variable of interest using a linear combination 
#of past values of the variable.

#The term autoregression indicates that it is a regression 
#of the variable against itself.
#AR handling a wide range of different time series patterns.

############2 Moving average models
#a moving average model uses past forecast errors in a 
#regression-like model.
#A moving average model is used for forecasting future 
#values while moving average smoothing is used for 
#estimating the trend-cycle of past values.


############3 Non-seasonal Autoregressive Moving average models
# The predictors act on the right hand side include both lagged values
# and lagged errors
arima.m<-arima.sim(list(order = c(0,0,12), ma = c(0.7,rep(0,10),0.9)), n = 200)
acf(arima.m)
plot(arima.m, type='h')

ar_pacf<-ARMAacf (ar = c(.6,0,0,0,0,0,0,0,0,0,0,.5,-.30),lag.max=30,pacf=T)
plot(ar_pacf, type='h')

x<-ts(dc[dc!=0])

graphics.off()
mod1<-sarima(x, 1,1,2,0,0,0,288)
x.forcast<-forecast.Arima(mod1,h=12,level=c(99.5))

ix<-1:162
ix<-ix*12
xx<-x[ix]
plot(xx,type='h')
#try a model
require(astsa)
mod1<-sarima(xx, 1,0,1,2,1,0,24)
mod1<-sarima(x, 1,1,2,1,1,0,24*12)

#using standard model
mod2<-Arima(xx,order=c(1, 0, 1),
            seasonal=list(order=c(2, 1, 0), period=24))
mod2
mod2<-Arima(xx,order=c(1, 0, 1),
            seasonal=list(order=c(2, 1, 0), period=24*12))
mod2

#see fit level
plot(mod2$x, col='red')
lines(fitted(mod2), col='blue')

#do forcast 
sarima.for(xx, 24, 1,0,1,2,1,0,24)

#native commend
predict(mod2, n.ahead=24)
