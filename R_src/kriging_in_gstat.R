require(gstat)

library(sp)
library(spacetime)
library(gstat)

#Constructs a spatio-temporal variogram with the given type and parameters
#set  spatio-temporal variogram model as "sumMetric"
#"stAni"-- A spatio-temporal anisotropy
sumMetricVgm <- vgmST("sumMetric",
                      space=vgm( 4.4, "Lin", 196.6,  3),
                      time =vgm( 2.2, "Lin",   1.1,  2),
                      joint=vgm(34.6, "Exp", 136.6, 12),
                      stAni=51.7)

data(air)

if (!exists("rural"))
  rural = STFDF(stations, dates, data.frame(PM10 = as.vector(air)))

rr <- rural[,"2005-06-01/2005-06-03"]
rr <- as(rr,"STSDF")

x1 <- seq(from=6,to=15,by=1)
x2 <- seq(from=48,to=55,by=1)

DE_gridded <- SpatialPoints(cbind(rep(x1,length(x2)), rep(x2,each=length(x1))), 
                            proj4string=CRS(proj4string(rr@sp)))
gridded(DE_gridded) <- TRUE
DE_pred <- STF(sp=as(DE_gridded,"SpatialPoints"), time=rr@time)

#krigeST using spatio-temporal covariance models "ordinary spatio-temporal (ST) kriging"
#kriging on point support Usage
#krigeST(formula, data, newdata, modelList, y, beta, nmax = Inf, stAni = NULL,
#        computeVar = FALSE,	fullCovariance = FALSE,
#        bufferNmax=2, progress=TRUE)

#formula defines the dependent variable as a linear model of independent variables
#ordinary and simple kriging
#'data'--- contain the dependent variable and independent variables
#'newdata'--- prediction/simulation locations in space and time
#'modelList'--- model type
DE_kriged <- krigeST(PM10~1, data=rr, newdata=DE_pred,
                     modelList=sumMetricVgm)
gridded(DE_kriged@sp) <- TRUE
stplot(DE_kriged)