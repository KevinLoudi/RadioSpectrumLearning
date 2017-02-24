# Author: Kevin
# Date: 24th Feb, 2017
# Propose: test methods in gstat 
# Environment: R 3.3
####################################
library(sp)
library(spacetime)
data(coalash)
summary(coalash)

data("DE_RB_2005")
str(DE_RB_2005)

#estimate spatio-temporal anisotropy
data(vv)
estiStAni(vv, c(10, 150))
estiStAni(vv, c(10, 150), "vgm", vgm(80, "Sph", 120, 20))

#set up a variogram models
sumMetricModel <- vgmST("sumMetric",
                        space=vgm(30, "Sph", 200, 6),
                        time =vgm(30, "Sph", 15, 7),
                        joint=vgm(60, "Exp", 84, 22),
                        stAni=100)
extractPar(sumMetricModel)
extractParNames(sumMetricModel)

#fit sample data to a given type of 
#spatio-temporal variogram model
# separable model: spatial and temporal sill will be ignored
# and kept constant at 1-nugget respectively. A joint sill is used.
## Not run:
separableModel <- vgmST("separable",
                        method = "Nelder-Mead", # no lower & upper needed
                        space=vgm(0.9,"Exp", 123, 0.1),
                        time =vgm(0.9,"Exp", 2.9, 0.1),
                        sill=100)
data(vv)
separableModel <- fit.StVariogram(vv, separableModel,
                                  method="L-BFGS-B",
                                  lower=c(10,0,0.01,0,1),
                                  upper=c(500,1,20,1,200))
plot(vv, separableModel)
                                  
                    
#spatial-only modeling and prediction 
library(sp)
library(spacetime)
library(gstat)
data(meuse)
coordinates(meuse) = ~x+y
data(meuse.grid)
gridded(meuse.grid) = ~x+y
m <- vgm(.59, "Sph", 874, .04)
# ordinary kriging:
x <- krige(log(zinc)~1, meuse, meuse.grid, model = m)
spplot(x["var1.pred"], main = "ordinary kriging predictions")
spplot(x["var1.var"], main = "ordinary kriging variance")
# simple kriging:
x <- krige(log(zinc)~1, meuse, meuse.grid, model = m, beta = 5.9)
# residual variogram:
m <- vgm(.4, "Sph", 954, .06)
# universal block kriging:
x <- krige(log(zinc)~x+y, meuse, meuse.grid, model = m, block = c(40,40))
spplot(x["var1.pred"], main = "universal kriging predictions")
# krige0, using user-defined covariance function and multiple responses in y:
# exponential variogram with range 500, defined as covariance function:
v = function(x, y = x) { exp(-spDists(coordinates(x),coordinates(y))/500) }
# krige two variables in a single pass (using 1 covariance model):
y = cbind(meuse$zinc,meuse$copper,meuse$lead,meuse$cadmium)
x <- krige0(zinc~1, meuse, meuse.grid, v, y = y)
meuse.grid$zinc = x[,1]
spplot(meuse.grid["zinc"], main = "zinc")
meuse.grid$copper = x[,2]
spplot(meuse.grid["copper"], main = "copper")
# the following has NOTHING to do with kriging, but --
# return the median of the nearest 11 observations:
x = krige(zinc~1, meuse, meuse.grid, set = list(method = "med"), nmax = 11)
# get 25%- and 75%-percentiles of nearest 11 obs, as prediction and variance:
x = krige(zinc~1, meuse, meuse.grid, nmax = 11,
          set = list(method = "med", quantile = 0.25))
# get diversity (# of different values) and mode from 11 nearest observations:
x = krige(zinc~1, meuse, meuse.grid, nmax = 11, set = list(method = "div"))

library(sp)
data(meuse)
coordinates(meuse) <- ~x+y
m <- vgm(.59, "Sph", 874, .04)
# five-fold cross validation:
x <- krige.cv(log(zinc)~1, meuse, m, nmax = 40, nfold=5)
bubble(x, "residual", main = "log(zinc): 5-fold CV residuals")
# multivariable; thanks to M. Rufino:
meuse.g <- gstat(id = "zn", formula = log(zinc) ~ 1, data = meuse)
meuse.g <- gstat(meuse.g, "cu", log(copper) ~ 1, meuse)
meuse.g <- gstat(meuse.g, model = vgm(1, "Sph", 900, 1), fill.all = TRUE)
x <- variogram(meuse.g, cutoff = 1000)
meuse.fit = fit.lmc(x, meuse.g)
out = gstat.cv(meuse.fit, nmax = 40, nfold = 5)
summary(out)
out = gstat.cv(meuse.fit, nmax = 40, nfold = c(rep(1,100), rep(2,55)))
summary(out)
# mean error, ideally 0:
mean(out$residual)
# MSPE, ideally small
mean(out$residual^2)
# Mean square normalized error, ideally close to 1
mean(out$zscore^2)
# correlation observed and predicted, ideally 1
cor(out$observed, out$observed - out$residual)
# correlation predicted and residual, ideally 0
cor(out$observed - out$residual, out$residual)