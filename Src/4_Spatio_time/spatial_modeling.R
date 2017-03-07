# Author: Kevin
# Date: 7th March, 2017
# Propose: kriging spatial modeling (pure spatial)
# Environment: R 3.3
####################################

#####using built-in 'meuse' dataset
library(sp)
data(meuse)

#assign grid-net and show details of data set
names(meuse) #variables in the dataset
coordinates(meuse)=~x+y #form a grid-net ('matrix') and promotes 'data.frame' into
  #spatialpoint.dataframe

coordinates(meuse)[1:5,] #show part of the grid-net
class(meuse)
summary(meuse) #summarize the dataset

#regular grids
data(meuse.grid)
summary(meuse.grid)
class(meuse.grid)
coordinates(meuse.grid) = ~x+y
class(meuse.grid)
gridded(meuse.grid) = TRUE
image(meuse.grid["dist"]) #demostrate the grid
title("distance to river (red = 0)")


#bubble ploting, x- y- are coordinates
bubble(meuse,"zinc",col=c("#00ff0088", "#00ff0088"),
       main = "zinc concentrations (ppm)")

library(gstat)
#inverse distance weighted methods
zinc.idw = idw(zinc~1, meuse, meuse.grid)
#larger concentrations are measured at locations 
#  close to the river
spplot(zinc.idw["var1.pred"], main = "zinc inverse distance weighted interpolations")

#linearize the concentrations with zinc log-transforming
# considering the distance to the river
plot(log(zinc)~sqrt(dist), meuse)
abline(lm(log(zinc)~sqrt(dist), meuse))

#calculate variograms
lzn.vgm = variogram(log(zinc)~1, meuse)
lzn.vgm
#fit emprical variogram to a variogram model
#'sph' model
lzn.fit = fit.variogram(lzn.vgm, model = vgm(1, "Sph", 900, 1))
plot(lzn.vgm, lzn.fit) #compare the emprical variogram and fitted variogram model
#'exp'
lznr.vgm = variogram(log(zinc)~sqrt(dist), meuse)
lznr.fit = fit.variogram(lznr.vgm, model = vgm(1, "Exp", 300, 1))
plot(lznr.vgm, lznr.fit)

#kriging prediction (ordinary kriging)
lzn.kriged = krige(log(zinc)~1, meuse, meuse.grid, model = lzn.fit)
spplot(lzn.kriged["var1.pred"])

#############
#calculate directional sample variogram
# in four direction angle
lzn.dir = variogram(log(zinc)~1, meuse, alpha = c(0, 45, 90, 135))
lzndir.fit = vgm(.59, "Sph", 1200, .05, anis = c(45, .4))
plot(lzn.dir, lzndir.fit, as.table = TRUE)

#using sqrt distance
lznr.dir = variogram(log(zinc)~sqrt(dist), meuse, alpha = c(0, 45, 90, 135))
plot(lznr.dir, lznr.fit, as.table = TRUE)

#draw a variogram maps
vgm.map = variogram(log(zinc)~sqrt(dist), meuse, cutoff = 1500, width = 100,
                    map = TRUE)
#show only semivariogram map values based on at least 5 point pairs
plot(vgm.map, threshold = 5)
