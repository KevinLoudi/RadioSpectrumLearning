# Author: Kevin
# Date: 7th March, 2017
# Propose: kriging spatial-temporal modeling
# Environment: R 3.3
####################################

#run script in shortcut
#library(gstat)
#demo(stkrige)
#demo(stkrige-prediction)
#demo(stkrige-crossvalidation)

#load large STSDF
library(gstat)
data(DE_RB_2005, package = "gstat")


dd<-DE_RB_2005@data
plot(dd$PM10)

paper <- FALSE

set.seed(123)
smplDays <- sort(sample(365,8))


# load German boundaries
data(air)
DE_NUTS1 <- spTransform(DE_NUTS1, CRS("+init=epsg:32632"))
if(!paper)
  plot(DE_NUTS1)

# station wise coverage
if(!paper)
  #plot days of each station
  barplot(sort(table(DE_RB_2005@index[,1])),
          main="reported days per station",
          ylab="number of days", xaxt="n")

# acf
if(!paper) {
    acf(DE_RB_2005[sample(68,1),,drop=F]@data)
    var(DE_RB_2005@data$PM10)
}
# DE_RB_2005[sample(68,1),,drop=F]@data #every 68 
# DE_RB_2005@data$PM10 #the spatial-temporal series

# a few daily snapshots
 if(paper)
     png("vignettes/figures/daily_means_PM10.png", width=9, height=6, "in", res=150)

 #DE_RB_2005[,smplDays]  #eight sample day
 
 stplot(as(DE_RB_2005[,smplDays],"STFDF"),
                col.regions=bpy.colors(120)[-(1:20)],
                 sp.layout = list("sp.polygons", DE_NUTS1), scales=list(draw=F), 
                 key.space="right", colorkey=T, cuts=0:70,
                 main="Eight sample day")
 
 if(paper)
      dev.off()
 
 # number of stations: 68 stations
  length(DE_RB_2005@sp)
  
  #calculate the empirical variogram: the first important step
  #formulation 
  empVgm <- variogramST(PM10~1, DE_RB_2005, tlags=0:6)
  DE_RB_2005@data #all available data
  (DE_RB_2005@sp@coords) #grids
  DE_RB_2005@time #time slot index
  
  #plot spatial-temporal empircal variogram
  if(!paper) {
       plot(empVgm, wireframe=T, scales=list(arrows=F))
       plot(empVgm)
     }

  # fit of theoretical purely spatial models #
  ############################################
  # take only variogram data of one step
  spEmpVgm <- empVgm[empVgm$timelag == 0,]
  
  class(spEmpVgm) <- c("gstatVariogram", "data.frame")
  
  spEmpVgm <- spEmpVgm[-1,1:3]
  
  spEmpVgm$dir.hor <- 0
  
  spEmpVgm$dir.ver <- 0
  
  #fit into theoretical models
  #vgm parameters: sill, model, range, nugget
  spVgmMod <- fit.variogram(spEmpVgm, vgm(80,"Exp",300000,20))
  print(spVgmMod,digit=3) #combined two model: nugget model + exponential model
  
  if(!paper)
      plot(spEmpVgm, spVgmMod)
  
  # fit of theoretical spatio-temporal models #
   #############################################
  #estimation the ansiotrpy
  linStAni <- estiStAni(empVgm, c(50000,200000))
  
  #demo variogram-distance relationship
  if(!paper) {
       plot(gamma~dist, empVgm[empVgm$timelag == 0,], ylim=c(0,100), xlim=c(0,800000))
       points(empVgm[empVgm$spacelag == 0,]$timelag*linStAni, empVgm[empVgm$spacelag == 0,]$gamma, col="red")
   }
  
  ##
  # rescale empVgm and linStAni to km for estimation
   empVgm$dist  <- empVgm$dist/1000
  
   empVgm$avgDist  <- empVgm$avgDist/1000
  
   empVgm$spacelag <- empVgm$spacelag/1000
  
   linStAni <- linStAni/1000
  
   # four categories of spatio-temporal model
   # separable
   separableModel <- vgmST("separable", 
                          space=vgm(0.9,"Exp", 200, 0.1),
                          time =vgm(0.9,"Sph", 3.5, 0.1),
                          sill=120)
   print(separableModel)
  
   #parameters: spatio-temporal sample variogram, fit model, fit method
   #  spatio-temporal anisotropy, 'L-BFGS-B' ?
   fitSepModel <- fit.StVariogram(empVgm, separableModel, fit.method = 7, 
                          stAni = linStAni, method = "L-BFGS-B", 
                          control = list(parscale=c(100,1,10,1,100)),
                          lower = c(10,0,.1,0,0.1), 
                          upper = c(2000,1,12,1,200))
  
   attr(fitSepModel, "optim.output")$value
   
   # Exp+Exp: 9.87, Exp+Sph: 6.82, Sph+Exp: 10.42, Sph+Sph: 7.50
   if(!paper)
        plot(empVgm, fitSepModel, wireframe=T, all=T, scales=list(arrows=F), zlim=c(0,135))
   
   # product-sum
    prodSumModel <- vgmST("productSum",
                                                  space=vgm(10, "Exp", 200, 1),
                                                  time= vgm(10, "Sph",   2, 1), 
                                                  k=2)
   
    fitProdSumModel <- fit.StVariogram(empVgm, prodSumModel, fit.method = 7, 
                                                  stAni = linStAni, method = "L-BFGS-B", 
                                                  control = list(parscale = c(1,10,1,1,0.1,1,10)),
                                                   lower = rep(0.0001,7))
    attr(fitProdSumModel, "optim.output")$value
    
    # Exp+Exp: 10.09, Exp+Sph: 6.91, Sph+Exp: 10.64, Sph+Sph: 7.59
    plot(empVgm, fitProdSumModel, wireframe=T, all=T, scales=list(arrows=F), zlim=c(0,135))
    
    # metric
    metricModel <- vgmST("metric",
                           joint = vgm(60, "Mat", 150, 10, kappa=0.6),
                           stAni = 60)
    
    fitMetricModel <- fit.StVariogram(empVgm, metricModel, fit.method = 7,
                           stAni = linStAni, method = "L-BFGS-B",
                           control = list(parscale = c(10,20,5,10)),
                           lower = c(80,50,5,50),
                           upper = c(200,1500,60,300))
    
    attr(fitMetricModel, "optim.output")$value
    
    # Exp: 10.25, Sph: 10.59,
    # Gau: 21.32, Mat 5: 18.20, Mat 2: 14.43, Mat 1.25: 12.04,
    # Mat 1: 11.07, Mat 0.75: 10.23, Mat 0.6: 10.05
    if(!paper)
       plot(empVgm, fitMetricModel, wireframe=T, all=T, scales=list(arrows=F), zlim=c(0,135))
    
    # simplified sumMetric model?
    sumMetricFromsimpleSumMetric <- function(vgm) {
      vgmST("sumMetric",
                space=vgm(vgm$space$psill, vgm$space$model, vgm$space$range, vgm$nugget/3),
                time =vgm(vgm$time$psill, vgm$time$model, vgm$time$range, vgm$nugget/3),
                joint=vgm(vgm$joint$psill, vgm$joint$model, vgm$joint$range, vgm$nugget/3),
                stAni=vgm$stAni)
      }
    
     simpleSumMetricModel <- vgmST("simpleSumMetric",
                space=vgm(120,"Sph", 150),
                time =vgm(120,"Exp", 10),
                joint=vgm(120,"Sph", 150),
                nugget=10, stAni=150)
    
    fitSimpleSumMetricModel <- fit.StVariogram(empVgm, simpleSumMetricModel,
                fit.method = 7, stAni=linStAni,
                method = "L-BFGS-B",
                lower = c(sill.s = 0, range.s = 10,
                                  sill.t = 0, range.t = 0.1,
                                  sill.st= 0, range.st= 10,
                                  nugget=0, anis = 40),
                          upper = c(sill.s = 200,  range.s = 500,
                                   sill.t = 200,  range.t = 20,
                                   sill.st= 200, range.st = 5000,
                                   nugget = 100, anis = 1000),
                          control = list(parscale = c(1,10,1,1,1,100,1,10)))
     
     
     attr(fitSimpleSumMetricModel, "optim.output")$value
     # Exp+Exp+Exp: 4.10 Exp+Sph+Exp: 3.60 Sph+Exp+Exp: 3.94 Sph+Sph+Exp: 3.32
     # Exp+Exp+Sph: 3.74 Exp+Sph+Sph: 3.98 Sph+Exp+Sph: 3.31 Sph+Sph+Sph: 3.56
       if(!paper)
           plot(empVgm,fitSimpleSumMetricModel, wireframe = T, scales = list(arrows = F), all = T , zlim=c(0,130))
     
     
     # sum-metric
     # sumMetricModel <- sumMetricFromsimpleSumMetric(fitSimpleSumMetricModel)

       sumMetricModel <- vgmST("sumMetric",
                   space = vgm(20, "Sph", 150, 1),
                    time = vgm(10, "Exp", 2, 0.5),
                    joint = vgm(80, "Sph", 1500, 2.5),
                    stAni = 120)
     
       fitSumMetricModel <- fit.StVariogram(empVgm, sumMetricModel, fit.method = 7, stAni=linStAni,
                    method = "L-BFGS-B", 
                    lower = c(sill.s = 0,  range.s = 10,  nugget.s = 0,
                                                      sill.t = 0,  range.t = 0.1,   nugget.t = 0,
                                                             sill.st= 0, range.st = 10, nugget.st = 0, 
                                                             anis = 40),
                    upper = c(sill.s = 200,  range.s = 1E3,  nugget.s = 20,
                                                      sill.t = 200,  range.t = 75,   nugget.t = 20,
                                                            sill.st= 200, range.st = 5E3, nugget.st = 20,
                                                             anis = 500),
                    control = list(parscale = c(1,100,1,1,0.5,1,1,100,1,100),
                                                      maxit=1e4))
       
       