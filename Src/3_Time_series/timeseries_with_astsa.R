# Author: Kevin
# Date: 6th March, 2017
# Propose: time series with astsa package
# Environment: R 3.3
####################################

library(astsa)             # then load it (has to be done at the start of each session)
dc <- read.csv(file="D:/Code/WorkSpace/ThesisCode/Src/3_Time_series/TimeSeriesDataset/dc_2_1740_1750.csv", header=FALSE)

#access methods for observed data
dc[1]
length(dc)
dim(dc)
nrow(dc)
ncol(dc)

dcm=as.matrix(dc)
dim(dcm)

plot(dc$V1, ylab="duty cycle", main="variation of duty cycle",type="o", col="blue", lty="dashed")
