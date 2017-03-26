#install.packages("fGarch",dependencies = TRUE)

require(fGarch)
## UNIVARIATE TIME SERIES INPUT:
# In the univariate case the lhs formula has not to be specified ...
# A numeric Vector from default GARCH(1,1) - fix the seed:
N = 200
x.vec = as.vector(garchSim(garchSpec(rseed = 1985), n = N)[,1])
garchFit(~ garch(1,1), data = x.vec, trace = FALSE)
# An univariate timeSeries object with dummy dates:
x.timeSeries = dummyDailySeries(matrix(x.vec), units = "GARCH11")
garchFit(~ garch(1,1), data = x.timeSeries, trace = FALSE)
## Not run:
# An univariate zoo object:
x.zoo = zoo(as.vector(x.vec), order.by = as.Date(rownames(x.timeSeries)))
garchFit(~ garch(1,1), data = x.zoo, trace = FALSE)
## End(Not run)
# An univariate "ts" object:
x.ts = as.ts(x.vec)
garchFit(~ garch(1,1), data = x.ts, trace = FALSE)
## MULTIVARIATE TIME SERIES INPUT: