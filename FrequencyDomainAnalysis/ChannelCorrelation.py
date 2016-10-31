# -*- coding: utf-8 -*-
"""
Created on Thu Oct 20 09:44:40 2016
Propose: Analysis correlation between two channel
@author: kevin,James
"""
"""
Macro's point is correct the proper way to compare 
for relationships between time series is by the cross-correlation
function (assuming stationarity). Having the same length is not 
essential. The cross correlation at lag 0 just computes a correlation 
like doing the Pearson correlation estimate pairing the data at the 
identical time points. If they do have the same length as you are assuming, 
you will have exact T pairs where T is the number of time points for each 
series. Lag 1 cross correlation matches time t from series 1 with time t+1 
in series 2. Note that here even though the series are the same length you 
only have T-2 pair as one point in the first series has no match in the second 
and one other point in the second series will not have a match from the first. 
Given these two series you can estimate the cross-correlation at several lags . 
If any of the cross correlations is statistically significantly different from 
0 it will indicate a correlation between the two series.
"""

from sklearn import linear_model

clf = linear_model.Lasso(alpha=0.1)
clf.fit([[0,0], [1, 1], [2, 2]], [0, 1, 2])
print(clf.coef_)
print(clf.intercept_)