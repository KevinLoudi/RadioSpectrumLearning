# -*- coding: utf-8 -*-
"""
Created on Sun Nov 06 09:15:18 2016
Propose: Fit power data into distribution
@author: kevin
Environment: Python 2.7
"""
#load level data
#current data 1710-1730MHz, step by 25kHz
import scipy.io as sio
import numpy as np
t=sio.loadmat('timestamp.mat')
time=t["dateStamp"]
l=sio.loadmat('level.mat')
level=l["dataLevel"].astype('float')

def select_part_by_freq(da,start,step,stop,low,high,col_or_row):
    if((low<start)or(high>stop)):
        return
    lowix=int((low-start)/step)
    highix=int((high-start)/step)
    if col_or_row=="col":
        return da[:,lowix:highix]
    else:
        return da[lowix:highix,:]

startf=1710
stopf=1740
stepf=0.025
level_part=select_part_by_freq(level,startf,stepf,stopf,1710,1730,"col")

#scale the loaded data into range [0,1]
from sklearn import preprocessing
min_max_scaler = preprocessing.MinMaxScaler()
#give an explicit range
min_max_scaler.feature_range=(0,1)
level_scaled_part = min_max_scaler.fit_transform(level_part)

#calculate channel state and show
import matplotlib.pyplot as plt
import pylab as pl
thre=0.35
cs=(level_scaled_part>thre)
im = plt.matshow(cs, cmap=pl.cm.hot, aspect='auto')
plt.colorbar(im)
plt.show()

#calculate time-occupy rate
def calc_occ(cs):
    [times,freqs]=cs.shape
    occ=np.arange(times)*0
    for i in np.arange(times):
      occ[i]=sum(cs[i,:])/freqs
    return occ
    
#occ=calc_occ(cs)
#plt.plot(occ)

#fit power data in find the best theoretical distribution model
import matplotlib.pyplot as plt
import scipy
import scipy.stats

dist_names = ['gamma', 'beta','rayleigh', 'norm', 'pareto']
#size of x and y expected to be indentical
#x=np.arange(1710,1730,0.025)
y=level_part[100,:]
size=y.size
x=np.arange(size)
#double-peed often indicate the existence of spectrum use
h = plt.hist(y, bins=range(size), color='w')
plt.show()

    #for dist_name in dist_names:
    #getattr(object, name[, default]) 
    #getattr(x, 'foobar') is equivalent to x.foobar
    #dist = getattr(scipy.stats, dist_name)
    #maximum likelihood estimation of distribution parameters
    #fit(self, data, *args, **kwds) 
    #Return MLEs for shape, location, and scale parameters from data.
dist_name=dist_names[2]
dist=getattr(scipy.stats,dist_name)
param = dist.fit(y)
pdf_fitted = dist.pdf(x, *param[:-2], loc=param[-2], scale=param[-1]) * size
plt.plot(pdf_fitted, label=dist_name)
plt.xlim(0,size-1)
print("Fit finished the ",dist_name)

plt.legend(loc='upper right')
plt.show()


#     |      Parameters
#     |      ----------
#     |      data : array_like
#     |          Data to use in calculating the MLEs.
#     |      args : floats, optional
#     |          Starting value(s) for any shape-characterizing arguments (those not
#     |          provided will be determined by a call to ``_fitstart(data)``).
#     |          No default value.
#     |      kwds : floats, optional
#     |          Starting values for the location and scale parameters; no default.
#     |          Special keyword arguments are recognized as holding certain
#     |          parameters fixed:
#     |      
#     |          - f0...fn : hold respective shape parameters fixed.
#     |            Alternatively, shape parameters to fix can be specified by name.
#     |            For example, if ``self.shapes == "a, b"``, ``fa``and ``fix_a``
#     |            are equivalent to ``f0``, and ``fb`` and ``fix_b`` are
#     |            equivalent to ``f1``.
#     |      
#     |          - floc : hold location parameter fixed to specified value.
#     |      
#     |          - fscale : hold scale parameter fixed to specified value.
#     |      
#     |          - optimizer : The optimizer to use.  The optimizer must take ``func``,
#     |            and starting position as the first two arguments,
#     |            plus ``args`` (for extra arguments to pass to the
#     |            function to be optimized) and ``disp=0`` to suppress
#     |            output as keyword arguments.
#     |      Returns
#     |      -------
#     |      shape, loc, scale : tuple of floats
#     |          MLEs for any shape statistics, followed by those for location and
#     |          scale.
#     |      
#     |      Notes
#     |      -----
#     |      This fit is computed by maximizing a log-likelihood function, with
#     |      penalty applied for samples outside of range of the distribution. The
#     |      returned answer is not guaranteed to be the globally optimal MLE, it
#     |      may only be locally optimal, or the optimization may fail altogether.
#     |      
#     |      
#     |      Examples
#     |      --------
#     |      
#     |      Generate some data to fit: draw random variates from the `beta`
#     |      distribution
#     |      
#     |      >>> from scipy.stats import beta
#     |      >>> a, b = 1., 2.
#     |      >>> x = beta.rvs(a, b, size=1000)
#     |      
#     |      Now we can fit all four parameters (``a``, ``b``, ``loc`` and ``scale``):
#     |      
#     |      >>> a1, b1, loc1, scale1 = beta.fit(x)
#     |      
#     |      We can also use some prior knowledge about the dataset: let's keep
#     |      ``loc`` and ``scale`` fixed:
#     |      
#     |      >>> a1, b1, loc1, scale1 = beta.fit(x, floc=0, fscale=1)
#     |      >>> loc1, scale1
#     |      (0, 1)
#     |      
#     |      We can also keep shape parameters fixed by using ``f``-keywords. To
#     |      keep the zero-th shape parameter ``a`` equal 1, use ``f0=1`` or,
#     |      equivalently, ``fa=1``:
#     |      
#     |      >>> a1, b1, loc1, scale1 = beta.fit(x, fa=1, floc=0, fscale=1)
#     |      >>> a1
#     |      1
#     |  
#     |  fit_loc_scale(self, data, *args)
#     |      Estimate loc and scale parameters from data using 1st and 2nd moments.
#     |      
#     |      Parameters
#     |      ----------
#     |      data : array_like
#     |          Data to fit.
#     |      arg1, arg2, arg3,... : array_like
#     |          The shape parameter(s) for the distribution (see docstring of the
#     |          instance object for more information).
#     |      
#     |      Returns
#     |      -------
#     |      Lhat : float
#     |          Estimated location parameter for the data.
#     |      Shat : float
#     |          Estimated scale parameter for the data.
#     |  
#     |  isf(self, q, *args, **kwds)
#     |      Inverse survival function (inverse of `sf`) at q of the given RV.
#     |      
#     |      Parameters
#     |      ----------
#     |      q : array_like
#     |          upper tail probability
#     |      arg1, arg2, arg3,... : array_like
#     |          The shape parameter(s) for the distribution (see docstring of the
#     |          instance object for more information)
#     |      loc : array_like, optional
#     |          location parameter (default=0)
#     |      scale : array_like, optional
#     |          scale parameter (default=1)
#     |      
#     |      Returns
#     |      -------
#     |      x : ndarray or scalar
#     |          Quantile corresponding to the upper tail probability q.
#     |  
#     |  logcdf(self, x, *args, **kwds)
#     |      Log of the cumulative distribution function at x of the given RV.
#     |      
#     |      Parameters
#     |      ----------
#     |      x : array_like
#     |          quantiles
#     |      arg1, arg2, arg3,... : array_like
#     |          The shape parameter(s) for the distribution (see docstring of the
#     |          instance object for more information)
#     |      loc : array_like, optional
#     |          location parameter (default=0)
#     |      scale : array_like, optional
#     |          scale parameter (default=1)
#     |      
#     |      Returns
#     |      -------
#     |      logcdf : array_like
#     |          Log of the cumulative distribution function evaluated at x
#     |  
#     |  logpdf(self, x, *args, **kwds)
#     |      Log of the probability density function at x of the given RV.
#     |      
#     |      This uses a more numerically accurate calculation if available.
#     |      
#     |      Parameters
#     |      ----------
#     |      x : array_like
#     |          quantiles
#     |      arg1, arg2, arg3,... : array_like
#     |          The shape parameter(s) for the distribution (see docstring of the
#     |          instance object for more information)
#     |      loc : array_like, optional
#     |          location parameter (default=0)
#     |      scale : array_like, optional
#     |          scale parameter (default=1)
#     |      
#     |      Returns
#     |      -------
#     |      logpdf : array_like
#     |          Log of the probability density function evaluated at x
#     |  
#     |  logsf(self, x, *args, **kwds)
#     |      Log of the survival function of the given RV.
#     |      
#     |      Returns the log of the "survival function," defined as (1 - `cdf`),
#     |      evaluated at `x`.
#     |      
#     |      Parameters
#     |      ----------
#     |      x : array_like
#     |          quantiles
#     |      arg1, arg2, arg3,... : array_like
#     |          The shape parameter(s) for the distribution (see docstring of the
#     |          instance object for more information)
#     |      loc : array_like, optional
#     |          location parameter (default=0)
#     |      scale : array_like, optional
#     |          scale parameter (default=1)
#     |      
#     |      Returns
#     |      -------
#     |      logsf : ndarray
#     |          Log of the survival function evaluated at `x`.
#     |  
#     |  nnlf(self, theta, x)
#     |      Return negative loglikelihood function.
#     |      
#     |      Notes
#     |      -----
#     |      This is ``-sum(log pdf(x, theta), axis=0)`` where `theta` are the
#     |      parameters (including loc and scale).
#     |  
#     |  pdf(self, x, *args, **kwds)
#     |      Probability density function at x of the given RV.
#     |      
#     |      Parameters
#     |      ----------
#     |      x : array_like
#     |          quantiles
#     |      arg1, arg2, arg3,... : array_like
#     |          The shape parameter(s) for the distribution (see docstring of the
#     |          instance object for more information)
#     |      loc : array_like, optional
#     |          location parameter (default=0)
#     |      scale : array_like, optional
#     |          scale parameter (default=1)
#     |      
#     |      Returns
#     |      -------
#     |      pdf : ndarray
#     |          Probability density function evaluated at x
#     |  
#     |  ppf(self, q, *args, **kwds)
#     |      Percent point function (inverse of `cdf`) at q of the given RV.
#     |      
#     |      Parameters
#     |      ----------
#     |      q : array_like
#     |          lower tail probability
#     |      arg1, arg2, arg3,... : array_like
#     |          The shape parameter(s) for the distribution (see docstring of the
#     |          instance object for more information)
#     |      loc : array_like, optional
#     |          location parameter (default=0)
#     |      scale : array_like, optional
#     |          scale parameter (default=1)
#     |      
#     |      Returns
#     |      -------
#     |      x : array_like
#     |          quantile corresponding to the lower tail probability q.
#     |  
#     |  sf(self, x, *args, **kwds)
#     |      Survival function (1 - `cdf`) at x of the given RV.
#     |      
#     |      Parameters
#     |      ----------
#     |      x : array_like
#     |          quantiles
#     |      arg1, arg2, arg3,... : array_like
#     |          The shape parameter(s) for the distribution (see docstring of the
#     |          instance object for more information)
#     |      loc : array_like, optional
#     |          location parameter (default=0)
#     |      scale : array_like, optional
#     |          scale parameter (default=1)
#     |      
#     |      Returns
#     |      -------
#     |      sf : array_like
#     |          Survival function evaluated at x
#     |  
