# -*- coding: utf-8 -*-
"""
Created on Thu Nov 17 09:20:35 2016
Propose: Recursive Thresholding techiques for spectrum signal/noise decision,
and introduce image segmentation method to solve the time-frequency 2-D threshold
making problem
@author: kevin
Environment: Python 2.7
"""

import matplotlib.pyplot as plt
from skimage import data
try:
    from skimage import filters
except ImportError:
    from skimage import filter as filters
from skimage import exposure

#load data
#load level data
#current data 1710-1730MHz, step by 25kHz
import scipy.io as sio
import numpy as np
t=sio.loadmat('timestamp.mat')
time=t["dateStamp"]
l=sio.loadmat('level.mat')
level=l["dataLevel"].astype('float')
#camera = data.camera()
#val = filters.threshold_otsu(camera)

#select part of the data by freq, presume da as a 2-d matrix
#col,row
def select_part_by_freq(da,start,step,stop,low,high,col_or_row):
    if((low<start)or(high>stop)):
        return
    lowix=int((low-start)/step)
    highix=int((high-start)/step)
    if col_or_row=="col":
        return da[:,lowix:highix]
    else:
        return da[lowix:highix,:]

#get part of level data
startf=1710
stopf=1740
stepf=0.025

starttime=0
stoptime=50
level_part=select_part_by_freq(level,startf,stepf,stopf,1710,1730,"col")

#a sliding window that will walk through all elements of level,
#and each return the window-inside part
def sliding_window_otsu(level,win_rowsize,win_colsize,move_step):
    [rownum,colnum]=level.shape
    channelstatus=np.ndarray((rownum,colnum),float)*0
    cs_inside=np.ndarray((win_rowsize,win_colsize),float)*0
    if(win_rowsize>=rownum)and(win_colsize>=colnum)and(move_step>=min(rownum,colnum)):
        print("too large step")
        return -1
    win_ix=0
    #while(str_rowix<rownum)and(str_colix<colnum):
    for str_rowix in np.arange(0,rownum,move_step): #start,stop,step
        for str_colix in np.arange(0,colnum,move_step):
            win_inside=level[str_rowix:(str_rowix+win_rowsize),str_colix:(str_colix+win_rowsize)]
            #calculate threshold for each window
            threshold_inside=filters.threshold_otsu(win_inside)
            cs_inside=(win_inside>threshold_inside)*1
            channelstatus[str_rowix:(str_rowix+win_rowsize),str_colix:(str_colix+win_rowsize)]+=cs_inside;
            win_ix+=1
    #make a decision, under how many windows' support can an element regard 
    #as signal
    channelstatus=(channelstatus>3)*1
    print("create",win_ix," temporal windows") 
    return channelstatus              
    
#calculation
#cs->val
cs = sliding_window_otsu(level_part,20,20,5)
plt.imshow(cs,cmap='gray',interpolation='nearest')
plt.savefig('cs.png',dpi=300,facecolor='w',edgecolor='w')
plt.hist(level_part[0,:],bins='auto')

#test ROHT algorithm
#t_tim=50
#t_frq=100
#data_part=level_part[0:t_tim,0:t_frq]
#data_vec=np.reshape(data_part,(t_tim*t_frq))
#plt.hist(data_vec,bins='auto')
#d_mean=np.mean(data_vec)
#d_std=np.std(data_vec)
##z_alph=1.645 as 95% probability 
#z_alph=1.645
#d_cutpoint=d_mean+z_alph*d_std #/np.sqrt(t_tim*t_frq) 
#data_list=data_vec.tolist()
##substract the signal part with 95% confidence
#data_list_filter=filter(lambda x: x<d_cutpoint, data_list) 

#recursive design
#recursive threshold to subtract signal part for noise/signal decision making
def recursive_oneside_hypthesis_testing(data_vec,max_ix):
   z_alph=1.645 #95% confidence
   d_mean=list()
   d_std=list()
   d_cutpoint=list()
   for ix in np.arange(max_ix):
       #collect mean and std of each iterate
       d_mean.append(np.mean(data_vec))
       d_std.append(np.std(data_vec))
       d_cutpoint.append(d_mean[ix]+z_alph*d_std[ix])
       #subtract the signal part with 95% confidence
       #and re-reshape to array
       data_list=data_vec.tolist()
       data_list_filter=filter(lambda x: x<d_cutpoint[ix], data_list)
       data_vec=np.asarray(data_list_filter)
       print ix
       print d_cutpoint[ix]
       #if(ix>0): 
           #if((d_std[ix]-d_std[ix-1])<0.2):
               #return [d_std,d_mean,d_cutpoint,ix]
   return [d_std,d_mean,d_cutpoint,ix]

#distribution type(z_alph), window-size(uniform noise, double-peak),
#def sliding_window_ROHT(more_data,win_tsize,win_fsize,move_step):
    
######1-D curve
t_tim=400
t_frq=400
data_part=level_part[0:t_tim,0:t_frq]
level_splic_time=level_part[0,:]
plt.plot(level_splic_time)
level_splic_freq=level_part[:,10]
plt.plot(level_splic_freq)
plt.savefig('level_time_slice.png')

#parse date time string
timestr=time[0]
from datetime import datetime
timeix=datetime.strptime(timestr,'%d-%b-%Y %H:%M:%S') #like '15-Dec-2015 19:00:00' 
    
######2-D image  

#scale the loaded data into range [0,255]
from sklearn import preprocessing
min_max_scaler = preprocessing.MinMaxScaler()
#give an explicit range
min_max_scaler.feature_range=(0,255)
data_image = min_max_scaler.fit_transform(data_part)
#filter
n=20
l=100
im = filters.gaussian_filter(data_image, sigma=l / (4. * n))
data_signal = im > im.mean()
plt.imshow(data_image,cmap='gray', interpolation='nearest')
plt.savefig('level_image.png')

from skimage import measure
#Label all connected components
all_labels = measure.label(data_signal)
blobs_labels = measure.label(data_signal, background=0)
plt.imshow(blobs_labels)
#Label only foreground connected components:
blobs_labels = measure.label(blobs, background=0)    

#################
data_vec=np.reshape(data_part,(t_tim*t_frq))
plt.hist(data_vec,bins='auto')

d_mean=list()
d_std=list()
d_cutpoint=list()
[d_std,d_mean,d_cutpoint,ix]=recursive_oneside_hypthesis_testing(data_vec,20)
threshold_ostu=filters.threshold_otsu(data_part)

print("OSTU:",threshold_ostu) #can distguish more weak signal
print("ROHT:",d_cutpoint[ix])


#applying image segementation method in threshold-making
check = np.zeros((9, 9))
check[::2, 1::2] = 1
check[1::2, ::2] = 1
import matplotlib.pyplot as plt
plt.imshow(check, cmap='gray', interpolation='nearest') 

from skimage import io
import os
filename = os.path.join(skimage.data_dir, 'camera.png')
camera = io.imread(filename)



##doc of savefig
#matplotlib.pyplot.savefig = savefig(*args, **kwargs)
#    Save the current figure.
#    
#    Call signature::
#    
#      savefig(fname, dpi=None, facecolor='w', edgecolor='w',
#              orientation='portrait', papertype=None, format=None,
#              transparent=False, bbox_inches=None, pad_inches=0.1,
#              frameon=None)
#    
#    The output formats available depend on the backend being used.
#    
#    Arguments:
#    
#      *fname*:
#        A string containing a path to a filename, or a Python
#        file-like object, or possibly some backend-dependent object
#        such as :class:`~matplotlib.backends.backend_pdf.PdfPages`.
#    
#        If *format* is *None* and *fname* is a string, the output
#        format is deduced from the extension of the filename. If
#        the filename has no extension, the value of the rc parameter
#        ``savefig.format`` is used.
#    
#        If *fname* is not a string, remember to specify *format* to
#        ensure that the correct backend is used.
#    
#    Keyword arguments:
#    
#      *dpi*: [ *None* | ``scalar > 0`` | 'figure']
#        The resolution in dots per inch.  If *None* it will default to
#        the value ``savefig.dpi`` in the matplotlibrc file. If 'figure'
#        it will set the dpi to be the value of the figure.
#    
#      *facecolor*, *edgecolor*:
#        the colors of the figure rectangle
#    
#      *orientation*: [ 'landscape' | 'portrait' ]
#        not supported on all backends; currently only on postscript output
#    
#      *papertype*:
#        One of 'letter', 'legal', 'executive', 'ledger', 'a0' through
#        'a10', 'b0' through 'b10'. Only supported for postscript
#        output.
#    
#      *format*:
#        One of the file extensions supported by the active
#        backend.  Most backends support png, pdf, ps, eps and svg.
#    
#      *transparent*:
#        If *True*, the axes patches will all be transparent; the
#        figure patch will also be transparent unless facecolor
#        and/or edgecolor are specified via kwargs.
#        This is useful, for example, for displaying
#        a plot on top of a colored background on a web page.  The
#        transparency of these patches will be restored to their
#        original values upon exit of this function.
#    
#      *frameon*:
#        If *True*, the figure patch will be colored, if *False*, the
#        figure background will be transparent.  If not provided, the
#        rcParam 'savefig.frameon' will be used.
#    
#      *bbox_inches*:
#        Bbox in inches. Only the given portion of the figure is
#        saved. If 'tight', try to figure out the tight bbox of
#        the figure.
#    
#      *pad_inches*:
#        Amount of padding around the figure when bbox_inches is
#        'tight'.
#    
#      *bbox_extra_artists*:
#        A list of extra artists that will be considered when the
#        tight bbox is calculated.
#fig=plt.imshow(cs,cmap='gray',interpolation='nearest')
#fig.show()
#fig.savefig('cs.png')

#val = filters.threshold_otsu(cs)
#hist, bins_center = exposure.histogram(cs)


#plot 
#plt.figure(figsize=(9, 4))
#plt.subplot(131)
#plt.imshow(cs, cmap='gray', interpolation='nearest')
#plt.axis('off')
#plt.subplot(132)
#plt.imshow(cs < val, cmap='gray', interpolation='nearest')
#plt.axis('off')
#plt.subplot(133)
#plt.plot(bins_center, hist, lw=2)
#plt.axvline(val, color='k', ls='--')
#
#plt.tight_layout()
#plt.show()