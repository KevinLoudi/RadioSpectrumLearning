# -*- coding: utf-8 -*-
"""
Created on Thu Nov 17 09:20:35 2016
Propose: Recursive Thresholding techiques for spectrum signal/noise decision
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