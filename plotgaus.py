import numpy as np
import matplotlib.pyplot as plt
import d2r
import patch_width
import los
import patch_center


def plotpatch(P, alpha, beta, hmin, hmax, npatch):
    """Function to plot the patches for a given height range. Using a 2d gaussian
    
       Args:
       -----
       P       : rotational period (seconds)
       alpha   : inclination angle (degrees)
       beta    : impact parameter (degrees)
       hmin    : minimum emission height (in km)
       hmax    : maximum emission height (in km)
       npatch  : number of emission patches
       
       Returns:
       --------
       A plot of the patches projected on to observational plane.
    
    """
    
    X = np.linspace(-180, 180, num=50, endpoint=True)
    Y = np.linspace(-180, 180, num=50, endpoint=True)
    patchwidths = patch_width.patch_width(P, hmin, hmax)
    
#   choose random patch widths (wp) depending on number of patches specified:
    np.random.shuffle(patchwidths)
    wp = patchwidths[0:npatch] 

#   wp == the spread (sigma)
#   project the sigmax and sigmay to the line of sight plane:
    xlos, ylos, thetalos = los.los(alpha, beta)
    sigmax = wp * d2r.sinD(thetalos)
    sigmay = wp * d2r.cosD(thetalos)
    
#   patchcenter or mean of the plot from patch center:
    cx, cy = patch_center.patch_center(P, hmin, hmax, npatch)
    
#   2d gaussian function
#    for i in np.arange(len(cx)):
#        gauss_2D = peak*np.exp(-((X - cx[i])**2 / 2*sigmax**2 + (Y - cy[i])**2 / 2*sigmay**2)) 
    peak = 5 # trial peak 
    for i in np.arange(len(cx)):
        gauss_2D = peak*np.exp(-((X - cx[i]) ** 2 / 2 * sigmax[i] ** 2 + (Y - cy[i]) ** 2 / (2 * sigmay[i] **2)))
 
    
    #plt.imshow(gauss_2D)
    plt.grid()
    plt.show()
