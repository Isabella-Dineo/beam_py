#!/usr/bin/env python

import numpy as np
import mapphi as mp

def los(alpha, beta):
    """Function to determine the line of sight cut across the beam.
    
       Args:
       -----
       alpha    : inclination angle
       beta     : impact parameter
       
       Returns:
       --------
       xlos     : the line of sight x-coordinates
       ylos     : the line of sight y-coordinates
       thetalos : the line of sight angle in degrees
    """
    
    phi = np.arange(-180,180)
    xp, yp = mp.mapphi(alpha, beta, phi)
    thetalos = np.arctan2(yp, xp) * (180 / np.pi) - 90.0
    
    for i in np.arange(len(thetalos)):
        if thetalos[i] < 0:
            thetalos[i] = -thetalos[i]
    
    return xp, yp, thetalos

##################### simple test #########################
if __name__ == "__main__":
    alpha = 90
    beta = 10
    xp, yp, thetalos = los(alpha, beta)
    print thetalos
