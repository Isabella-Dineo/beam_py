#!/usr/bin/env python 

import numpy as np
from scipy import constants
import emission_height

def rho(P, hmin, hmax, npatch):
    """Function to determine the opening angle rho given the rotational 
       period and emission height. Allowed emision heights for young pulsars 
       range between [950, 1000] and between [20, 1000] for old pulsars.
    
       Args:
       -----
       P      : rotational period (seconds)
       hmin   : minimum emission height (in km).
       hmax   : maximum emission height (in km).
       npatch : integer number of patches 
       
       Returns:
       --------
       rho    : the opening angle (degrees)

    """

    H = emission_height.emission_height(P, hmin, hmax, npatch)
    rho = np.degrees(np.sqrt((9 * np.pi * H * 10**3) / (2 * constants.c * P)))
    
    return rho

#################### simple test #####################################
if __name__ == "__main__":
    P = 0.16
    hmin = 30
    hmax = 990
    npatch = 3
    opa = rho(P, hmin, hmax, npatch)
    print opa
