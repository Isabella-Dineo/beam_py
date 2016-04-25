#!/usr/bin/env python 

import numpy as np
from scipy import constants
import emission_height

def rho(P, hmin, hmax):
    """Function to determine the opening angle rho given the rotational 
       period and emission height. Allowed emision heights for young pulsars 
       range between [950, 1000] and between [20, 1000] for old pulsars.
    
       Args:
       -----
       P   : rotational period (seconds)
       hmin : minimum emission height (in km).
       hmax : maximum emission height (in km).
       
       Returns:
       --------
       rho : the opening angle (degrees)

    """

    H = emission_height(P, hmin, hmax)
    rho = np.degrees(np.sqrt((9 * np.pi * H) / (2 * constants.c * P)))
    
    return rho

#################### simple test #####################################
if __name__ == "__main__":
    hmin = 955
    hmax = 975
    P = 0.015
    opa = rho(P, hmin, hmax)
    print opa
