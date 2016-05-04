#!/usr/bin/env python

import numpy as np
import emission_height

def patch_width(P, hmin, hmax, npatch):
    """Function to calculate the width of a patchy emission region 
       within a pulsar beam at a given height.
    
       Args:
       -----
       P             : rotational period (seconds).
       hmin          : minimum emission height (in km).
       hmax          : maximum emission height (in km).
       npatch        : integer number of emission patches.
       
       Returns:
       --------
       patchwidths   : the width of the patchy emission region (degrees).

    """   
    
    H = emission_height.emission_height(P, hmin, hmax, npatch)
    patchwidths = 2.45 * 0.2 * np.sqrt(H / ( 10 * P)) 
    
    return patchwidths

#################### simple test #################################
if __name__ == "__main__":
    P = 0.16
    hmin = 50
    hmax = 950
    npatch = 3
    patchwidths = patch_width(P, hmin, hmax, npatch)
    print patchwidths

