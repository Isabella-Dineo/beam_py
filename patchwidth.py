#!/usr/bin/env python

import numpy as np
import emission_height

def patch_width(P, hmin, hmax):
    """Function to calculate the width of a patchy emission region 
       within a pulsar beam at a given height.
    
       Args:
       -----
       P    : rotational period (seconds).
       hmin : minimum emission height (in km).
       hmax : maximum emission height (in km).
       
       Returns:
       --------
       wp   : the width of the patchy emission region (degrees).

    """   
    
    H = emission_height(P, hmin, hmax)
    wp = 2.45 * 0.2 * np.sqrt(H / ( 10 * P)) 
    
    return wp

#################### simple test #########################
if __name__ == "__main__":
    hmin = 955
    hmax = 975
    P = 0.015
    wp = patch_width(P, hmin, hmax)
