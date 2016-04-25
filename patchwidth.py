#!/usr/bin/env python

import numpy as np
import emission_height

def patcheswidth(P, hmin, hmax):
    """Function to calculate the width of a patchy emission region within a radio pulsar beam at a given the
       rotational period of the pulsar.
    
       Args:
       -----
       P    : rotational period (seconds).
       
       Returns:
       --------
       wp   : the width of the patchy emission region (degrees).

    """   
    
    H = emission_height(P, hmin, hmax)
    wp = 2.45 * 0.2 * np.sqrt(H / ( 10 * P)) 
    
    return wp
