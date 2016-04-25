#!/usr/bin/env python

import numpy as np

def emission_height(P, hmin, hmax):
    """Function to determine the emission height range given hmin and hmax.
    
       Args:
       -----
       P    : rotational period.
       hmin : minimum emission height (in km).
       hmax : maximum emission height (in km).
       
       Returns:
       --------
       H    : random emission height (km).

    """
 
    num_H = np.random.randint(3,10)  # random number of discrete emission height 
#   np.random.seed(0)  # makes the random numbers predictable; the same numbers appear every time.   

#   Height for a short period pulsar:
    if P < 0.15:
        H = np.random.uniform(hmin, hmax) # only one distinct emission height 
    
        
#   Height for a longer period pulsar:
    if P > 0.15:
        H = np.random.uniform(hmin, hmax, size=num_H)        
        
    return H
