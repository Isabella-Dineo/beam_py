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
 
    num_H = np.random.randint(3,8)  # random number of discrete emission height 
#   np.random.seed(0)  # makes the random numbers predictable; the same numbers appear every time.   

#   Height for a short period pulsar:
    if P < 0.15:
        if hmin >= 950 and hmax <= 1000:
            H = np.random.uniform(hmin, hmax) # only one emission height 
           
        else: print "error, emission range not allowed for pulse period P < 0.15 seconds"
        
#   Height for a longer period pulsar:
    if P > 0.15:
        if hmin >= 20 and hmax <= 1000:
            H = np.random.uniform(hmin, hmax, size=num_H)
        
        else: print "error, emission range not allowed for pulse period P > 0.15 seconds"
        
    return H

########################### simple test #############################
if __name__ == "__main__":
    hmin = 955
    hmax = 975
    P = 0.015
    H = emission_height(P, hmin, hmax)
    print H    
