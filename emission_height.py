import numpy as np 

def emission_height(P, hmin, hmax, npatch):
    """Function to determine the emission height range given hmin and hmax. Allowed emision heights for young 
       pulsars range between [950, 1000] and between [20, 1000] for old pulsars.
    
       Args:
       -----
       P      : rotational period.
       hmin   : minimum emission height (in km).
       hmax   : maximum emission height (in km).
       npatch : integer number of emission patches.
       
       Returns:
       --------
       H      : random emission height.
    """
    
    #num_H = np.random.randint(3,8)  # random number of discrete emission height 
    num_H = npatch
    #np.random.seed(0)  # makes the random numbers predictable; the same numbers appear every time.
    
#   emission height for a short period pulsar: only one emission height 
    if P <= 0.15:
        if hmin >= 950 and hmax <= 1000:
            H = np.random.uniform(hmin, hmax)
            
        else: print "error, emission range not allowed for pulse period P < 0.15 seconds"
            
#   emission height for a long period pulsar:        
    if P > 0.15:
        if hmin >= 20 and hmax <= 1000:
            H = np.random.uniform(hmin, hmax, size=num_H)
        
        else: print "error, emission range not allowed for pulse period P > 0.15 seconds"
        
    return H


########################### simple test #############################
if __name__ == "__main__":
    P = 0.16
    hmin = 50
    hmax = 950
    npatch = 3
    H = emission_height(P, hmin, hmax, npatch)
    print H  
