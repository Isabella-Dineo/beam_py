#!/usr/bin/python

import numpy as np

def deriv(inar):
    """
    Function to calculate the central difference
    as an approximation of the derivative.

    Args:
    -----

        inar     : the input array

    Returns:
    --------

        outar    : the output array 

    outar[1] and outar[-1] will use forward and backward difference
    respectively.
    """

    outar = np.zeros(inar.shape[0], dtype=float) #assume 1-D input array
    outar[2:-2] = 0.25*(inar[4:] - inar[:-4])
    outar[1] = 0.5*(inar[2] - inar[0])
    outar[0] = outar[1]
    outar[-2] = 0.5*(inar[-1] - inar[-3])
    outar[-1] = inar[-1] - inar[-2]
    
    return outar
    
#########################################################################
#simple test statement

if __name__ == "__main__":
    #inar = np.random.random(100)
    inar = np.arange(100)
    outar = deriv(inar)
    print np.allclose(outar, 1.)
