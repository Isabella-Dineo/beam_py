#!/usr/bin/python
import numpy as np

def cendif(inar,jo,lag):
    """
    Function to calculate the central difference
    as an approximation of the derivative
     
    Args:
    -----

        inar     : the input array
        dim      : the dimension of the two arrays (integer)
        jo       : takes in the desired threshold 
        lag      : writes out the first lag 
     
    Returns:
    --------

        outar    : the output array

    outar[1] and outar[dim] will use forward and backward difference
    respectively.
    """

    dim = inar.shape[0] #assuming 1-D array
    outar = np.zeros(dim)
    outar[1:-1] = .5*(inar[:-2] - inar[2:])
    outar[0] = outar[1]
    outar[-1] = outar[-2]

    for i in np.arange(1, dim-2):
        if outar[i] > jo and outar[i + 1] > jo:
           lag = i
           continue

    return outar, lag

#########################################################################
#simple test statement

if __name__ == "__main__":
    inar = np.array([1.,2.,3.,10.,50.,6.])
    print cendif(inar,0.5,0.5)

