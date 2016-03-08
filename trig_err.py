#!/usr/bin/python

import numpy as np

# Functions to compute trig calculations in radians 

def cosD(ang):
    """Returns the cosine of the input ang in radians
    ang: float, in degrees.
    """
    angR = np.deg2rad(ang)
    return np.cos(angR)

def sinD(ang):
    """Returns the sine of the input ang in radians
    ang: float, in degrees.
    """
    return np.sin(np.deg2rad(ang))

def tanD(ang):
    """Returns the tan of the input ang in radians
    ang: float, in degrees.
    """
    return np.tan(np.deg2rad(ang))

def acosD(ang):
    """Returns the arccos of the input ang in radians
    ang: float, in degrees.
    """
    return np.arccos(np.deg2rad(ang))

def asinD(ang):
    """Returns the arcsin of the input ang in radians
    ang: float, in degrees.
    """
    return np.arcsin(np.deg2rad(ang))

def atanD(ang):
    """Returns the arctan of the input ang in radians
    ang: float, in degrees.
    """
    return np.arctan(np.deg2rad(ang))

#########################################################

# correct for round-off errors 

def correct(x):
    """Function to correct for round-off errors to 7 
       decimal places.
    """
    tol = 1.0e-07
    sign = np.sign(x)
    y = np.abs(x)

    if np.abs(y-1.0) < tol:
        y = 1.0
    elif np.abs(y-0.0) < tol:  
          y = 0.0

    return y * sign

##########################################################

# simple test statement

if __name__ == "__main__":
	
    print cosD(90.)
    print sinD(90.)
    print tanD(90.)
    print asinD(0.5)
    print acosD(0.5)
    print atanD(0.5)
    print correct(1.0e-7)

