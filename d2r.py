#!/usr/bin/python 
import numpy as np

# Functions to compute trig calculations in radians 

def cosD(ang):
    """
    Calculate the cosine of the input ang in radians

    Args:
    -----
       ang: angle in degrees.
    """
    angR = np.deg2rad(ang)
    return np.cos(angR)

def sinD(ang):
    """
    Calculate the sine of the input ang in radians

    Args:
    -----
       ang: angle in degrees.
    """
    return np.sin(np.deg2rad(ang))

def tanD(ang):
    """
    Calculate the tan of the input ang in radians

    Args:
    -----
       ang: angle in degrees.
    """
    return np.tan(np.deg2rad(ang))

def acosD(ang):
    """
    Calculate the inverse cosine of the input ang in radians

    Args:
    -----
       ang: angle in degrees.
    """
    return np.arccos(np.deg2rad(ang))

def asinD(ang):
    """
    Calculate the inverse sine of the input ang in radians

    Args:
    -----
       ang: angle in degrees.
    """
    return np.arcsin(np.deg2rad(ang))

def atanD(ang):
    """
    Calculate the inverse tan function of the input ang in radians

    Args:
    -----
       ang: angle in degrees.
    """
    return np.arctan(np.deg2rad(ang))

#####################################################################

# correct for round-off errors 

def correct(x):
    """Function to correct for round-off errors to 7 
       decimal places; the input x should be an ndarray.
    """
    tol = 1.0e-07
    sign = np.sign(x)
    y = np.abs(x)

    for i in np.arange(np.size(x)):
        if np.abs(y[i] - 1.0) < tol:
            y[i] = 1.0
        elif np.abs(y[i] - 0.0) < tol:
            y[i] = 0.0

    return y * sign

####################################################################

# simple test statement

if __name__ == "__main__":
	
    print cosD(90.)
    print sinD(90.)
    print tanD(90.)
    print asinD(0.5)
    print acosD(0.5)
    print atanD(0.5)
    print correct(1.0e-7)


