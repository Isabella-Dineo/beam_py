#!/usr/bin/python 

# Code to map the rotational phase a pulsar, given the geometry of the pulsar. 
# Input angles are the inclination angle, the impact parameter, and the
# polarisation position angle phi. The function takes input angle in degrees and
# return an angle R in radians, which in turn used to calculate the x and y 
# coordinates of the plane of rotation 

import numpy as np
import d2r

def mapphi(alpha,beta,phi):
    """Maps the rotational phase (pulse longitude) phi.
     
       Args: 
         alpha   : inclination angle of the magnetic axis from the rotational axis. 
         beta    : impact parameter (the closest approach of the magnetic axis to the LOS.)
         phi     : polarisation position angle.

       Return:
         The coordinates of the plane of rotation xp and yp
    """
    cosR = d2r.cosD(alpha+beta) * d2r.cosD(alpha) + \
           d2r.sinD(alpha+beta) * d2r.sinD(alpha) * d2r.cosD(phi)

    cosR_corr = d2r.correct(cosR)
    R = d2r.acosD(cosR_corr)

    # problems with precision for 180 degrees
    
    if int(R*100.0) == 180.0:
        R = int(R*100.0)/100.0
    
    if R!=0.0 and R != 180.0 and alpha > 0.0:
        cosgamma = (d2r.cosD(alpha+beta) - d2r.cosD(alpha) * cosR) \
                  /(d2r.sinD(alpha) * d2r.sinD(R))
    else:
         cosgamma = 0.0

    cosgamma_corr = d2r.correct(cosgamma)
    gamma = d2r.acosD(cosgamma_corr)
    xp = R * d2r.sinD(gamma)
    
    if phi > 0.0:
        xp = -xp
    
    yp = -R * cosgamma

    return yp, xp

###########################################################################################

# simple test statement

if __name__ == "__main__":
   print mapphi(90,90,90)
