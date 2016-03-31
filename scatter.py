#!/usr/bin/python

import numpy as np

def scatter(prof,dm,period,freq):
    """
    Function to perform convolution for a scattered profile. If dm = 0, the function
    returns the profile without convolution. 
    
    Args:
    -----

        prof   : an array of the profile
        dm     : dispersion measure in cm^-3 pc.
        period : rotational period in seconds.
        freq   : observing frequency in GHz.

    Returns:
    --------

        conv   : convolution of the profile.
    """

    if dm != 0:
        
        npatches = float(prof.shape[0])
        tbin = (period * 1000) / npatches    # tbin = time of each bin in ms.
        log_tau = -6.46 + 0.154 * np.log10(dm) + 1.07 * (np.log10(dm))**2 - 3.86 * np.log10(freq)
#       tau = scattering time scale as in Bhat et al. (2004) in ms.
        tau = 10**log_tau
        tcrit = -tau / np.log(0.01)
        bincrit = tcrit / tbin

        verbose = raw_input('verbose (y/n)? : ')
        if verbose is 'y':
           print "Scattering for tau = " + str(tau) + ", dm = " + str(dm) + \
           ", period = " + str(period) + ", freq = " + str(freq) + ", tbin = " + str(tbin)
        

        if bincrit > 10 * npatches: bincrit = 10 * npatches

        if bincrit < npatches: bincrit = npatches

        if verbose is 'y': print "Critical bin: " + str(bincrit)
        
#       Convolution:
#       1. Transpose the profile:
        
        tprof = prof[::-1]

#       2. Slide along the time axis for as long as the scatter response is still high enough
#          (tb = straight bin number, tbf = folded bin number)
        tb = 0
        tbf = 0
        
        conv = np.zeros(int(npatches))
        scr = []
        while tb < bincrit + npatches - 1:

            #print tbf, conv.shape, npatches, bincrit+npatches-1
            for i in np.arange(int(npatches) - 1):
               
                if tb + i - 2 > npatches:
                    frac = (tb + i - 2 - npatches) * tbin / tau
                    scr.append(np.exp(-frac))
                else: 
                    scr.append(0.)

                conv[tbf] += tprof[i] * scr[i]

            tb += 1
            tbf += 1
            if tbf > int(npatches) - 1: tbf = 0 

        conv = conv[:]         
    
    if dm == 0.: 
        print('No scattering for dm == 0; returning the profile without convolution!')

    return conv

###############################################################################################
#simple test statement

if __name__ == "__main__":
    
    prof = np.arange(12)
    print scatter(prof,2.,0.5,1.4)
    
    
