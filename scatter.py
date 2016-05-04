#!/usr/bin/python

import numpy as np

def scatter(prof,dm,period,freq):
    """
    Function to perform convolution for a scattered profile. If dm = 0, the function
    returns the input profile without convolution. 
    
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
        
        npts = float(prof.shape[0])
        tbin = (period * 1000) / npts    # tbin = time of each bin in ms.
        log_tau = -6.46 + 0.154 * np.log10(dm) + 1.07 * (np.log10(dm))**2 - 3.86 * np.log10(freq)
#       tau = scattering time scale as in Bhat et al. (2004) in ms.
        tau = 10**log_tau
        tcrit = -tau / np.log(0.01)
        bincrit = tcrit / tbin

        verbose = raw_input('verbose (y/n)? : ')
        if verbose is 'y':
           print "Scattering for tau = " + str(tau) + ", dm = " + str(dm) + \
           ", period = " + str(period) + ", freq = " + str(freq) + ", tbin = " + str(tbin)
        

        if bincrit > 10 * npts: bincrit = 10 * npts

        if bincrit < npts: bincrit = npts

        if verbose is 'y': print "Critical bin: " + str(bincrit)
        
#       Convolution:
#       1. Transpose the profile:
        
        tprof = prof[::-1]

#       2. Slide along the time axis for as long as the scatter response is still high enough
#          (tb = straight bin number, tbf = folded bin number)
        tb = 0
        tbf = 0
        
        conv = np.zeros(int(npts))
        #scr = np.zeros(int(npts))
        while tb < bincrit + npts - 1:

            #print tbf, conv.shape, npts, bincrit+npts-1
            #the loops slows down the program
            for i in np.arange(int(npts) - 1):
               
                if tb + i - 1 > npts:
                    frac = (tb + i - 1 - npts) * tbin / tau
                    scr = np.exp(-frac)
                else: 
                    scr = 0.

                #conv[tbf] = tprof[i] * scr[i]
                conv[tbf] += tprof[i] * scr

            tb += 1
            tbf += 1
            if tbf > int(npts) - 1: tbf=1 
       
    
    if dm == 0.: 
        print('No scattering for dm == 0; returning the profile without convolution!')
        conv = prof 
    return conv

###############################################################################################
#simple test statement

if __name__ == "__main__":
    
    prof = np.arange(200)
    aa = scatter(prof,200,0.5,0.3)
    print aa
    print len(aa)
    
