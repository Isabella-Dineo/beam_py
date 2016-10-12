#!/usr/bin/env python 

import numpy as np
import os

for i in range(1):
    alpha = np.rad2deg(np.arccos(np.random.uniform()))  # inclination angle in degrees
    beta = np.random.uniform(-20, 20)                   # impact parameter in degrees
    P = np.random.uniform(0.15, 1.5)                    # Period in seconds
    fmin = 0.065                                        # minimum frequency in GHz
    bw = 0.01                                           # channel bandwidth in GHz
    nch = 2                                           # Number of channels
    nc = 4					        # Number of components
    npatch = 4                                          # Number of active patches
    iseed = np.random.randint(0,9000)
    os.system('generate_profile.py -alpha %.3f -beta %.3f -p %.3f -min_freq %.3f -chbw %.3f -nch %d -nc %d -npatch %d -snr 10000 -iseed %d' \
               % (alpha, beta, P, fmin, bw, nch, nc, npatch, iseed ))
    os.system('plot_profile.py')
