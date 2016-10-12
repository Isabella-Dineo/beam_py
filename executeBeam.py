#!/usr/bin/env python 

import numpy as np
import os, argparse
parser = argparse.ArgumentParser()
parser.add_argument('-f')
parser.add_argument('-it', type = int, default = 1)
parser.add_argument('-p', type = float, default = 1.0)
args = parser.parse_args()
filename = args.f
iterations = args.it
period = args.p
for i in range(iterations):
    alpha = np.rad2deg(np.arccos(np.random.uniform()))  # inclination angle in degrees
    beta = np.random.uniform(-20, 20)                   # impact parameter in degrees
#    P = np.random.uniform(0.15, 1.5)                    # Period in seconds
    P = period
    fmin = 0.065                                        # minimum frequency in GHz
    bw = 0.01                                           # channel bandwidth in GHz
    nch = 10                                            # Number of channels
    nc = 4					        # Number of components
    npatch = 4                                          # Number of active patches
    iseed = np.random.randint(0,90012410)
    os.system('generate_profile.py -alpha %.3f -beta %.3f -p %.3f -min_freq %.3f -chbw %.3f -nch %d -nc %d -npatch %d -snr 10000 -iseed %d -outfile %s' \
               % (alpha, beta, P, fmin, bw, nch, nc, npatch, iseed, filename ))
    os.system('plot_profile.py -outfile %s' % filename)
