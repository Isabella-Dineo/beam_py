#!/usr/bin/env python

import numpy as np
import argparse, os
import time

parser = argparse.ArgumentParser(prog='Run the beam code multiple times.', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('-f', metavar='<fileName>', type=str, default='outFile')
parser.add_argument('-it', metavar='<iterarions>', type=int, default=1, help='Number of iterations.')
parser.add_argument('-p', metavar='<period>', type=float, default=None, help='Period in seconds.')
parser.add_argument('-nch', metavar='<nch>', type=float, default=10, help='Number of channels.')
parser.add_argument('-min_freq', metavar='<minFreq>', type=float, default=0.03, help='Start frequency (GHz).')
parser.add_argument('-bw', metavar='<bandwidth>', type=float, default=0.01, help='Channel bandwidth (GHz).')
parser.add_argument('-snr', metavar='<SNR>', type=int, default=100, help='Signal to noise ratio.')
parser.add_argument('-iseed', metavar='<seed>', type=int, default=None, help='Seed for the random munber generator.')
parser.add_argument('-dmFile', metavar='<fileName>', type=str, default='psrcatdm.dat', help='File containing known dm values from psrcat (used for scattering).')
parser.add_argument('--getPlot', action="store_true", default=True, help='Option to produce the beam plot.')
parser.add_argument('--scatter', action="store_true", default=False, help='Option to include scattering effects.')
parser.add_argument('--doFan', action="store_true", help='Option to use the fanBea model.')
parser.add_argument('--doHC', action="store_true", help='Option to use the hollow cone beam model.')
args = parser.parse_args()
filename = args.f
iterations = args.it
P = args.p
fmin = args.min_freq
bw = args.bw
nch = args.nch
snr = args.snr
for i in range(iterations):
    if args.iseed == None:
        iseed = int(time.time())
    else:
        iseed = args.iseed
    np.random.seed(iseed) 
    if P == None:
        P = np.random.uniform(0.1, 2)
  
    alpha = np.rad2deg(np.arccos(np.random.uniform()))
    nc = 4
    npatch = 4 
    os.system('generateBeam.py -alpha %.3f -p %.3f -min_freq %.3f -chbw %.3f -nch %d -nc %d -npatch %d \
               -snr %d --outfile --doHC --getPlot --diagnostic --randombeta --scatter'\
               % (alpha, P, fmin, bw, nch, nc, npatch, snr))
