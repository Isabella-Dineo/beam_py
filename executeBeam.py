#!/usr/bin/env python

import numpy as np
import argparse, os

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
parser.add_argument('--doFan', action="store_true", help='Option to produce the beam plot.')
parser.add_argument('--scatter', action="store_true", default=False, help='Option to include scattering effects.')
args = parser.parse_args()
filename = args.f
iterations = args.it
P = args.p
fmin = args.min_freq
bw = args.bw
nch = args.nch
snr = args.snr
print args.getPlot
print args.scatter
print args.doFan
for i in range(iterations):
    if args.iseed == None:
        iseed = np.random.randint(0, 90012410)
    else:
        iseed = args.iseed
    
    if P == None:
        P = np.random.uniform(0.1, 2)
  
    beta = np.random.uniform(-20, 20)
    alpha = np.rad2deg(np.arccos(np.random.uniform()))
    nc = 4
    npatch = 4 
    os.system('generateBeam.py -alpha %.3f -beta %.3f -p %.3f -min_freq %.3f -chbw %.3f -nch %d -nc %d -npatch %d \
               -snr %d -iseed %d -outfile %s --getPlot --doFan'\
               % (alpha, beta, P, fmin, bw, nch, nc, npatch, snr, iseed, filename))
