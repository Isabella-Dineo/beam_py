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
parser.add_argument('-getPlot', metavar='<boolean>', type=int, default=None, help='Option to produce the beam plot.')
parser.add_argument('-doFan', metavar='<boolean>', type=str, default=None, help='Option to produce the beam plot.')
parser.add_argument('-scatter', metavar='<boolean>', default=None, type=int, help='Option to include scattering effects.')
parser.add_argument('-dmFile', metavar='<fileName>', type=str, default='psrcatdm.dat', help='File containing known dm values from psrcat (used for scattering).')
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
        iseed = np.random.randint(0, 90012410)
    else:
        iseed = args.iseed
    
    if P == None:
        P = np.random.uniform(0.001, 2)
  
    beta = np.random.uniform(-20, 20)
    alpha = np.rad2deg(np.arccos(np.random.uniform()))
    nc = 4
    npatch = 4 
    os.system('generateBeam.py -alpha %.3f -beta %.3f -p %.3f -min_freq %.3f -chbw %.3f -nch %d -nc %d -npatch %d \
               -snr %d -iseed %d -outfile %s -scatter %s '\
               % (alpha, beta, P, fmin, bw, nch, nc, npatch, snr, iseed, filename, args.scatter))
