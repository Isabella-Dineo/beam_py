#!/usr/bin/env python

import numpy as np
import argparse, os

parser = argparse.ArgumentParser(prog='Run the beam code multiple times.')
parser.add_argument('-f', type=str, default='outFile')
parser.add_argument('-it', type=int, default=1, help='Number of iterations.')
parser.add_argument('-p', type=float, default=None, help='Period in seconds.')
parser.add_argument('-nch', type=float, default=10, help='Number of channels.')
parser.add_argument('-min_freq', type=float, default=0.03, help='Start frequency.')
parser.add_argument('-bw', type=float, default=0.01, help='Channel bandwidth.')
parser.add_argument('-snr', type=int, default=100, help='Signal to noise ratio.')
parser.add_argument('-iseed', type=int, default=None, help='Seed for the random munber generator.')
#parser.add_argument('-getPlot', type=str, default=None, help='Option to produce the beam plot.')
parser.add_argument('-doFan', type=str, default=None, help='Option to produce the beam plot.')
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
    os.system('generateBeam.py -alpha %.3f -beta %.3f -p %.3f -min_freq %.3f -chbw %.3f -nch %d -nc %d -npatch %d -snr %f \
               -iseed %d -outfile %s -doFan %s'\
               % (alpha, beta, P, fmin, bw, nch, nc, npatch, snr, iseed, filename, args.doFan ))
