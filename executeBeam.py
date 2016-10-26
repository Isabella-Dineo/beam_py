#!/usr/bin/env python

import numpy as np
import argparse, os

parser = argparse.ArgumentParser(prog='Run the beam code multiple times.')
parser.add_argument('-f', type=str, default='outFile')
parser.add_argument('-it', type=int, default=1)
parser.add_argument('-p', type=float, default=1.0)
parser.add_argument('-nch', type=float, default=10)
parser.add_argument('-min_freq', type=float, default=0.03)
parser.add_argument('-bw', type=float, default=0.01)
parser.add_argument('-snr', type=int, default=100)
args = parser.parse_args()
filename = args.f
iterations = args.it
P = args.p
fmin = args.min_freq
bw = args.bw
nch = args.nch
snr = args.snr

for i in range(iterations):
    alpha = np.rad2deg(np.arccos(np.random.uniform()))
    iseed = np.random.randint(0, 90012410)
    beta = np.random.uniform(-20, 20)
    nc = 4
    npatch = 4
    os.system('generateBeam.py -alpha %.3f -beta %.3f -p %.3f -min_freq %.3f -chbw %.3f -nch %d -nc %d -npatch %d -snr %f -iseed %d -outfile %s' \
               % (alpha, beta, P, fmin, bw, nch, nc, npatch, snr, iseed, filename ))

