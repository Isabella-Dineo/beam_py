#!/usr/bin/env python

import numpy as np
import argparse, os
import time

parser = argparse.ArgumentParser(prog='Run the beam code multiple times.', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('-it', metavar='<iterarions>', type=int, default=1, help='Number of iterations.')
parser.add_argument('-p', metavar='<period>', type=float, default=None, help='Period in seconds.')
parser.add_argument('-nch', metavar='<nch>', type=float, default=10, help='Number of channels.')
parser.add_argument('-min_freq', metavar='<minFreq>', type=float, default=0.03, help='Start frequency (GHz).')
parser.add_argument('-bw', metavar='<bandwidth>', type=float, default=0.01, help='Channel bandwidth (GHz).')
parser.add_argument('-snr', metavar='<SNR>', type=int, default=100, help='Signal to noise ratio.')
args = parser.parse_args()
iterations = args.it
P = args.p
fmin = args.min_freq
bw = args.bw
nch = args.nch
snr = args.snr
for i in range(iterations):
    nc = 4
    npatch = 4 
    alpha = np.rad2deg(np.arccos(np.random.uniform()))
    iseed = time.time()
    os.system('generateBeam.py -alpha %f -p %f -min_freq %f -chbw %f -nch %d -nc %d -npatch %d \
               -snr %f -iseed %d --outfile --diagnostic --random_beta --template_matching --getPlot'\
               % (alpha, P, fmin, bw, nch, nc, npatch, snr, iseed))
    command = 'generateBeam.py -alpha %f -p %f -min_freq %f -chbw %f -nch %d -nc %d -npatch %d -snr %f -iseed %d --outfile --diagnostic --random_beta --template_matching --getPlot' % (alpha, P, fmin, bw, nch, nc, npatch, snr, iseed)
    f = open('executed_commands.txt', 'a')
    f.write(' '.join([str(command)]) + ' \n')
