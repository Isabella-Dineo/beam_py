#!/usr/bin/env python

# Program to fit gaussians to profiles 
# from the epn database.

import numpy as np
import csv
import matplotlib.pyplot as plt
import argparse
import os, sys, glob, time
from scipy.optimize import curve_fit, leastsq

#----------------------
# Define functions:
#----------------------
def getProfile(txtFile):
    profile = []
    with open(txtFile) as f:
        csvFile = csv.reader(f, delimiter=' ', skipinitialspace=True)
        for line in csvFile:
            profile.append(line[3])
    x = np.linspace(-1,1,num=len(profile))
    profile = np.asarray(profile)
    return x, profile

def gaussFit(data, amp, dataMean, dataVar, offset):
    return  amp * np.exp(-(data - dataMean)**2/(2*dataVar**2)) + offset

def gaussian1(data, amp, center, width, offset):
    return amp * np.exp(-(data - center)**2/(2*width**2)) + offset     

def gaussian3(data, h1, c1, w1, h2, c2, w2, h3, c3, w3, offset):
    g1 = gaussian1(data, h1, c1, w1, offset=0)
    g2 = gaussian1(data, h2, c2, w2, offset=0)
    g3 = gaussian1(data, h3, c3, w3, offset=0)
    sumGauss = g1 + g2 + g3 + offset
    return sumGauss

def gaussian2(data, h1, c1, w1, h2, c2, w2, offset):
    gauss = gaussian3(data, h1, c1, w1, h2, c2, w2, 0, 0, 1, offset)
    return gauss

errorFunc1 = lambda p, x, y: (gaussian1(x, *p) - y)**2
errorFunc2 = lambda p, x, y: (gaussian2(x, *p) - y)**2
errorFunc3 = lambda p, x, y: (gaussian3(x, *p) - y)**2
    
#----------------------
# Command Line parser:
#----------------------
parser = argparse.ArgumentParser(description='Program to find relationship between pulsar profile widths and frequency for profiles from the EPN database.',
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('--txt', default=None, type=str, help='Give a text file with stokes I values.')
parser.add_argument('--direc', default=None, type=str, help='Or give location of the txt files.')
parser.add_argument('--height', nargs='+', default=None, type=float, help='Give 1 or more peak value of the gaussians to fit.')
parser.add_argument('--center', nargs='+', default=None, type=float, help='Give 1 or more centers of the peaks.')
parser.add_argument('--width', nargs='+', default=None, type=float, help='Give 1 or more widths of the profiles.')
parser.add_argument('--freq', type=str, default=None, help='Give freq of the profile.')
args = parser.parse_args()

#----------------------
# Import the text file:
#----------------------

if args.txt != None:
    if args.direc !=None:
        files = np.array(args.direc + args.txt)
    
    else:
        files = args.txt
else:
    if args.direc !=None:
        files = np.array(glob.glob(args.direc + '/*.txt'))
        if files.size > 1:
            print 'Select a file from the files below: \n', files
            files = raw_input()
    else:
        print 'Give a text file with stokes parameters!'
        sys.exit(1)

if args.freq == None:
    print 'Give the frequency of the profile!'
    sys.exit(1)

# Guess the peaks and perform least square fitting
# Find mean, sigma and amplitude of each profile.
phase, prof = getProfile(files)
prof = prof.astype(np.float)
profileMean = np.mean(prof)
varience = np.var(prof)
FWHM = 2.355 * np.sqrt(varience)
#guess1 = np.asarray([1.0, 0.01, 0.1, 0])#, 0.81, 0.02, FWHM, 0])
#guess2 = np.asarray([1.0, 0.01, 0.1, 0.81, 0.02, 0.01, 0])

if args.height == None:
    print "Give the height, center, and width of a gaussian to fit!"
    sys.exit(1)

else:
#    phase, prof = getProfile(files)
#    prof = prof.astype(np.float)
#    profileMean = np.mean(prof)
#    varience = np.var(prof)
#    FWHM = 2.355 * np.sqrt(varience) 
    h = args.height
    c = args.center
    w = args.width
    fig = plt.figure()
    plt.plot(phase, prof, lw=3, c='b', label='proflile')
    plt.xlabel('Phase')
    plt.ylabel('Flux density')
    plt.legend()

    if len(h) == 1:
        print "Fitting 1 Gaussian to data..."
        guess1 = np.asarray([h[0], c[0]], w[0], 0)
        popt1, pcov1 = curve_fit(gaussian1, phase, prof, p0=guess1)
        fit = gaussian1(phase, *popt1)
        residual = prof - fit 
        plt.plot(phase, fit, lw=1, c='r', ls='--', label='Fit of 1 Gaussian.')
        fig = plt.figure()
        plt.plot(phase, prof, lw=3, c='b', label='proflile')
        plt.xlabel('Phase')
        plt.ylabel('Flux density')
        plt.legend()
        plt.figure()
        plt.plot(phase, residual)
        plt.title('Residual plot')

    if len(h) == 2:
        print "Fitting %d Gaussians to data..." %len(h)
        #for hid in range(len(h)):
        guess2 = np.asarray([h[0], c[0], w[0], h[1], c[1], w[1], 0])
        popt2, pcov2 = curve_fit(gaussian2, phase, prof, p0=guess2)
        fit = gaussian2(phase, *popt2)
        residual = prof - fit
        plt.plot(phase, fit, lw=1, c='r', ls='--', label='Fit of 2 Gaussians.')
        plt.xlabel('Phase')
        plt.ylabel('Flux density')
        plt.legend()
        fig.savefig(str(args.freq) + '_' + str(time.time()) + '.png')
        plt.figure()
        plt.plot(phase, residual)
        plt.title('Residual plot')

    if len(h) == 3:
        print "Fitting a multi-peak gaussian..."
        popts = []
        guess3 = np.asarray([h[0], c[0], w[0], h[1], c[1], w[1], h[2], c[2], w[2], 0])
        popt3, pcov3 = curve_fit(gaussian3, phase, prof, p0=guess3)
        fit = gaussian3(phase, *popt3)
        residual = prof - fit 
        plt.plot(phase, fit, lw=1, c='r', ls='--', label='Fit of 1 Gaussian')
        plt.xlabel('Phase')
        plt.ylabel('Flux density')
        plt.legend()
        fig.savefig(str(args.freq) + '_' + str(time.time()) + '.png')
        plt.figure()
        plt.plot(phase, residual, 'g:')
        plt.title('Residual plot')
       
        #popts.append(popt)
        
#        plt.figure() 
#        plt.plot()
#        popts_sum = np.sum(popts, axis=0)
#        plt.plot(phase, prof, lw=5, c='b', label='proflile')
        
        
    #guess1 = np.asarray([1.0, 0.01, 0.1, 0])#, 0.81, 0.02, FWHM, 0])
    #guess2 = np.asarray([1.0, 0.01, 0.1, 0.81, 0.02, 0.01, 0])
    #popt1, pcov = curve_fit(gaussian1, phase, prof, p0=guess1)
    #popt2, pcov = curve_fit(gaussian2, phase, prof, p0=guess2)
    #plt.plot(phase, prof, lw=5, c='b', label='proflile')
    #plt.plot(phase, gaussian1(phase, *popt1), lw=1, c='r', ls='--', label='Fit of 1 Gaussian')
    #plt.plot(phase, gaussian2(phase, *popt2), lw=3, c='g', ls='-', label='Fit of 2 Gaussians')
    
    #plt.figure()
    #optim2, success = leastsq(errorFunc2, guess2[:], args=(phase[:], prof[:]))
    #optim1, success = leastsq(errorFunc1, guess1[:], args=(phase[:], prof[:]))
    #plt.plot(phase, prof, lw=5, c='g', label='Profile')
    #plt.plot(phase, gaussian2(prof, *guess2), lw=3, c='b', ls='-', label='Fit of 2 gaussians')
    #plt.plot(phase, gaussian1(prof, *optim1), lw=1, c='r', ls='--', label='Fit of 1 gaussian') 
    plt.show()
