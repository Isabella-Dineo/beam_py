#!/usr/bin/env python

# Program to plot profiles from the EPN database
# and fit a gausssian to see the relationship
# between the  

import numpy as np
import csv
import matplotlib.pyplot as plt
import argparse
import os, sys, glob, time

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
    return profile, x

#----------------------
# Command Line parser:
#----------------------
parser = argparse.ArgumentParser(description='Program to find relationship between pulsar profile widths and frequency for profiles from the EPN database.', 
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('--txt', metavar='<.txt file>', default=None, type=str, help='Give a text file with stokes I values.')
parser.add_argument('--direc', metavar='<directory>', default=None, type=str, help='Or give location of the txt files.')
parser.add_argument('--anim', default=None, help='Convert the profiles to gif plot.')
parser.add_argument('--title', metavar='<title>', type=str, help='Give title to the profile plot.')
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

# Find mean, sigma and amplitude of each profile.


prof, phase = getProfile(files)
prof = prof.astype(np.float)
profileMean = np.mean(prof)
print 'Mean: ', profileMean
varience = np.var(prof)
print 'Var: ', varience
amplitude = np.max(prof)
print 'Amplitude: ', amplitude
fig = plt.figure(figsize=(20,20)) 
plt.title(args.title)
plt.xlabel('phase')
plt.ylabel('Flux density')
plt.plot(phase.astype(np.float), prof)
fig.savefig(args.title + '.png')

if args.anim != None:
    os.system("convert -delay 100 -loop 0 *.png animation.gif")
    os.system("gnome-open animation.gif")
    print 'Created animation "animation.gif"'
