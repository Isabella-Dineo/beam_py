#!/usr/bin/env python

# Determine the mean and standard deviation of the DMs

import numpy as np
import matplotlib.pyplot as plt
import argparse
from scipy import stats

parser = argparse.ArgumentParser(prog='Run the beam code multiple times.')
parser.add_argument('-f', nargs='+', type=str, help='files containing dm in 1 column only.')
parser.add_argument('-o', default='outFile.stats', type=str, help='outfile to write mean and variance.')
parser.add_argument('-model', default='patchy', type=str, help='Model used for the dm simulation.')
parser.add_argument('-fcen', type=float, help='Centre frequency (GHz).')
parser.add_argument('-outfile', type=str, help='png file name title.')
args = parser.parse_args()


files = args.f
mean = []
varience = []
for fid in range(len(files)):
    dat = np.loadtxt(files[fid], delimiter=' ')
    nu = np.mean(dat)
    mean.append(nu)
    var = np.var(dat)
    varience.append(var)
    skewness = stats.skew(dat)
    kurtosis = stats.kurtosis(dat) 
    # write out important parameters
    mean_var = np.asarray([nu, var, skewness, kurtosis])
    f = open(args.o, 'a')
    f.write(' '.join([str(item) for item in mean_var]) + ' \n')
    # Plot a dm distribution
    fig = plt.figure(figsize=(10, 10))
    plt.hist(dat, bins=100)
    plt.yscale('log', nonposy='clip')
    plt.title('$\Delta$DM distribution: ' + str(args.o) + ' band', fontsize=18)
    plt.xlim(-0.075, 0.075)
    #if abs(np.max(dat)) > abs(np.min(dat)):
    #    plt.xlim(-np.max(dat), np.max(dat))
    #elif abs(np.max(dat)) < abs(np.min(dat)):
    #    plt.xlim(-abs(np.min(dat)), abs(np.min(dat)))
    plt.xlabel(r'$\Delta$DM', fontsize=18)
    plt.ylabel('log(samples)', fontsize=18)
    fig.savefig(args.o + str(fid) + '.png')
    plt.show()
print 'Mean = ', mean
print 'Varience = ', varience

