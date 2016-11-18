#!/usr/bin/env python

# Determine the mean and standard deviation of the DMs

import numpy as np
import matplotlib.pyplot as plt
import argparse

parser = argparse.ArgumentParser(prog='Run the beam code multiple times.')
parser.add_argument('-f', nargs='+', type=str, help='files containing dm in 1 column only.')
parser.add_argument('-o', default='outFile.stats', type=str, help='outfile to write mean and variance.')
parser.add_argument('-model', default='patchy', type=str, help='Model used for the dm simulation.')
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
    fig = plt.figure(figsize=(10,10))
    plt.hist(dat, bins=100)
    plt.title('Distribution of DMs_%s: Mean %.3f, Var %.3f' %(files[fid], nu, var))
    fig.savefig(files[fid] + str(fid) + '.png')
    plt.show()
print 'Mean = ', mean
print 'Varience = ', varience


bandCenters = np.array([110, 200, 580, 1275])
fig2 = plt.figure(figsize=(10,10))
plt.plot(bandCenters, varience, 'r*')
plt.xlabel('Center frequency (MHz)')
plt.ylabel('DM varience')
fig2.savefig('Var_Distribution_with_frequency %s .png' % str(args.model))
plt.show()

