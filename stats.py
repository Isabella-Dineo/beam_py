#!/usr/bin/env python

# Determine the mean and standard deviation of the DMs

import numpy as np
import matplotlib.pyplot as plt
import argparse

parser = argparse.ArgumentParser(prog='Run the beam code multiple times.')
parser.add_argument('-f', nargs='+', type=str, help='files containing dm in 1 column only.')
parser.add_argument('-o', default='outFile.stats', type=str, help='outfile to write mean and variance.')
parser.add_argument('-model', default='patchy', type=str, help='Model used for the dm simulation.')
parser.add_argument('-fcen', type=float, help='Centre frequency (GHz).')
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
    # write out important parameters
    mean_var = np.asarray([nu, var, args.fcen])
    f = open(args.o, 'a')
    f.write(' '.join([str(item) for item in mean_var]) + ' \n')
    # Plot a dm distribution
    fig = plt.figure(figsize=(10,10))
    plt.hist(dat, bins=100)
    plt.title('Distribution of DMs_%s : Mean %.5f, Var %.5f' %(files[fid], nu, var))
    fig.savefig(files[fid] + str(fid) + '.png')
    plt.show()
print 'Mean = ', mean
print 'Varience = ', varience

