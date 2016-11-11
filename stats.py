#!/usr/bin/env python

# Determine the mean and standard deviation of the DMs

import numpy as np
import matplotlib.pyplot as plt
import argparse

parser = argparse.ArgumentParser(prog='Run the beam code multiple times.')
parser.add_argument('-f', type=str, default='input file')
args = parser.parse_args()

dat = np.loadtxt(args.f, delimiter=' ')
nu = np.mean(dat)
var = np.var(dat)

fig = plt.figure(figsize=(10,10))
plt.hist(dat, bins=200)
plt.title('Distribution of DMs_ ' + args.f+'_mean_'+str(nu)+'_var_'+str(var) )
fig.savefig(args.f+'.png')
plt.show()

