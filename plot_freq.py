#!/usr/bin/env python 
import numpy as np
import matplotlib.pyplot as plt
import argparse

parser = argparse.ArgumentParser(description="Plot the delta dm standard deviation as a funstion of frequency")
parser.add_argument('-f', type=str, help="A file containing a mean, var, skew, kurtosis in column for each frequency.")
args = parser.parse_args()

# Load the txt:
stats = np.loadtxt(args.f)
# Define the center frequencies
freq = stats[:, 0]
var = stats[:, 2]
print 'freq', stats[:,0] 
print 'var', stats[:,2]

fig = plt.figure()
plt.plot(np.log10(stats[:,0]), np.log10(stats[:,2]), 'k.')
plt.xlabel(r'$log(\nu_{c})$', fontsize=18)
plt.ylabel(r'$log (\sigma^2_{\Delta DM}$)', fontsize=18)
plt.savefig('var_freq_hc_noscr.png')
plt.show()

