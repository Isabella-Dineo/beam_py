#!/usr/bin/env python 
import numpy as np
import matplotlib.pyplot as plt
import argparse

parser = argparse.ArgumentParser(description="Plot the delta dm standard deviation as a funstion of frequency")
parser.add_argument('-f', type=str, help="A file containing a mean, var, skew, kurtosis in column for each frequency.")
args = parser.parse_args()

# Load the txt:
m, v, s, k = np.loadtxt(args.f, unpack=True)
print v
# Define the center frequencies
fc = np.array([50, 150, 200, 1400])
#fc = np.array([150,50,1400, 200])

fig = plt.figure()
#plt.errorbar(fc, v, yerr=np.sqrt(v), fmt='o')
plt.plot(np.log10(fc), np.log10(np.sqrt(v)), 'ko')
plt.xlabel(r'$log(\nu_{c})$')
plt.ylabel(r'$log (\Delta DM_{var}$)')
#plt.xscale('log', nonposy='clip')
#plt.title(r'Mean $\Delta DM$ as a funtion of center frequency')
plt.savefig('meandm_freq_pc_noscr.png')
plt.show()

