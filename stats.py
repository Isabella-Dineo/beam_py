#!/usr/bin/env python

# Determine the mean and standard deviation of the DMs

import numpy as np
import matplotlib.pyplot as plt
import argparse
from scipy import stats
from scipy.optimize import curve_fit 

#def gauss(dat, mu, sigma, A):
#    """A gaussian function from dat with given mean, sigma and amplitude"""
#    return (A*np.exp(-(dat-mu)**2/sigma**2))
#def bimodal(dat, mu1, sigma1, A1, mu2, sigma2, A2):
#    """Function to sum two gaussian functions"""
#    return (A1*np.exp(-(dat-mu1)**2/sigma1**2) + A2*np.exp(-(dat-mu2)**2/sigma2**2))

parser = argparse.ArgumentParser(prog='Run the beam code multiple times.')
parser.add_argument('-f', nargs='+', type=str, help='files containing dm in 1 column only.')
parser.add_argument('-o', default='outFile.stats', type=str, help='outfile to write mean and variance.')
parser.add_argument('-model', default='patchy', type=str, help='Model used for the dm simulation.')
parser.add_argument('-fcen', type=float, help='Centre frequency (GHz).')
parser.add_argument('-outfile', type=str, help='png file name title.')
parser.add_argument('--fit', action='store_true', help='option to fit for a bimodal distribution')
args = parser.parse_args()


files = args.f
mean = []
varience = []
for fid in range(len(files)):
    dat = np.loadtxt(files[fid], delimiter=' ')
    mu = np.mean(dat)
    mean.append(mu)
    var = np.var(dat)
    varience.append(var)
    skewness = stats.skew(dat)
    kurtosis = stats.kurtosis(dat) 
    # write out important parameters
    mean_var = np.asarray([mu, var, skewness, kurtosis])
    f = open(args.o, 'a')
    f.write(' '.join([str(item) for item in mean_var]) + ' \n')
    # Plot a dm distribution
    fig = plt.figure(figsize=(10, 10))
    x,y,_=plt.hist(dat, bins=100)
    plt.yscale('log', nonposy='clip')
    plt.title('$\Delta$DM distribution: ' + str(args.o) + ' band', fontsize=18)
    plt.xlim(-21, 21)
    plt.title('$\Delta$DM distribution: ' + str(args.o) + ' band')
    #if abs(np.max(dat)) > abs(np.min(dat)):
    #    plt.xlim(-np.max(dat), np.max(dat))
    #elif abs(np.max(dat)) < abs(np.min(dat)):
    #    plt.xlim(-abs(np.min(dat)), abs(np.min(dat)))
    plt.xlabel(r'$\Delta$DM', fontsize=18)
    plt.ylabel('log(samples)', fontsize=18)
    plt.xlabel(r'$\Delta$DM')
    plt.ylabel('log(samples)')
 
    # Fitting a bimodal distribution:
#    if args.fit:
#        expected=(mu, np.sqrt(var), np.max(x), mu, np.sqrt(var), np.max(x)) 
#        print len(x), len(y)
#        params, cov = curve_fit(bimodal, x, y[:-1],expected) 
#        sigma = np.sqrt(np.diag(cov))
#        plt.plot(x, (bimodal(x,*params)), color='red', lw=3, label='model')
#        plt.legend()
    fig.savefig(args.o + str(fid) + '.png')
    plt.show()
print 'Mean = ', mean
print 'Varience = ', varience

