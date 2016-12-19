#!/usr/bin/env python 

import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as stats


# Load a file with known dm value from psrcat
dm_rv = np.loadtxt('catdm.dat')

# Compute the probabilities (Normalize)
probs = [i/sum(dm_rv) for i in dm_rv]

# Define a probability distribution of function
iseed = 41 #
normdiscrete = stats.rv_discrete(values=(dm_rv, probs), seed=iseed, name='normdiscrete')
print 'mean = %6.4f, varience = %6.4f, skew = %6.4f, kurtosis = %6.4f ' \
%(normdiscrete.stats(moments='mvsk'))

# Draw a random dm to use for scatteing:
rand_dm = normdiscrete.rvs(size=1)
print 'Drawn sample: random dm value = %.4f' % rand_dm
