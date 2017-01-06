#!/usr/bin/env python 

import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as stats


# Load a file with known dm value from psrcat
psrcatDM = np.loadtxt('psrcatdm.dat')

# Compute the probabilities (Normalize)
probabilities = [i/sum(psrcatDM) for i in psrcatDM]

print sum(probabilities)

#------------------------------------------------------------------------------
# Define stats
#------------------------------------------------------------------------------
# Define a probability distribution of function
normdiscrete = stats.rv_discrete(values=(psrcatDM, probabilities), name='normdiscrete')
print 'mean = %6.4f, varience = %6.4f, skew = %6.4f, kurtosis = %6.4f ' \
%(normdiscrete.stats(moments='mvsk'))
# Draw a random dm to use for scatteing:
rand_dm = normdiscrete.rvs()
print 'Drawn sample: random dm value = %.4f' % rand_dm
#-------------------------------------------------------------------------------
# Use numpy random choice 
#-------------------------------------------------------------------------------
iseed = 41 #
np.random.seed(iseed)
rand_dm_choice = np.random.choice(psrcatDM, p=probabilities, size=1)
print 'Using random choice : random dm value = %.4f' % rand_dm_choice
