#!/usr/bin/env python 

import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as stats


# Load a file with known dm values from psrcat
psrcatdm = np.loadtxt('psrcatdm.dat')
# Create a histogram:
hist, bin_edges = np.histogram(psrcatdm, bins=9) #10bins forces empty bins ()
# Compare the hist values with matplotlib values (just because I don't trust my own judgement!)
plt.bar(bin_edges[:-1], hist, width=bin_edges[1]-bin_edges[0], color='red', alpha=0.5, label='numpy histogram' )
plt.hist(psrcatdm, bins=9, alpha=0.5, label='matplotlib histogram' )
plt.title('PSRCAT dm distribution')
plt.xlabel('dm (pc cm^-3)')
plt.legend()
# Arrange the dm into arrays containing values that fall in each of the bin intervals respectively
dm_arranged = []
for i in range(len(bin_edges) - 1): # N_bin_edges == 1 + N_hist
    dm_arranged.append(np.array(psrcatdm[np.where((psrcatdm > bin_edges[i]) & (psrcatdm <= bin_edges[i+1]))]))
    # dm_arranged will not include zero dm values that were used to replace null values from psrcat
# Find probabilities of each of the values in that interval:
probs = []
for j in range(len(bin_edges) - 1):
    probs.append([dm/sum(dm_arranged[j]) for dm in dm_arranged[j]])
    # probs will contain arrays of probabilities divided into intervals = bin numbers
# Same probability distribution pattern??:
plt.figure()
plt.hist(probs, bins=9)
plt.title('Probability distribution')
plt.xlabel('probabilities')
for l in range(len(probs)):
    print sum(probs[l])
#------------------------------------------------------------------------------
# Define an arbitrary distribution
#probs = np.hstack(np.asarray(probs))
rand_int = np.random.randint(0, len(dm_arranged))
normdiscrete = stats.rv_discrete(values=(dm_arranged[rand_int], probs[rand_int]))
rand_dm = normdiscrete.rvs(size=500)
plt.figure()
plt.hist(rand_dm, alpha=0.5, bins=9)
plt.title('Randomly samples dm values from rv_discrete')
plt.xlabel('dm (pc cm^-3)')
plt.show()

"""#-------------------------------------------------------------------------------
for k in range(len(dm_arranged))
    rand_int = np.random.randint(0, len(dm_arranged)) 
    print 'Interval: ', rand_int
    normdiscrete = stats.rv_discrete(values=(dm_arranged[k], probs[k]))
    print 'mean = ', (normdiscrete.stats(moments='m'))
    print 'varience = ', (normdiscrete.stats(moments='v'))
    print 'skew = ', (normdiscrete.stats(moments='s'))
    print 'kurtosis = ', (normdiscrete.stats(moments='k'))
    rand_dm = normdiscrete.rvs(size=500)
    print 'Drawn sample: random dm value = ' ,rand_dm
    plt.hist(rand_dm, alpha=0.5)
plt.show()

# Use numpy random choice 
#-------------------------------------------------------------------------------
iseed = 41 #
np.random.seed(iseed)
rand_dm_choice = np.random.choice(psrcatdm, p=probabilities, size=1)
print 'Using random choice : random dm value = %.4f' % rand_dm_choice"""
