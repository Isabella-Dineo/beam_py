#!/usr/bin/env python

# Determine the mean and standard deviation of the DMs

import numpy as np

dataFile = np.loadtxt('DM.dat', delimiter=' ')
print np.mean(dataFile)
print np.var(dataFile)

