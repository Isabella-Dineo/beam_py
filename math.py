import numpy as np

# Functions to compute trig calculations in radians 

def cosD(ang):
	angR = np.deg2rad(ang)
	return np.cos(angR)

def sinD(ang):
	return np.sin(np.deg2rad(ang))

def tanD(ang):
	return np.tan(np.deg2rad(ang))

def acosD(ang):
	return np.arccos(np.deg2rad(ang))

def asinD(ang):
	return np.arccos(np.deg2rad(ang))

def atanD(ang):
	return np.arctan(np.deg2rad(ang))

#########################################################

# correct for round-off errors (tolerence???)

def correct(x):
	tol = 1.0*np.exp(-7)
	np.sign = 1.0
	
	sign = np.sign(x) # checks the sign of x, (Can I do that to replace the "if" statement originally used?)

	y = np.abs(x)

	if np.abs(y-1.0) < tol:
		y = 1.0
	elif np.abs(y-0.0) < tol:  # why "y-0.0" and not just "y" only?
		y = 0.0

	return y*sign

if __name__ == "__main__":
	
	print cosD(90.)
	print cosD(60.)

