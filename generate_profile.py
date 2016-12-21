#!/usr/bin/env python 

import numpy as np 
import matplotlib.pyplot as plt 
import matplotlib.cm as cm
from scipy import constants
import argparse
import pickle
import os
#######################
# DEFINE FUNCTIONS.
#######################

#====================================================================================================================================================
#							PRECISION ERRORS		
#====================================================================================================================================================
def correct(x):
    """Function to correct for precision errors.
       
       Args:
       -----
       x    : ndarray.
       
       Returns:
       --------
       y    : ndarray correction of x.
    """

#   Tolarance:
    tol = 1.0e-07
    sign = np.sign(x)
    y = np.abs(x)
    
#   correct for precision:
    for i in np.arange(np.size(x)):
        if np.abs(y[i] - 1.0) < tol:
            y[i] = 1.0
        elif np.abs(y[i] - 0.0) < tol:
            y[i] = 0.0

    return y * sign

#====================================================================================================================================================
#							ROTATIONAL AXES		
#====================================================================================================================================================
def mapphi(alpha, beta, phi):
    """Function to map the rotational axis:
       
       Args:
       -----
       alpha       : magnetic inclination angle w.r.t the rotatinal axis (degrees)
       beta        : line of sight closest approach to the magnetic axis (degrees)
       
       Returns:
       --------
       xlos, ylos  : The coordinates of the rotatinal plane; both the size of phi.
       
    """
    cosR = np.cos(np.deg2rad(alpha+beta)) * np.cos(np.deg2rad(alpha))
    cosR += np.sin(np.deg2rad(alpha+beta)) * np.sin(np.deg2rad(alpha)) * np.cos(np.deg2rad(phi))

    R = np.arccos(correct(cosR))

    # problems with precision for 180 degrees
    cosgamma = np.zeros_like(R)

    for i in np.arange(len(R)):
        if int(R[i]*100.0) == 180.0:
            R[i] = int(R[i]*100.0)/100.0
        if R[i] != 0.0 and R[i] != 180.0 and alpha > 0.0:
            cosgamma[i] = (np.cos(np.deg2rad(alpha+beta)) - np.cos(np.deg2rad(alpha)) * cosR[i]) \
                  /(np.sin(np.deg2rad(alpha)) * np.sin(R[i]))
        else:
             cosgamma[i] = 0.0

    cosgamma_corr = correct(cosgamma)
    gamma = np.arccos(cosgamma_corr)
    xp = R * np.sin(gamma)
    
    for i in np.arange(len(phi)):
        if phi[i] > 0.0:
            xp[i] = -xp[i]

    yp = -R * np.cos(gamma)
 
    return np.rad2deg(xp)[::-1], np.rad2deg(yp)[::-1] #clockwise los
#====================================================================================================================================================
#							LINE OF SIGHT
#====================================================================================================================================================
def los(alpha, beta, res):
    """Function to determine the line of sight cut across the beam.
    
       Args:
       -----
       alpha       : Inclination angle (degrees).
       beta        : Impact parameter (degrees).
       
       Returns:
       --------
       xlos        : The line of sight x-coordinates (degrees).
       ylos        : The line of sight y-coordinates (degrees).
       thetalos    : The line of sight angle in degrees (degrees).
    """

#   rotational phase:
    phi = np.linspace(-180, 180, num=res, endpoint=True)
    
#   line of sight x,y plane:
    xlos, ylos = mapphi(alpha, beta, phi)
    thetalos = np.arctan2(ylos, xlos) * (180 / np.pi) - 90.0
    thetalos = np.abs(thetalos)
    #for i in np.arange(len(thetalos)):
    #    if thetalos[i] < 0:
    #        thetalos[i] = -thetalos[i]       
            
    return xlos, ylos, thetalos

#====================================================================================================================================================
#							EMISSION HEIGHTS
#====================================================================================================================================================
def emission_height(P, ncomp, iseed, hmin, hmax):
    """Function to determine emission heights given the period. If no emission height range
       is specified, default used for P < 0.15 is between [950, 1000] and between [20, 1000] 
       for P > 0.15.
    
       Args:
       -----
       P      : Rotational period (seconds). 
       ncomp  : Integer number of component.
       iseed  : Integer seed for a pseudo-random number generator.
       hmin   : Minimum emission height (in km).
       hmax   : Maximum emission height (in km).

       
       Returns:
       --------
       H      : Emission heights (in km).
    """

    np.random.seed(iseed)
    num_H = ncomp # number of discrete emission height
    
#   If height range is not specified:
    if hmin == None and hmax == None:

        # emission height for a short period pulsar: only one emission height (KJ07) 
        if P <= 0.15:
            hmin = 950
            hmax = 1000
            H = np.random.uniform(hmin, hmax, size=1)

        elif P > 0.15:
            hmin = 20
            hmax = 500
            H = np.random.uniform(hmin, hmax, size=num_H)

#   For specified height range:
    else: H = np.random.uniform(hmin, hmax, size=num_H)

    return H
#======================
# Frequency dependence:
#======================
def height_f(H, freq):
    """Function to determine frequency dependent emission heights.

       Args:
       -----
       H     : heights (km)
       freq  : frequency (GHz)

       Returns
       -------
       H_mu  : Frequency dependent height (km)
    """
    #gamma = 0.83 # with rho \prop mu^-0.43 (error +/- 0.06 ref: fig.12 Hassall et al. 2012.)

    #H_mu = 0.6*H * (freq)**(-gamma) + 0.4*H # frequency dependence on height (KJ07 eqn.4/beam code)
    H_mu = H * (9 * freq**(-0.95) + 41)/(9 + 41) 
    return H_mu
#====================================================================================================================================================

#====================================================================================================================================================
#							OPENING ANGLE
#====================================================================================================================================================
def rho(P, heights):
    """Function to determine the opening angle of the beam given the rotational period and emission height.
    
       Args:
       -----
       P         : Rotational period (seconds).
       heights   : Emission heights (km).
       
       Returns:
       --------
       rho       : The opening angle (degrees).
       
    """

#   opening angle (eqn 3.29, Lorimer and Kramer 2005):
    rho = np.rad2deg(np.sqrt((9 * np.pi * heights) / (2 * (constants.c / 1000) * P)))   
    
    return rho

#====================================================================================================================================================
#							PATCH WIDTH
#====================================================================================================================================================
def patch_width(P, heights):
    """Function to calculate the width of a patchy emission region 
       within a pulsar beam at a given height.
    
       Args:
       -----
       P             : rotational period (seconds).
       heights       : emission heights (km).
      
       
       Returns:
       --------
       patchwidths   : the width of the patchy emission region (degrees).
       
    """   
    
#   width of the patch (eqn 3, KJ2007):
    patchwidths = 2.45 * 0.2 * np.sqrt(heights / ( 10 * P))    

    return patchwidths

#====================================================================================================================================================
# 							PATCH CENTER:
#====================================================================================================================================================
def patch_center(P, heights, npatch):
    """Function find centres of the patches
       
       Args:
       -----
       P       : rotatinal period
       heights : emission heights (in km).
       
       
       Returns:
       --------
       centerx : the patch center projection on the x-axis 
       centery : the patch center projection on the y-axis 
    """
    
#   opening angle:    
    opa = rho(P, heights) 
    centerx = []
    centery = []
    np.random.seed(iseed)
#   Fan beam model:
    if fanBeam == None:
        theta = 2 * np.pi * np.random.random(len(heights) * npatch)
    else:
        theta = 2 * np.pi * np.random.random(npatch)
    for j in range(len(heights)): #for each emission height (comp!)
#       find the center of the patch
        tempCenterX = []
        tempCenterY = []
        
        if fanBeam == None:
            for i in np.arange(npatch):
                tempCenterX.append(opa[j] * np.sin(theta[j*npatch + i]))
                tempCenterY.append(opa[j] * np.cos(theta[j*npatch + i]))
                
        else:
            for i in np.arange(npatch):
                tempCenterX.append(opa[j] * np.sin(theta[i]))
                tempCenterY.append(opa[j] * np.cos(theta[i]))


        centerx.append(tempCenterX)
        centery.append(tempCenterY)
        
    return centerx, centery


#====================================================================================================================================================
#							    POLARIZATION:
#====================================================================================================================================================
def rvm(alpha, beta, prof):
    """Function to determine polarization swing.
      
       Args:
       -----
       alpha   : inclination angle (degrees)
       beta    : impact parameter (degrees) 
       prof    : one dimentional profile

       Return:
       -------
       pa      : polarization position angle (degrees).

    """
    
 #  predict the position angle (psi) swing through the observer's sight line:
    zeta = alpha + beta
    points = range(len(prof))

    pa = []
    phi0 = 0.
    for point in points:
        phi = np.deg2rad(prof[point])
        numer = np.sin(np.deg2rad(alpha)) * np.sin(np.deg2rad(phi - phi0))
        denom = np.sin(np.deg2rad(zeta)) * np.cos(np.deg2rad(alpha)) - np.cos(np.rad2deg(zeta)) * np.sin(np.rad2deg(alpha)) * np.cos(np.deg2rad(phi - phi0))

        psi = np.arctan(numer/denom)

        # restrict psi between [-pi/2, pi/2]
        if psi < -np.pi/2.:
            psi = psi + np.pi

        if psi > np.pi/2.:
              psi = psi - np.pi

        # Convert psi back to degrees
        pa.append(np.rad2deg(psi))

       # phi0 = phi0 + 1
#    print "pa", pa
    return pa

#====================================================================================================================================================
#							ABERRATION:
#====================================================================================================================================================
def aberration(heights):
    """Function to determine the aberration ofset due to the curvature of the magnetic field.
 
       Args:
       -----
       heights     : emission heights 

       Returns:
       --------
       ab_ofsetx   : ofset of the projected x coordinates
       ab_ofsety   : ofset of the projected y coordinates

    """
#   aberration time scale:
    ab_time = heights / (constants.c / 1e3)
    ab_deg = (ab_time / P) * 360
    ab_xofset, ab_yofset = mapphi(alpha, 0.0, ab_deg)
    
    #print "aberration angle: " ab_deg
    return ab_xofset, ab_yofset

#====================================================================================================================================================
#							SCATTERING:
#====================================================================================================================================================
def sc_time(freq, dm, iseed):
    """Function to determine the scattering time scale as in Bhat et al. (2004).
               
       Args:
       -----
       freq   :   frequency (in GHz) 
       dm     :   dispersion measure (pc cm^-3)
    
       Return:
       -------
       tau     :   the scattering time (in sec)
       
    """
#   tau = scattering time scale as in Bhat et al. (2004)
    np.random.seed(iseed)     
    log_tau = -6.46 + 0.154 * np.log10(dm) + 1.07 * (np.log10(dm))**2 - 3.86 * np.log10(freq) +  np.random.uniform(-1,1) # scattering time with added noise term
    tau = 10**log_tau / 1e3 # (time scale in seconds)
    
    return tau

def pulsetrain(npulses, numberofbins, prof):
    """Function to create a train of pulses given a single pulse profile.

       Args:
       -----
       npulses        : number of pulses
       numberofbins   : number of bins (res; def = 1e3)
       prof           : pulse profile

       Return:
       -------
       train          : a train of pulse profiles (size = size(profile)*npulses)
    """

    binsrange = np.linspace(1, numberofbins, num=numberofbins, endpoint=True)
    nbins = np.max(binsrange)
    train = np.zeros(npulses * int(nbins))
    for i in range(npulses):
        startbin = i * nbins
        train[startbin:startbin + nbins] = prof
    
    return train


def extractpulse(sc_train, pulsesfromend, binsperpulse):
    """Function that takes the output convolution
       
       Args:
       -----
       sc_train     : a scattered train of pulse profiles 
       pulsefromend : number position of pulse to extract (from the last pulse)
       binsperpulse : number of bins per pulse

       Returns:
       --------
       pulse        : a single pulse profile
    """
    
    if pulsesfromend == 0:
        start = 0
        end = binsperpulse
        #zerobpulse = train[start:int(end)] - np.min(train[start:int(end)])
        #rectangle = np.min(train[start:int(end)])*binsperpulse
        #flux = np.sum(train[start:int(end)]) - rectangle

    else:
        start = -pulsesfromend*binsperpulse
        end = start + binsperpulse
        #zerobpulse = train[start:int(end)]-np.min(train[start:int(end)])
        #rectangle = np.min(train[start:inte(end)])*binsperpulse
        #flux = np.sum(train[start:int(end)]) - rectangle

    pulse = sc_train[int(start):int(end)]

    return pulse


def broadening(tau, P, res):
    """Function to determine the broadening function of the profile due to scattering.
       
       Args:
       -----
       tau         : scattering time (in seconds)
       P           : period (in seconds)

       Return:
       -------
       broad_func  : broadening function
    """
    t = np.linspace(0, P, num=res, endpoint=True)
    broad_func = 1/tau * np.exp(-(t / tau))
    #print "tau" + str(tau)
    return broad_func

def scatter(train, bf):
    """Function to scatter a pulse profile / a train of pulse profiles. Returns a convolution of the profile with the scattering function.

       Args:
       -----
       pulse   : pulse profile (extracted from a function extractpulse())
       bf      : broadening function

       Returns:
       -------
       conv    : scattered profile 
    """
    conv = np.convolve(train, bf)
    # normalise the profile:
    profint = np.sum(train) # integral / area of the profile 
    convint = np.sum(conv) # integral / area of the scattered profile
    sc_prof = conv * (profint / convint)
    out = sc_prof[0 : len(train) + 1] 

    return out


def find_peak(prof):
    """Function that finds a maximum peak of a profile.
       
       Args:
       -----
       prof  : pulse profile

       Return:
       -------
       peak  : peak value of the profile
    """
    
    peak = np.max(prof)
    
    return peak

def find_width(profile):
    peak = find_peak(profile)
    left = 0
    right = 0
    for i in range(len(profile)):
        if profile[i] > 0.1 * peak:
            left = i
            break
    for i in range(len(profile)):
        if profile[-i] > 0.1 *peak:
            right = len(profile) - i
            break
    w10 =(float(right -left)/float(len(profile)))* 360.
    return w10
 
#====================================================================================================================================================
#						   		ADD NOISE:
#====================================================================================================================================================
def noise_rms(snr, peak):
    """Function to determine the noise level given a signal to noise
       and the peak of the profile. Detemines the rms that will give maximum
       signal to noise. Assumes 0 baseline.
       
       Args:
       -----
       snr   : signal to noise ratio 
       peak  : peak of the profile

       Return:
       -------
       rms   : noise rms
    """
    
    rms = peak / snr

    return rms

def signal_to_noise(peak, rms):
    """Function to determine signal to noise ratio for each profile.
       Uses the previously determined noise level from function 
       "noise_rms" to determine the signal to noise ratio for a number 
       of profiles.

       Args:
       ---- 
       peak   : peak of the profile
       rms    : noise rms 

       Return:
       -------
       snr    : signal to noise of each profile
    """
    snr = (peak / rms )
    
    return snr


def add_noise(prof, rms, iseed, res):
    """Function that add noise to a profile. Finds a signal to noise of 
       a profile and determine the noise level to add for that specific
       profile. Adds a gaussian noise to a profile using a fixed noise
       rms, assuming a normalised profile.

       NB: If 'prof' is the scattered profile from scatterering function
       'scatter()', then the profile is normalised!

       Args:
       -----
       prof       : pulse profile
       rms        : noise rms
       iseed      : seed for random number generator 
       
       Returns:
       --------
       noisy_prof : a profile with added noise

    """
    peak = find_peak(prof)
    noise = np.random.normal(0, rms, res)
    noisy_prof = np.asarray(prof).T
    for i in range(np.shape(prof)[0]):
        noisy_prof[:,i] += noise
    '''
    plt.plot(noisy_prof)
    plt.show()
    '''
    noisy_prof = noisy_prof.T
    return noisy_prof

#====================================================================================================================================================
# 							BEAM PLOT:
#====================================================================================================================================================
def generateBeam(P, alpha, beta, freq, dm, heights, npatch, snr, do_ab):
    """Function to plot the patches for a given rotation period.
    
       A rgs:
       -----
       P       : rotational period (seconds)
       alpha   : inclination angle (degrees)
       beta    : impact parameter (degrees)
       heights : emission heights (in km)
       centerx : the patch center projection on the x-axis 
       centery : the patch center projection on the y-axis
       snr     : signal to noise ratio       
       Returns:
       --------
       prof    : an array containing 1D profile
       Z       : a 2D beam 
    
    """    
    
#   initialize parameters:
    xmin = -180.
    xmax = 180.
    res = 1e3 #resolution
    ymin = -180.
    ymax = 180.
    dx = (xmax - xmin)/res
    dy = (ymax - ymin)/res
    x = np.linspace(xmin, xmax, num=res, endpoint=True)
    y = np.linspace(ymin, ymax, num=res, endpoint=True)
    X,Y = np.meshgrid(x,y)
    gauss = np.zeros(len(x))
    Z = np.zeros_like(X)
    prof_los = np.zeros_like(gauss)

#   find the width of the patches:
    patchwidths = patch_width(P, heights)
#   An arbitrary peak of the profile:
    peakAmp = 1.   
#   Get the line of sight:
    xlos, ylos, thetalos = los(alpha, beta, res)
#   Get the centre of the emission patches on the xy-plane
    centerx, centery = patch_center(P, heights, npatch)
#   Get the ofset due to abberation:
    ab_xofset, ab_yofset = aberration(heights)

#   Find the 1D and 2D profile:
    for cid, comp in enumerate(heights):
#       widths for circular patches:        
        sigmax = patchwidths[cid]
        sigmay = patchwidths[cid]
    
#       center of the patch:
        patchCenterX = centerx[cid]
        patchCenterY = centery[cid]
        
#       2D patch (including aberation):
        for pc in zip(patchCenterX, patchCenterY):
            #if distance of current box from the patch center is
            #larger than 3 times the patchwidth, I do not want any
            #power contributed
            # first compute a grid of distance from the center of the patch, in units of the width
#            print "Patches: ", pc[0], pc[1], sigmax
            distance = (np.sqrt((X - pc[0])**2 + (Y - pc[1])**2))/sigmax
            distance[np.where(distance > 3.0)] = 0.0
            distance[np.where(distance != 0.0)] = peakAmp
            if do_ab == None:
                Z += distance * np.exp(-((X - pc[0])**2 / (2 * sigmax**2) + (Y - pc[1])**2 / (2 * sigmay**2)))
            else:
                Z += distance * np.exp(-((X - pc[0] - ab_xofset[cid])**2 / (2 * sigmax**2) + (Y - pc[1] - ab_yofset[cid])**2 / (2 * sigmay**2)))
            
#   1D profile from 2D patch, closest to the line of sight (select nearest neighbors):
    ZxIdx = np.array((xlos-xmin)/dx, dtype=int) # x index
    ZyIdx = np.array((ylos-ymin)/dy, dtype=int) # y index
    prof = Z[ZxIdx, ZyIdx]

    return prof, Z 

#                                            =============================================
#                                                      Command Line Parser
#                                            =============================================
#==================================================================================================================================================
#			    		      Initialise parameters from the command line
#==================================================================================================================================================

parser = argparse.ArgumentParser(description='Plot the patchy emission region as well as the line of sight profile. Running the file without specified argument will produce an output beam and profile from default parameters.')
parser.add_argument('-alpha', metavar="<alpha>", type=float, default='45', help='inclination angle in degrees (default = 45)')
parser.add_argument('-beta', metavar="<beta>", type=float, default='5', help='impact parameter in degrees (default = 5)')
parser.add_argument('-p', metavar="<p>", type=float, default='0.16', help='period in seconds (default = 0.16 s)')
parser.add_argument('-hmin', metavar="<hmin>", type=float, default=None, help='minimum emission height in km (default = {20 km for P > 0.15 s}, and {950 km for P < 0.15 s})')
parser.add_argument('-hmax', metavar="<hmax>", type=float, default=None, help='maximum emission height in km (default = 1000 km)')
parser.add_argument('-nc', metavar="<ncomp>", type=int, default='4', help='integer number of components (default = 4)')
parser.add_argument('-npatch', metavar="<npatch>", type=int, default='10', help='number of emission patches (default= 10)' )
parser.add_argument('-min_freq', metavar="<minfreq>", type=float, default='0.2', help='min frequency in GHz (default = 0.2 GHz)')
parser.add_argument('-chbw', metavar="<chanbw>", type=float, default='0.8', help='channel bandwidth in GHz (default = 0.8 GHz)')
parser.add_argument('-nch', metavar="<nch>", type=int, default='5', help='number of channels (default = 5)')
parser.add_argument('-iseed', metavar="<iseed>", type=int, default='4', help='integer seed for a pseudo-random number generator (default = 4)')
parser.add_argument('-snr', metavar="<snr>", type=float, default=None, help='signal to noise ratio (default = None)')
parser.add_argument('-dm', metavar="<dm>", type=float, default=1, help='dispersion measure in cm^-3 pc (default = 1)')
parser.add_argument('-outfile', metavar="<output file>", help="Write to file.")
parser.add_argument('-do_ab', default=None, help='include aberration ofset (default = None)')
parser.add_argument('-scatter', default=None, help='include scattering (default = None)')
parser.add_argument('-doFan', default=None, help='Fan beam - default: patchy beam')
args = parser.parse_args()
P = args.p
ncomp = args.nc
npatch = args.npatch
iseed = args.iseed
hmin = args.hmin
hmax = args.hmax
alpha = args.alpha
beta = args.beta
snr = args.snr
dm = args.dm
do_ab = args.do_ab
nch = args.nch
min_freq = args.min_freq
chbw = args.chbw
scr = args.scatter
fanBeam = args.doFan
fileName = args.outfile
beamFile = fileName+'_profiles'
paramsFile = fileName+'_parameters'
#====================================================================================================================================================
#				                         MAIN BODY: 
#====================================================================================================================================================
#=====================
# initialize params:
#=====================
beam = []
prof = []
peaks = []
w10 = []
res = 1e3
t_res = P/res
phase = np.linspace(-180, 180, num=res)
max_freq = (nch - 1) * chbw + min_freq
freq = np.linspace(min_freq, max_freq, nch) #channel frequency in GHz!!!

#=======================================
#     1. Find the emission height:
#=======================================
H = emission_height(P, ncomp, iseed, hmin, hmax)
#========================================
#     2. Get profile at each frequency:
#========================================
for i in np.arange(len(freq)):
    heights = height_f(H, freq[i]) # frequency dependent H
    pr, Z = generateBeam(P, alpha, beta, freq[i], dm, heights, npatch, snr, do_ab)
    w10.append(find_width(pr))
    prof.append(pr)
    beam.append(Z)


#==========================================
#     4. Scatter the line of sight profile: 
#==========================================
train = []
bf = []
tau = sc_time(freq, dm, iseed)
if scr == None:
    sc_prof = prof # returns the profile without scattering and store that in sc_prof

else:
    sc_prof = []
    for pid in np.arange(len(prof)):
        train.append(pulsetrain(3, res, prof[pid]))

    for fid in np.arange(len(freq)):
        tau = sc_time(freq[fid], dm, iseed)
        bf = broadening(tau, P, res)
        sc_train = scatter(train[fid], bf)
        sc_prof.append(extractpulse(sc_train, 2, res)) #scattered pulse profile

#===========================================
#     5. Add noise
#===========================================
for j in np.arange(len(prof)):
    peaks.append(find_peak(sc_prof[j]))
 
if snr == None:
    profile = sc_prof
else:
    rms = noise_rms(snr, np.max(peaks)) 
    profile = add_noise(sc_prof, rms, iseed, res) 

#========================================
#2.1 Write out the profile into a file:
#=======================================
# pickle "prof" using protocol 0 (ASCII)
profFile = open(beamFile, 'wb')
#profFile.write(' '.join([str(item) for item in prof]) + ' ')
pickle.dump(profile, profFile)
pickle.dump(phase, profFile)
pickle.dump(beam, profFile)
profFile.close()

pulsarParams = np.asarray([P, alpha, beta, w10[0], w10[-1], iseed])
f = open(fileName, 'a')
f.write(' '.join([str(item) for item in pulsarParams]) + ' ')
f.close()

params = open(paramsFile, 'wb')
pickle.dump(freq, params)
pickle.dump(P, params)
pickle.dump(iseed, params)
pickle.dump(res, params)
pickle.dump(t_res, params)
params.close()

#===========================================
# SET A ZERO BASELINE AND PLOT THE PROFILE:
#===========================================
fig = plt.figure()
for k in np.arange(len(profile)):
    plt.plot(phase, profile[k] + k, label='frequency = %0.2f GHz' %freq[k])
plt.title('A sequence of %i pulse profile' %nch)
plt.xlabel('phase (degrees)')
plt.xlim(-180,180)
plt.grid()

#    2D emission region:
meanBeam = np.mean(beam, axis=0)
xlos, ylos, thetalos = los(alpha, beta, res)
#for i in range(len(beam)):
'''
=======
>>>>>>> 5ff6357
plt.figure(figsize=(10,5))
plt.subplot(1, 2, 1)
plt.plot(xlos, ylos, '+r')
plt.imshow(beam[0].T, extent=[-180, 180, 180, -180])
plt.xlabel('X (degrees)')
plt.ylabel('Y (degrees)')
# find zoomed extent for plot
nonZero = np.where(beam[0]!= 0.0)
nzx_min = np.min(np.amin(nonZero,0))
nzy_min = np.min(np.amin(nonZero,1))
nzx_max = np.max(np.amax(nonZero,0))    
nzy_max = np.max(np.amax(nonZero,1))    
x1 = phase[nzx_min]
x2 = phase[nzx_max]
y1 = phase[nzy_min]
y2 = phase[nzy_max]
plt.xlim(x1,x2)
plt.ylim(y1,y2)
plt.colorbar()
plt.subplot(1, 2, 2)
<<<<<<< HEAD
#plt.plot(phase, averageP[0]) # average profile using first DM?
plt.plot(phase, profile[0]) # average profile using first DM?
plt.xlim(-180, 180)#lt.title('Profile at freq = %.4f GHz' % freq[-1])
'''

#==========================================
#   Scatter the profile:

#plt.figure(figsize=(10,5))
#plt.subplot(1, 2, 2)
#for pid in np.arange(len(sc_prof)):
#    plt.plot(phase, sc_prof[pid] + pid)
#plt.plot(phase, avg_prof(sc_prof))
#plt.title("Scattered Emission profile from LOS")
#plt.xlabel("phase")
#plt.ylabel("intensity")
#plt.xlim(-180, 180)
#plt.tight_layout()

#   patchy emission region:

#plt.subplot(1, 2, 1)
#plt.plot(xlos, ylos, '--r')
#plt.imshow(Z, extent=[-np.amax(Z),np.amax(Z),-np.amax(Z),np.amax(Z)])#, cmap=cm.gray)
#plt.title('Scattered Patchy emission region')
#plt.xlabel('X (degrees)')
#plt.ylabel('Y (degress)')
#plt.colorbar()
#==========================================


#   profile: with scattering
#plt.figure()
#for spid in np.arange(len(sc_prof)):
#    plt.plot(phase, sc_prof[spid] + spid)
    #plt.plot(phase, sc_prof[spid])
#plt.title("Scattered Emission profile from LOS")
#plt.xlabel("phase")
#plt.ylabel("intensity")
#plt.xlim(-180, 180)
#plt.tight_layout()
#plt.plot(x, prof_i)
#plt.figure()
#plt.plot(phase, avg_prof(sc_prof))
#   savefigure:

#================================
#  TEST THE DM_RANGE FUNCTION:
#===============================
#print "==================================="
#print "TESTING DM TRIAL...."
#print "==================================="
#for n in np.arange(len(prof)):
#    b,p = dm_trial(prof[n])
#    print "peaks", p
#    print "bin", b
#plt.show()
plt.plot(phase, profile[0]) 
plt.xlim(-180, 180)
plt.xlabel('Phase (degrees)')

plt.show()
