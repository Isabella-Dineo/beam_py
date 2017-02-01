#!/usr/bin/env python

# A module containing function used to generate pulsar beam and profile.

import numpy as np
from scipy import constants
import scipy.stats as stats

#====================================================================================================================================================
#                                                       PRECISION ERRORS                
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
#                                                       ROTATIONAL AXES         
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
#                                                       LINE OF SIGHT
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
#                                                       EMISSION HEIGHTS
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
#                                                       OPENING ANGLE
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
#                                                       PATCH WIDTH
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
#                                                       PATCH CENTER:
#====================================================================================================================================================
def patch_center(P, heights, npatch, iseed, fanBeam):
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
#                                                           POLARIZATION:
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

    return pa

#====================================================================================================================================================
#                                                       ABERRATION:
#====================================================================================================================================================
def aberration(heights, P, alpha):
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

    return ab_xofset, ab_yofset

#====================================================================================================================================================
#                                                       SCATTERING:
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
        train[int(startbin):int(startbin + nbins)] = prof

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

def find_width(prof):
    """ Function to find the width at 10% of the peak intensity.

       Args:
       -----
       prof    : profile (1D array)

       Returns:
       --------
       w10     : width at 10% of the peak intensity
    """

    peak = find_peak(prof)
    left = 0
    right = 0
    for i in range(len(prof)):
        if prof[i] > 0.1 * peak:
            left = i
            break
    for i in range(len(prof)):
        if prof[-i] > 0.1 *peak:
            right = len(prof) - i
            break
    w10 =(float(right -left)/float(len(prof)))* 360.
    return w10

def getadm(psrcatdm, iseed, nbins, n):
    """Function to randomly select a dm value from the known pulsar dms in the catalogue.

       Creates a distribution of the dm values given their probabilities, and randomly select
       a dm value to use for scattering.
       
       Args:
       -----
       psrcatdmfile : a file containing psrcat dm values in 1 column (nan values replaced with zeros)
       iseed        : seed for the random number generator [int].
       nbins        : number of bins.
       n            : size of the samples to draw
       
       Returns:
       --------
       rand_dm      : randomly selected dm value (pc cm^-3). 
    """
    dm_file_name = str(psrcatdm)
    dm_file = np.loadtxt(dm_file_name)                  # Load the txt file containing the DM
    dm_dat = dm_file[np.where(dm_file > 0)]             # Exclude the zero dms used to replace null values from psrcat
    hist, bin_edges = np.histogram(dm_dat, bins=nbins)  # creates a histogram distribution
    probs = hist/float(len(dm_dat))                     # Compute probabilities
    dm_range = np.linspace(np.min(dm_dat), np.max(dm_dat), endpoint=True, num=len(probs)) 
    normdiscrete = stats.rv_discrete(values=(dm_range, probs), seed=iseed) # Find an arbitrary distribution
    rand_dm = normdiscrete.rvs(size=n)                         # draw a sample of size n
    
    return rand_dm
#====================================================================================================================================================
#                                                               ADD NOISE:
#====================================================================================================================================================
def noise_rms(snr):
    #noise_rms(snr, peak):
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
    # assumes the beam is normalized such that the LOS cut at center produce a beam (2d gaussian with max amplitude 1)
    # The peak of the profile = 10 times    
    #rms = peak / snr
    rms = 10/snr  
    
    return rms

def add_noise(prof, rms, res):
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
       
       Returns:
       --------
       noisy_prof : a profile with added noise

    """
    peak = find_peak(prof)
    noise = np.random.normal(0, rms, int(res))
    noisy_prof = np.asarray(prof).T
    for i in range(np.shape(prof)[0]):
        noisy_prof[:,i] += noise
    noisy_prof = noisy_prof.T
    return noisy_prof


def signal_to_noise(peak, rms):
    """Function to determine signal to noise ratio for each profile.
       Uses the previously determined noise level from function 
       'noise_rms' to determine the signal to noise ratio for each 
       profile.

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
#====================================================================================================================================================
#                                                       DM CURVE FITTING:
#====================================================================================================================================================
def find_phase_bin(prof):
    """Function to find a bin closest to the peak of a profile.
       
       Args:
       -----
       profile   : a profile (1D numpy array)

       Returns:
       --------
       phase_bin : a bin where the profile peaks

    """
    # Find the bin where the profile is max:

    peak = np.max(prof)
    for prof_id, prof_val in enumerate(prof):
        if prof[prof_id] == peak:
            phase_bin = prof_id

    return phase_bin


# Find a phase/time corresponding to the peak of the profile
def find_delta_dm(P, prof, phase, phase_bin0, phase_bin1, freq_ref, freq, nch):
    """Function to determine a range of dispersion measures 
       
       Args:
       -----
       P          : rotation period (seconds)
       prof       : profile (numpy array)
       phase      : rotation phase (numpy array)
       freq_ref   : reference frequecy (in GHz)
       freq       : frequency to shift (in GHz)
       phase_bin0 : phase bin corresponding to the peak of the profile at max frequency (or reference frequency) (integer)
       phase_bin1 : phase bin corresponding to the peak of the profile at min frequency (integer)

       Returns:
       --------
       delta_dm   : a range of dispersion measures to try (numpy array)
       
    """
    # Find the delta phase shift between the min and max frequency profile:
    phase_at_peak0 = phase[phase_bin0] # peak at reference profile (/max freq)
    phase_at_peak1 = phase[phase_bin1]
    D = 4.148808 * 1e3 # +/- 3e-6 MHz^2 pc^-1 cm^3 s
    range_dm = 0.01 # dm units
    if phase_bin0==phase_bin1:
        dm = 0
    else:
        delta_phase = phase_at_peak1 - phase_at_peak0
        # Convert the phase to time and find a corresponding delta_dm:
        delta_t = delta_phase/360. * P
        dm = delta_t / (D * ((freq_ref * 1e3)**(-2) - (freq * 1e3)**(-2)))
#    delta_dm = np.linspace(-dm - range_dm, -dm + range_dm, num=20)# try only 20 for now
    delta_dm = np.linspace(- 5. * dm, 3. * dm , num=20)# try only 20 for now
    return delta_dm


def delay(freq_ref, freq , delta_dm, t_res):
    """Function to determine the delay of the profiles as a function of freq.
       Assumes the the profiles are already de-despersed; this is an additional
       delay due to the dm variation with profile
       
       Args:
       -----
       freq_ref : reference frequecy (in GHz)
       freq     : frequency to shift (in GHz)
       delta_dm : dispersion measures to try
       t_res    : time per bin (Period/resolution in seconds)

       Return:
       -------
       bin_shift: relative shift of freq w.r.t the reference freq (in bins)
       
       Function use the maximum frequency as the reference frequency. 
       Find the dispersion measure that give the maximum S/N.

    """
    D = 4.148808 * 1e3 # +/- 3e-6 MHz^2 pc^-1 cm^3 s
    delta_t = D * ((freq_ref * 1e3)**(-2) - (freq * 1e3)**(-2)) * delta_dm # in seconds; eqn. 5.2 (handbook)
    if np.isnan(delta_t):
        delta_t = np.nan_to_num(delta_t)
    else:
        delta_t = delta_t

    bin_shift = np.int(delta_t * (1/t_res))    # phase shift in bins  

    return bin_shift


def avg_prof(prof):
    """Function to average profiles
    """
    profile = np.asarray(prof)
    averageP = np.average(prof, 0)

    return averageP

