#!/usr/bin/env python
#plotbeam.py , last mod 03/02/2016 [IDR]

import numpy as np 
import matplotlib.pyplot as plt 
import matplotlib.cm as cm
from scipy import constants
import argparse

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
 
    return np.rad2deg(xp), np.rad2deg(yp)
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
    for i in np.arange(len(thetalos)):
        if thetalos[i] < 0:
            thetalos[i] = -thetalos[i]       
            
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
            hmax = 1000
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
    gamma = 0.86 # with rho \prop mu^-0.43 (error +/- 0.06 ref: fig.12 Hassall et al. 2012.)

    H_mu = 0.6*H * freq**(-gamma) + 0.4*H # frequency dependence on height (KJ07 eqn.4/beam code)
    
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
    
#   initialize the array:
    centerx = []
    centery = []
    np.random.seed(iseed)
    npatch = np.random.randint(2,10+1)
    
    for comp in opa: #for each emission height (comp!)
#       find the center of the patch
        tempCenterX = []
        tempCenterY = []
        theta = 2 * np.pi * np.random.random(npatch)

        for i in np.arange(npatch):
            tempCenterX.append(comp * np.sin(theta[i]))
            tempCenterY.append(comp * np.cos(theta[i]))
                

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
#							SCATTERING TIME:
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
#    print "tau: %f sec" %tau
    
    return tau

#====================================================================================================================================================
#                					PERIODIC PROFILE:
#====================================================================================================================================================
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

#===================================================================================================================================================
#							EXTRACT A PULSE
#===================================================================================================================================================

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
#====================================================================================================================================================
#						     BROADENING FUNCTION:
#====================================================================================================================================================
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

#====================================================================================================================================================
#						        SCATTERING:
#====================================================================================================================================================
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

#===========================================================================================================================================================
def dm_trial(prof):
    """Function to find dm range for dm trial.
       
       Args:
       -----
       width     : pulse width at 10% of the peak( w10 in degrees)
       freq0     : reference frequency (GHz)
       freq      : frequency to shift (GHz) 

       Returns:
       --------
       delta_dm  : differential dispersion measure (pc cm^-3).

    """
    # Find the bin where the profile is max:
    peak = np.max(prof)
    for prof_id, prof_val in enumerate(prof):
        if prof[prof_id] == peak:
            phase_bin = prof_id
    # Find the width at 10% of the peak
    peak10 = peak * 0.1 #10% of the peak
    # width 

# Find the corresponding time shift:
#    tim = phase_at_peak/360 * P
#    tim_bin = int(tim * (1/t_res)) # index
    

 #   D = 4.148808 * 1e3 # +/- 3e-6 MHz^2 pc^-1 cm^3 s
 #   delta_t = phase_at_peak * P  # time in seconds
 #   delta_dm = delta_t / (D * (freq0**(-2) - freq**(-2)))
 #   dm_range = np.linpsace(-delta_dm, delta_dm, endpoint=True, num=10) # try only 10 dm trials (FOR NOW!!!!)

    return phase_bin, peak

#============================================================================================================================================================
#						DELAY
#============================================================================================================================================================
def delay(freq0, freq , delta_dm, t_res):
    """Function to determine the delay of the profiles as a function of freq.
       Assumes the the profiles are already de-despersed; this is an additional
       delay due to the dm variation with profile
       
       Args:
       -----
       freq0    : reference frequecy (in GHz)
       freq     : frequency to shift (in GHz)
       delta_dm : dispersion measures to try
       t_res    : time per bin (Period/resolution in seconds)

       Return:
       -------
       delta_t  : dispersive delay (in seconds)
       
       Function use the maximum frequency as the reference frequency. 
       Find the dispersion measure that give the maximum S/N.

    """
    D = 4.148808 * 1e3 # +/- 3e-6 MHz^2 pc^-1 cm^3 s
    delta_t = D * ((freq0 * 1e3)**(-2) - (freq * 1e3)**(-2)) * delta_dm # in seconds; eqn. 5.2 (handbook)
    bin_shift = np.int(delta_t * (1/t_res))    # phase shift in bins  

    return bin_shift
#====================================================================================================================================================
# 							AVERAGE PROFILE:
#====================================================================================================================================================
def avg_prof(prof):
    """Function to average profiles
    """
    profile = np.asarray(prof)
    averageP = np.average(prof, 0)

    return averageP
#====================================================================================================================================================
#							FIND PROFILE PEAK:
#====================================================================================================================================================
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

#====================================================================================================================================================
#						   DETERMINE THE NOISE LEVEL:
#====================================================================================================================================================
def noise_rms(snr, peak):
    """Function to determine the noise level given a signal to noise
       and the peak of the profile. Detemines the rms that will give maximum
       signal to noise
       
       Args:
       -----
       snr   : signal to noise ratio 
       peak  : peak of the profile

       Return:
       -------
       rms   : noise rms
    """
    
    rms = np.sqrt(peak**2 / snr)

    return rms

#====================================================================================================================================================
#						     SIGNAL TO NOISE:
#====================================================================================================================================================
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
#====================================================================================================================================================
#							ADD NOISE:
#====================================================================================================================================================
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
    noisy_prof = prof + noise
    
    return noisy_prof

#====================================================================================================================================================
# 							BEAM PLOT:
#====================================================================================================================================================
def plotpatch(P, alpha, beta, freq, dm, heights, npatch, snr, do_ab):
    """Function to plot the patches for a given rotation period.
    
       Args:
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
       A plot of the patches projected on to observational plane.
    
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
    peak = 10. 
#    peak = 1.   
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
            if do_ab == None:
                Z += peak * np.exp(-((X - pc[0])**2 / (2 * sigmax**2) + (Y - pc[1])**2 / (2 * sigmay**2)))
            else:
                Z += peak * np.exp(-((X - pc[0] - ab_xofset[cid])**2 / (2 * sigmax**2) + (Y - pc[1] - ab_yofset[cid])**2 / (2 * sigmay**2)))
            
#   1D profile from 2D patch, closest to the line of sight (select nearest neighbors):
    
    ZxIdx = np.array(xlos/dx, dtype=int) - int(res/2) # x index
    ZyIdx = np.array(ylos/dy, dtype=int) - int(res/2) # y index
    prof = Z[ZxIdx, ZyIdx]

#   Scattering:
    tau = sc_time(freq, dm, iseed)
    bf = broadening(tau, P, res)
    train =  pulsetrain(3, res, prof)
    sc_train = scatter(train, bf)
    sc_prof = extractpulse(sc_train, 2, res) #scattered pulse profile
    
#   patchy emission region:
    
#    plt.figure(figsize=(10,5))
#    plt.subplot(1, 2, 1)
#    plt.plot(xlos, ylos, '--r')
#    plt.imshow(Z, extent=[-np.amax(Z),np.amax(Z),-np.amax(Z),np.amax(Z)])#, cmap=cm.gray)
#    plt.title('Scattered Patchy emission region %.3f GHz' %freq)
#    plt.xlabel('X (degrees)')
#    plt.ylabel('Y (degress)')
#    plt.colorbar()


#   profile: with scattering
#    plt.subplot(1, 2, 2)
#    plt.title("Scattered Emission profile from LOS")
#    plt.xlabel("phase")
#    plt.ylabel("intensity")
#    plt.xlim(xmin, xmax)
#    plt.tight_layout()
    #plt.plot(x, prof_i)
    #plt.figure()
#    plt.plot(x, sc_prof)
#   savefigure:
#
#    plt.show()
    #time.ctime(time.time())
    #file_num = time.time()
 #   plt.savefig('beam_F_' + str(freq) + '_dm_' + str(dm) + '_' + str(file_num) + '.pdf', format='pdf')

    return prof, Z # the unscattered profile

#                                            =============================================
#                                                      Command Line Parser
#                                            =============================================
#==================================================================================================================================================
#			    		      Initialise parameters from the command line
#==================================================================================================================================================
# argparse is a complete argument processing library. Arguments can trigger different actions, specified by the "action" argument 
# to add_argument(). 
# Supported actions include storing the argument (singly, or as part of a list), storing a constant value when the argument is
# encountered (including special handling for true/false values for boolean switches), counting the number of times
# an argument is seen, and calling a callback. 
# The default action is to store the argument value. In this case, if a type is provided, the value is converted to that type 
# before it is stored. 
# If the dest argument is provided, the value is saved to an attribute of that name on the Namespace object returned when the 
# command line arguments are parsed.
#
# Once all of the arguments are defined, you can parse the command line by passing a sequence of argument strings to parse_args(). 
# By default, the arguments are taken from sys.argv[1:], but you can also pass your own list. The options are processed using the 
# GNU/POSIX syntax, so option and argument values can be mixed in the sequence.
# The return value from parse_args() is a Namespace containing the arguments to the command. The object holds the argument values 
# as attributes, so if your argument dest is "myoption", you access the value as args.myoption. 
# (source: https://pymotw.com/2/argparse/)

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
#parser.add_argument('-o','--outfile', metavar="<name_suffix>", dest='output', action='store', type=argparse.FileType('w'), help="Write to file.")
parser.add_argument('-do_ab', default=None, help='include aberration ofset (default = None)')
#parser.add_argument('-d', metavar="<dir>", default='/home/', help='Directory to save plots.')
#parser.add_argument('-x', "--x11", action='store_true', help='X11 window plot override switch.')

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
#output = args.output
#plotwindow = args.x11
#dir = args.dir
#====================================================================================================================================================
#				                   plot the beam 
#====================================================================================================================================================
#=====================
# initialize params:
#=====================
beam = []
prof = []
peaks = []
res = 1e3
t_res = P/res
phase = np.linspace(-180, 180, num=res)
max_freq = (nch - 1) * chbw + min_freq
freq = np.linspace(min_freq, max_freq, nch) #channel frequency in GHz!!!

#===========================
# Find the emission height:
#===========================
H = emission_height(P, ncomp, iseed, hmin, hmax)

#================================================
# Get profile for each frequency and their peaks:
#================================================
for i in np.arange(len(freq)):
    heights = height_f(H, freq[i]) # frequency dependent H
    print "heights : " + str(heights) + "km"
    pr, Z = plotpatch(P, alpha, beta, freq[i], dm, heights, npatch, snr, do_ab)
    prof.append(pr)
    beam.append(Z)
    print "freq: %f (GHz)" %freq[i]
#   find the phase shift per frequency channel:
#    phase_shift.append(delay(freq_c, chbw, dm, t_res))   

for j in np.arange(len(prof)):
    peaks.append(find_peak(prof[j]))


#========================================
# plot the two 2D profile with frequency:
#========================================
#plt.subplot(1, 2, 1)
plt.figure()
plt.title('Beam')
plt.xlabel('phase bin number')
plt.ylabel('channels')
plt.imshow(prof, aspect='auto', origin='lower')

# 2D emission region:
# Plot the pulse profiles for each frequency:
#plt.subplot(1, 2, 2)
#ax = fig.gca(projection='3d')
#z = np.arange(nch)


#shifted_phase = []
#for shif in np.arange(len(phase_shift)):
#    for phs in np.arange(len(phase)):
#        shifted_phase.append(phase[phs] - phase_shift[shif])
#print len(phase), len(phase_shift)


#ax = fig.add_subplot(111, projection='3d')
#phase = np.linspace(-180, 180, num=1e3)

#===================
# add noise:
#===================
profile = []
if snr == None:
    profile = prof
else:
    rms = noise_rms(snr, np.max(peaks)) 
    print "noise rms: %.4f" %rms
    profile = add_noise(prof, rms, iseed, res) 
    #print np.shape(profile), np.shape(prof)

#=======================================================
# Determine width of the profiles at 10% of the profile:
#=======================================================
# Find the peak of these profiles:
#peak_at_minf = peaks[0]
#peak_at_maxf = peaks[nch -1]

# Find the width of the profile at 10% of the peak:
#peaks10 = peaks * 10/100.
#peak_max10 = peak_at_maxf * 10/100.

# Find the width of the profile at 10% of the peak:
# (width using the "noise-less" profile prof; scattered profile sc_prof in due time):
#w_for_minf = [] # min and max phase for min_freq profile(to be used to find the width) 
#for pid1, p1 in enumerate(prof[0]):
#    if p1 >= int(peak_min10) - 10/100. and p1 <= int(peak_min10) + 10/100.:
#        w_for_minf.append(phase[pid1])
        #print "phase for min: ", phase[pid1]

#w_for_maxf = [] # min and max phase for max_freq profile:
#for pid2, p2 in enumerate(prof[nch - 1]):
#    if p2 >= int(peak_max10) - 10/100. and p2 <= int(peak_max10) + 10/100.:
#        w_for_maxf.append(phase[pid2])
        #print "phase for max: ", phase[pid2]
#width_at_min = w_for_minf[len(w_for_minf) - 1] - w_for_minf[0]
#width_at_max = w_for_maxf[len(w_for_maxf) - 1] - w_for_maxf[0]

#print "width at min_freq: ", width_at_min
#print "width at max_freq: ", width_at_max

#===============================================
# Find a range for dm trial and phase shift:
#===============================================
#delta_dm = [] 
#dm_range = [] # this will be an array(shape : (nch,10)) of arrays.
#phase_shift = []
#freq0 = max_freq # shift the profiles to max freq
#for frequency in np.arange(nch):
#    delta_dm.append(dm_trial(width_at_minf, freq0, freq[frequency]))
#    dm_range.append(np.linspace(-delta_dm, delta_dm, num=10)) # try 10 dm's (for now!)
#best_dm = find_dm(dm_range, prof)
#phase_shift.append(delay(freq0, freq, delta_dm))

# create an xy-plane:
#prof_num = np.arange(1,len(noisy_prof) + 1)
#Xp, Yp = np.meshgrid(phase, prof_num) 
#ax = fig.gca(projection = '3d')

#=========================================
# set zero a baseline and plot profile:
#=========================================
fig = plt.figure()
for k in np.arange(len(profile)):
    prof_num = np.arange(1, len(profile) + 1)
    #phase = np.linspace(-180, 180, num=1e3)
    profile[k] -= np.min(profile[k])
    #Xp, Yp, Zp = np.meshgrid(phase, profile[k], prof_num) 

    #ax.plot(phase, profile[k][k])
    plt.plot(phase, profile[k] + k, label='frequency = %0.2f GHz' %freq[k])
    #plt.plot(shifted_phase, noisy_prof[j] + j)

#for k in np.arange(len(prof)):
#    prof[j]-=np.min(prof[j]) 
    #ax.plot(phase, prof[j], j)
plt.title('A sequence of %i pulse profile' %nch)
plt.xlabel('phase (degrees)')
plt.xlim(-180, 180)
plt.ylabel('profile number')
plt.grid()

#====================
# AVERAGE PROFILE:
#====================
averageP = avg_prof(prof)

# 2D emission region:
xlos, ylos, thetalos = los(alpha, beta, res)
plt.figure(figsize=(10,5))
plt.subplot(1, 2, 1)
plt.plot(xlos, ylos, '--r')
plt.imshow(beam[0], extent=[-np.amax(beam[0]),np.amax(beam[0]),-np.amax(beam[0]),np.amax(beam[0])])#, cmap=cm.gray)
plt.title('Patchy emission')
plt.xlabel('X (degrees)')
plt.ylabel('Y (degress)')
plt.colorbar()
plt.subplot(1, 2, 2)
#plt.plot(phase, profile[0], label='freq : %.4f MHz' %min_freq)   # profile at min freq
#print len(profile), len(profile[0]) 
plt.plot(phase, averageP) # average profile
plt.legend()
plt.xlim(-180, 180)
plt.title('Average pulse profile')
#print dm, type(dm), freq, type(freq)
#print "delta t %i" %len(delta_t)
#print delta_t
#print "nchan %i" %nch
#plt.figure()
#plt.plot(np.linspace(0, nch, num=len(delta_t)), delta_t)

#=============================
#   Scattering: dublicate-ish
#=============================
train = []
bf = []
sc_prof = []
tau = sc_time(freq, dm, iseed)
for pid in np.arange(len(prof)):
    train.append(pulsetrain(3, res, prof[pid]))

plt.figure(figsize=(10,5))
plt.subplot(1, 2, 2)
for fid in np.arange(len(freq)):
    tau = sc_time(freq[fid], dm, iseed)
    bf = broadening(tau, P, res)
    sc_train = scatter(train[fid], bf)
    sc_prof.append(extractpulse(sc_train, 2, res)) #scattered pulse profile
    #print len(sc_prof), len(phase)
    #plt.plot(phase, sc_prof)
#print len(avg_prof(sc_prof))
plt.plot(phase, avg_prof(sc_prof))
plt.title("Scattered Emission profile from LOS")
plt.xlabel("phase")
plt.ylabel("intensity")
plt.xlim(-180, 180)
plt.tight_layout()
#plt.plot(x, prof_i)
#plt.figure()
#plt.plot(phase, avg_prof(sc_prof))

#   patchy emission region:

plt.subplot(1, 2, 1)
plt.plot(xlos, ylos, '--r')
plt.imshow(Z, extent=[-np.amax(Z),np.amax(Z),-np.amax(Z),np.amax(Z)])#, cmap=cm.gray)
plt.title('Scattered Patchy emission region')
plt.xlabel('X (degrees)')
plt.ylabel('Y (degress)')
plt.colorbar()

#   profile: with scattering
plt.figure()
for spid in np.arange(len(sc_prof)):
    plt.plot(phase, sc_prof[spid] + spid)
plt.title("Scattered Emission profile from LOS")
plt.xlabel("phase")
plt.ylabel("intensity")
plt.xlim(-180, 180)
plt.tight_layout()
#plt.plot(x, prof_i)
#plt.figure()
#plt.plot(phase, avg_prof(sc_prof))
#   savefigure:
plt.show()

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
