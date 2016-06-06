#!/usr/bin/env python
#plotbeam.py , last modified 25/05/2016 [ID]

import numpy as np 
import time
import matplotlib.pyplot as plt 
import matplotlib.cm as cm
from scipy import constants
import argparse

#                                                 ============================
#                                                       DEFINE FUNCTIONS:
#                                                 ============================
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
    """Function to project the beam onto the rotational axes.
       
       Args:
       -----
       alpha       : Magnetic inclination angle w.r.t the rotatinal axis (degrees).
       beta        : Line of sight closest approach to the magnetic axis (degrees).
       
       Returns:
       --------
       xlos, ylos  : The coordinates of the rotatinal plane (degrees).
       
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
            cosgamma[i] = (np.cos(np.deg2rad(alpha+beta)) - np.cos(np.deg2rad(alpha)) * cosR[i])/(np.sin(np.deg2rad(alpha)) * np.sin(R[i]))
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
def los(alpha, beta):
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
    phi = np.linspace(-180, 180, num=1e3, endpoint=True)
    
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
def emission_height(P, freq,  ncomp, iseed, hmin, hmax):
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
    num_H = ncomp # number of discrete emission heights
    
#   If height range is not specified:
    if hmin == None and hmax == None:
        
        # emission height for a short period pulsar: only one emission height 
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
             
        #if P > 0.15: 
        #    H = np.random.uniform(hmin, hmax, size=num_H) 

#   Frequency dependence:
    gamma = 0.86 # with rho \prop mu^-0.43 (error +/- 0.06 ref: fig.12 Hassall et al. 2012.)
      
    H_mu = 0.6*H * freq**(-gamma) + 0.4*H # frequency dependence on height (KJ07 eqn.4/beam code)
    #H_mu = H * freq**(-gamma) + H0
    print "H_mu",H_mu            
        
    return H

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
       heights : emission heights (in km)
       npatch  : integer number of emission patches
       
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
    #npatch = np.random.randint(2,10+1)
    
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
        
    #print "acive patches: " ,npatch
        
    return centerx, centery

#====================================================================================================================================================
#							    RVM:
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
    
    return ab_xofset, ab_yofset

#====================================================================================================================================================
#							SCATTERING TIME:
#====================================================================================================================================================
def sc_time(freq, dm):
    """Function to determine the scattering time scale as in Bhat et al. (2004).
               
       Args:
       -----
       freq   :   frequency (in GHz) (array)
       dm     :   dispersion measure (pc cm^-3).
    
       Return:
       -------
       tau     :   the scattering time in sec (array).
       
    """
#   tau = scattering time scale as in Bhat et al. (2004)     
    log_tau = -6.46 + 0.154 * np.log10(dm) + 1.07 * (np.log10(dm))**2 - 3.86 * np.log10(freq) * np.random.rand() # random time scale along ISM
    tau = 10**log_tau * 1e3 # (time scale in seconds)
    
    return tau

#====================================================================================================================================================
#						     BROADENING FUNCTION:
#====================================================================================================================================================
def broadening(tau, P):
    """Function to broaden the profile for scattering.
       
       Args:
       -----
       tau         : scattering time (in seconds)
       P           : period (in seconds)

       Return:
       -------
       broad_func : broadening function
    """
    t = np.linspace(0, 5*P, num=1e3, endpoint=True)
    broad_func = 1/tau * np.exp(-(t / tau))

    return broad_func

#====================================================================================================================================================
#						        SCATTERING:
#====================================================================================================================================================
def scatter(prof, bf):
    """Function to scatter a pulse profile. Returns a convolution of the profile with the scattering function.

       Args:
       -----
       prof    : profile
       bf      : broadening function

       Returns:
       -------
       conv    : scattered profile 
    """
    conv = np.convolve(prof, bf)
    # normalise the profile:
    profint = np.sum(prof) # integral / area of the profile 
    convint = np.sum(conv) # integral / area of the scattered profile
    sc_prof = conv * (profint / convint)
    out = sc_prof[0:1e3] 

    return out
#====================================================================================================================================================
# 							BEAM PLOT:
#====================================================================================================================================================
def plotpatch(P, alpha, beta, heights, centerx, centery, snr, do_ab):
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
    xnum = 1e3
    ymin = -180.
    ymax = 180.
    ynum = 1e3
    dx = (xmax - xmin)/xnum
    dy = (ymax - ymin)/ynum
    x = np.linspace(xmin, xmax, num=xnum, endpoint=True)
    y = np.linspace(ymin, ymax, num=ynum, endpoint=True)
    X,Y = np.meshgrid(x,y)
    gauss = np.zeros(len(x))
    Z = np.zeros_like(X)
    prof_los = np.zeros_like(gauss)

#   find the width of the patches:
    patchwidths = patch_width(P, heights)
      
#   An arbitrary peak of the profile:
    peak = 10. 
    
#   Get the line of sight:
    xlos, ylos, thetalos = los(alpha, beta)

#   Get the ofset due to abberation:
    ab_xofset, ab_yofset = aberration(heights)

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
    
    ZxIdx = np.array(xlos/dx, dtype=int) - int(xnum/2) # x index
    ZyIdx = np.array(ylos/dy, dtype=int) - int(ynum/2) # y index
    prof = Z[ZxIdx, ZyIdx]

#   Scattering:
#   1. Find the scattering time in seconds:
    tau = sc_time(freq, dm)
#   2. Determine the function that prodens the profile:
    bf = broadening(tau, P)
#   3. scatter the profile:
    sc_prof = scatter(prof, bf) # scattered profile.


#   Polarized emission:
    #fractional pol:
    #frac_c = peak * (np.random.random(1) - 0.5) / 2
#   polarization position angle
#   pa = rvm(alpha, beta, prof)

#   add noise: (gaussian noise)
    if snr == None:
        prof_i = prof # profile without noise
        
    else:
        sigma_s = sigmax #std dev of the profile
        sigma_n = sigmax/np.sqrt(snr) #std dev of the noise
        mean_n = (sigmax**2) * np.sqrt(snr)     
        noise = np.random.normal(mean_n, sigma_n, 1e3)
        prof_i = prof + noise
        #print sigma_n, mean_n
    
#   patchy emission region:
    
    plt.figure(figsize=(10,5))
    plt.subplot(1, 2, 1)
    plt.plot(xlos, ylos, '--r')
    plt.imshow(Z, extent=[-np.amax(Z),np.amax(Z),-np.amax(Z),np.amax(Z)])#, cmap=cm.gray)
    plt.title('Patchy emission region of the beam')
    plt.xlabel('X (degrees)')
    plt.ylabel('Y (degress)')
    plt.colorbar()


#   profile: with scattering
    plt.subplot(1, 2, 2)
    plt.title("Emission profile from LOS")
    plt.xlabel("phase")
    plt.ylabel("intensity")
    plt.xlim(xmin, xmax)
    plt.tight_layout()
#   plt.plot(x, prof_i)
#   plt.figure()
    plt.plot(x, sc_prof)
#   savefigure:
#
#    plt.savefig('', dpi=None, facecolor='w', edgecolor='w',
#        orientation='portrait', papertype=None, format='pdf',
#        transparent=False, bbox_inches=None, pad_inches=0.1,
#        frameon=None)
#    plt.show()
    #time.ctime(time.time())
    file_num = time.time()
    plt.savefig('beam' + str(file_num) + '.pdf', format='pdf')

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

# Once all of the arguments are defined, you can parse the command line by passing a sequence of argument strings to parse_args(). 
# By default, the arguments are taken from sys.argv[1:], but you can also pass your own list. The options are processed using the 
# GNU/POSIX syntax, so option and argument values can be mixed in the sequence.
# The return value from parse_args() is a Namespace containing the arguments to the command. The object holds the argument values 
# as attributes, so if your argument dest is "myoption", you access the value as args.myoption. 
# (source: https://pymotw.com/2/argparse/)

parser = argparse.ArgumentParser(description='Plot the patchy emission region as well as the line of sight profile. Running the file without specified argument will produce an output beam and profile from default parameters.')
parser.add_argument('--alpha', metavar="<alpha>", type=float, default='45', help='inclination angle in degrees (default = 45)')
parser.add_argument('--beta', metavar="<beta>", type=float, default='5', help='impact parameter in degrees (default = 5)')
parser.add_argument('--freq', metavar="<freq>", type=float, default='1.4', help='frequency in GHz (default = 1.4 GHz)')
parser.add_argument('-p', metavar="<p>", type=float, default='0.16', help='period in seconds (default = 0.16 s)')
parser.add_argument('-nc', metavar="<ncomp>", type=int, default='4', help='integer number of components (default = 4)')
parser.add_argument('--iseed', metavar="<iseed>", type=int, default='4', help='integer seed for a pseudo-random number generator (default = 4)')
parser.add_argument('--hmin', metavar="<hmin>", type=float, default=None, help='minimum emission height in km (default = {20 km for P > 0.15 s}, and {950 km for P < 0.15 s})')
parser.add_argument('--hmax', metavar="<hmax>", type=float, default=None, help='maximum emission height in km (default = 1000 km)')
parser.add_argument('--snr', metavar="<snr>", type=float, default=None, help='maximum emission height in km (default = None)')
parser.add_argument('-dm', metavar="<dm>", type=float, default=0.0, help='dispersion measure in cm^-3 pc (default = 0.0)')
#parser.add_argument('-o','--outfile', metavar="<name_suffix>", dest='output', action='store', type=argparse.FileType('w'), help="Write to file.")
parser.add_argument('--npatch', metavar="<npatch>", type=int, default='10', help='number of emission patches (default=10)' )
parser.add_argument('--do_ab', default=None, help='include aberration ofset (default = None)')
#parser.add_argument('-d', metavar="<dir>", default='/home/', help='Directory to save plots.')
#parser.add_argument('-x', "--x11", action='store_true', help='X11 window plot override switch.')

args = parser.parse_args()

P = args.p
ncomp = args.nc
npatch = args.npatch
iseed = args.iseed
hmin = args.hmin
hmax = args.hmax
freq = args.freq
alpha = args.alpha
beta = args.beta
snr = args.snr
dm = args.dm
do_ab = args.do_ab
#output = args.output
#plotwindow = args.x11
#dir = args.dir
#====================================================================================================================================================
#						Get the emission heights
#====================================================================================================================================================
heights = emission_height(P, freq, ncomp, iseed, hmin, hmax)
print heights
#====================================================================================================================================================
#						Find the opening angle
#====================================================================================================================================================
opa = rho(P, heights)

#====================================================================================================================================================
#					Determine the width of the emission patches
#====================================================================================================================================================
patchwidths = patch_width(P, heights)

#====================================================================================================================================================
#					   Determine the center of each patch
#====================================================================================================================================================
cx, cy = patch_center(P, heights, npatch)

#====================================================================================================================================================
# 						Find the line of sight	
#====================================================================================================================================================
xlos, ylos, thetalos = los(alpha, beta)

#====================================================================================================================================================
#					Plot the beam and line of sight profile
#====================================================================================================================================================
fig = plotpatch(P, alpha, beta, heights, cx, cy, snr, do_ab)
plt.show()

#with open('fig', 'w') as output_file:
#    output_file.write("%s\n" % fig)
#if output:
#    plt.savefig(fig)
##    if output:
#        plt.savefig(filename, format='pdf')
#else:
#    plt.savefig(filename, format='pdf')

#if opts.outfile:
#plt.savefig(filename, format='pdf')

