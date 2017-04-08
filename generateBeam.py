#!/usr/bin/env python

# Program to generate a 2D beam and the corresponding 1D line 
# of sight profiles at varying frequency.

import beamModel as bm
import numpy as np
import argparse
import time
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import scipy.signal as sci_sig
import scipy.stats as stats
from matplotlib.ticker import MultipleLocator

#==============================================================================================================================================
#                                          IMPORTANT FUNCTION:
#==============================================================================================================================================
def generateBeam(P, alpha, beta, freq, heights, npatch, snr, do_ab, iseed, fanBeam):
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
       iseed   : seed for the random number generator
       do_ab   : option to include abberration effects
       fanBeam : option to use fan beam model 

       Returns:
       --------
       prof    : an array containing 1D profile
       Z       : a 2D beam 
    
    """    
    
#   initialize parameters:
    xmin = -180.
    xmax = 180.
#    res = 1e4 #resolution
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
    patchwidths = bm.patch_width(P, heights)
#   An arbitrary peak of the profile:
    peakAmp = 1.
#   Get the line of sight:
    xlos, ylos, thetalos = bm.los(alpha, beta, res)
#   Get the centre of the emission patches on the xy-plane
    centerx, centery = bm.patch_center(P, heights, npatch, iseed, fanBeam)
#   Get the ofset due to abberation:
    ab_xofset, ab_yofset = bm.aberration(heights, P, alpha)
    if fanBeam:
        comp = 1
    else:
        comp = comp

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
            distance = (np.sqrt((X - pc[0])**2 + (Y - pc[1])**2))/sigmax
            distance[np.where(distance > 3.0)] = 0.0
            distance[np.where(distance != 0.0)] = peakAmp
            if not do_ab:
                Z += distance * np.exp(-((X - pc[0])**2 / (2 * sigmax**2) + (Y - pc[1])**2 / (2 * sigmay**2)))
            else:
                Z += distance * np.exp(-((X - pc[0] - ab_xofset[cid])**2 / (2 * sigmax**2) + (Y - pc[1] - ab_yofset[cid])**2 / (2 * sigmay**2)))
#   1D profile from 2D patch, closest to the line of sight (select nearest neighbors):
    ZxIdx = np.array((xlos-xmin)/dx, dtype=int) # x index
    ZyIdx = np.array((ylos-ymin)/dy, dtype=int) # y index
    prof = Z[ZxIdx, ZyIdx]
    
    return prof, Z

#===============================================================================================================================================
#                                       MAIN:
#===============================================================================================================================================

parser = argparse.ArgumentParser(description='Plot the patchy emission region as well as the line of sight profile.\
                                 Running the file without specified argument will produce an output beam and profile from default parameters.')
parser.add_argument('-alpha', metavar="<alpha>", type=float, default='45', help='inclination angle in degrees (default = 45)')
parser.add_argument('-beta', metavar="<beta>", type=float, default='5', help='impact parameter in degrees (default = 5)')
parser.add_argument('-p', metavar="<p>", type=float, default='0.16', help='period in seconds (default = 0.16 s)')
parser.add_argument('-hmin', metavar="<hmin>", type=float, default=None, help='minimum emission height in km\
                    (default = {20 km for P > 0.15 s}, and {950 km for P < 0.15 s})')
parser.add_argument('-hmax', metavar="<hmax>", type=float, default=None, help='maximum emission height in km (default = 1000 km)')
parser.add_argument('-nc', metavar="<ncomp>", type=int, default='4', help='integer number of components (default = 4)')
parser.add_argument('-npatch', metavar="<npatch>", type=int, default='10', help='number of emission patches (default= 10)' )
parser.add_argument('-min_freq', metavar="<minfreq>", type=float, default='0.2', help='min frequency in GHz (default = 0.2 GHz)')
parser.add_argument('-chbw', metavar="<chanbw>", type=float, default='0.8', help='channel bandwidth in GHz (default = 0.8 GHz)')
parser.add_argument('-nch', metavar="<nch>", type=int, default='5', help='number of channels (default = 5)')
parser.add_argument('-iseed', metavar="<iseed>", type=int, default='4', help='integer seed for a pseudo-random number generator (default = 4)')
parser.add_argument('-snr', metavar="<snr>", type=float, default=None, help='signal to noise ratio (default = None)')
parser.add_argument('-dmFile', metavar="<psrcat file>", default='psrcatdm.dat', type=str, help='A file containing PSRCAT dm values.')
parser.add_argument('-outfile', metavar="<output file>", help="Write to file.")
parser.add_argument('--do_ab', action="store_true", help='include aberration ofset (default = None)')
parser.add_argument('--scatter', action="store_true", help='include scattering (default = None)')
parser.add_argument('--doFan', action="store_true", help='Fan beam - default: patchy beam')
parser.add_argument('--getPlot', action="store_true", help='Option plot and save the beam / profiles')
args = parser.parse_args()
P = args.p
ncomp = args.nc
npatch = args.npatch
#iseed = args.iseed
iseed = time.time() # Using current time as seed to avoid repeating the same random number generation
hmin = args.hmin
hmax = args.hmax
alpha = args.alpha
beta = args.beta
snr = args.snr
do_ab = args.do_ab
nch = args.nch
min_freq = args.min_freq
chbw = args.chbw
scr = args.scatter
fanBeam = args.doFan
fileName = args.outfile

#====================================================================================================================================================
#                                                        MAIN BODY: 
#====================================================================================================================================================
#=====================
# initialize params:
#=====================
beam = []
prof = []
w10 = []
res = 1e3
t_res = P/res
phase = np.linspace(-180, 180, num=res)
max_freq = (nch - 1) * chbw + min_freq
freq = np.linspace(min_freq, max_freq, nch) #channel frequency in GHz!!!

#=======================================
#     1. Find the emission height:
#=======================================
H = bm.emission_height(P, ncomp, iseed, hmin, hmax)
print "Height: ", H
#========================================
#     2. Get profile at each frequency:
#========================================
for i in np.arange(len(freq)):
    heights = bm.height_f(H, freq[i]) # frequency dependent H
    pr, Z = generateBeam(P, alpha, beta, freq[i], heights, npatch, snr, do_ab, iseed, fanBeam)
    w10.append(bm.find_width(pr)) # Width at 10% of the peak 
    prof.append(pr)               # Profile for that frequency
    beam.append(Z)                # 2D beam 

#==========================================
#     3. Scatter the line of sight profile: 
#==========================================
train = []
bf = []
if not scr:
    sc_prof = prof # returns the profile without scattering 

    rand_dm = 0.0   # random Dm for scattering (time_scale : No scattering)
else:
    sc_prof = []
    rand_dm = bm.getadm(args.dmFile, iseed, nbins=20, n=1) # random dm value from a dist. of known psr dm
    # Follow the scattering routine:
    for pid in np.arange(len(prof)):
        # Compute a train of pulses
        train.append(bm.pulsetrain(3, res, prof[pid]))

    for fid in np.arange(len(freq)):
        #Find the scattering timescale for the random dm
        tau = bm.sc_time(freq[fid], rand_dm, iseed)
        #Determine a broadening function
        bf = bm.broadening(tau, P, res)
        #scatter the train of pulses with this function
        sc_train = bm.scatter(train[fid], bf)
        #Extract a pulse profile
        sc_prof.append(bm.extractpulse(sc_train, 2, res))


#===========================================
#     4. Add noise:
#===========================================
# Find amplitudes of the profiles.
#peaks = []
#for j in np.arange(len(prof)):
#    peaks.append(bm.find_peak(sc_prof[j]))

SN = []
if snr == None:
    profile = sc_prof
else:
    #rms = bm.noise_rms(snr, np.max(peaks))
    rms = bm.noise_rms(snr)                            # Determine the noise rms
    profile = bm.add_noise(sc_prof, rms, res)          # add noise to each profile
    for p in profile:
       SN.append(bm.signal_to_noise(np.max(p), rms))   # snr for each of the profiles
# YOU MAY NEED TO CHECK THE here before doing dm search instead of putting that in a loop/under condition
#==================================================================
#      5. Fit a DM Curve:
#==================================================================
#if any(SN > 10):
highres_phase = np.linspace(-180,180,10*res)
resampled = np.zeros((int(nch),int(10*res)))
for nfr in range(len(freq)):
    resampled[nfr] = sci_sig.resample(profile[nfr], int(10*res))

if all(i > 10 for i in SN):
    average_profile = []
    peaks_of_average = []
#    phase_bin0 = bm.find_phase_bin(profile[nch - 1])
#    phase_bin1 = bm.find_phase_bin(profile[0])
    phase_bin0 = bm.find_phase_bin(resampled[nch - 1])
    phase_bin1 = bm.find_phase_bin(resampled[0])
#    dm_range = bm.find_delta_dm(P, profile, phase, phase_bin0, phase_bin1, freq[nch - 1], freq[0], nch)
    dm_range = bm.find_delta_dm(P, resampled, highres_phase, phase_bin0, phase_bin1, freq[nch - 1], freq[0], nch)
    for dm_id in range(len(dm_range)):
        shifted_profile = []
        """fig1 = plt.figure(figsize=(10,5))
        st = fig1.suptitle("Dm trial, delta DM = %.5f" %(dm_range[dm_id]), fontsize="x-large")
        plt.title('DM trial, delta DM = %.5f' %dm_range[dm_id])
        plt.xlabel('phase (degrees)')
        plt.ylabel('Intensity')"""
        for freq_id in range(nch-1):
            bin_shift = bm.delay(freq[nch - 1], freq[freq_id], dm_range[dm_id], t_res/10.)
#            shifted_profile.append(np.roll(profile[freq_id], bin_shift))
            shifted_profile.append(np.roll(resampled[freq_id], bin_shift))
            """plt.subplot(1,2,1)
            plt.plot(highres_phase, shifted_profile[freq_id])
        plt.xlim(-180,180)
        #plt.ylim(0, 5)
        plt.grid()"""
        average_profile.append(bm.avg_prof(shifted_profile))
        peaks_of_average.append(bm.find_peak(average_profile[dm_id]))
        """plt.subplot(1,2,2)
        plt.plot(highres_phase, average_profile[dm_id])
        plt.xlim(-180, 180)
        #plt.ylim(0, 5)
        plt.grid()
        fig1.savefig('Dm_trial_DM_%d_%.5f.png' %(dm_id, dm_range[dm_id]))"""
    # Create a snr vs dm plot for visualization
    snrfig = plt.figure()
    plt.plot(dm_range, peaks_of_average)
    plt.title('Best dm trial')
    plt.xlabel('delta dm (pc cm^-3)')
    plt.ylabel('SNR')
    snrfig.savefig('SNR_DM.png')

    # Find the best dm (dm that maximises SNR)
    for i in range(len(peaks_of_average)):
        if peaks_of_average[i] == np.max(peaks_of_average):
            best_dm = dm_range[i]
    
    # Write out important parameters into a file    
    pulsarParams = np.asarray([P, alpha, beta, w10[0], w10[-1], iseed, rand_dm, best_dm])
    f = open(fileName, 'a')
    f.write(' '.join([str(item) for item in pulsarParams]) + ' \n')

#==================================================================
#                     PRODUCE PLOTS:
#==================================================================
if args.getPlot:
#===========================================
# SET A ZERO BASELINE AND PLOT THE PROFILE:
#===========================================
    fig2, ax2 = plt.subplots()
    for k in np.arange(len(profile)):
        colormap = plt.cm.gist_ncar
        #plt.gca().set_color_cycle([colormap(k+1)])
        ax2.plot(phase, profile[k] + k, label='frequency = %0.2f GHz' %freq[k])
        #ax2.plot(phase, profile[k], label='frequency = %0.2f GHz' %freq[k])
	plt.title('A sequence of %i pulse profile' %nch)
	plt.xlabel('Phase (degrees)')
        plt.ylabel('Profile number')
        plt.xlim(-180,180)
        plt.grid('on')
    fig2.savefig('sequence.png')
    #============================================
    #    2D emission region:
    #============================================
    xlos, ylos, thetalos = bm.los(alpha, beta, res)
    fig3 = plt.figure(figsize=(10,5))
    ax3 = fig3.add_subplot(1,2,1)
    plt.plot(xlos, ylos, '+r')
    plt.imshow(beam[0].T, extent=[-180, 180, 180, -180])
    plt.xlabel('X (degrees)')
    plt.title('Beam')
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
    # 1D plot
    ax4 = fig3.add_subplot(1,2,2)
    plt.plot(phase, profile[0])
    plt.title('Profile at %.3f GHz' % freq[0])
    plt.xlim(-180, 180)
    plt.xlabel('Phase ')
    plt.ylabel('Intensity')
    fig3.savefig('beam_%.3f_GHz_.png' %freq[0]) 
    """for bid in range(nch):
        xlos, ylos, thetalos = bm.los(alpha, beta, res)
        fig3 = plt.figure(figsize=(10,5))
        ax3 = fig3.add_subplot(1,2,1)
        plt.plot(xlos, ylos, '+r')
        plt.imshow(beam[bid].T, extent=[-180, 180, 180, -180])
        plt.xlabel('X (degrees)')
        plt.title('Radio pulsar beam')
        plt.ylabel('Y (degrees)')
        # find zoomed extent for plot
        plt.xlim(-40, 40)
        plt.ylim(-40, 40)
        plt.colorbar()
        # 1D profile plot
        ax4 = fig3.add_subplot(1,2,2)
        plt.plot(phase, profile[bid])
        plt.title('Profile at %.3f GHz' % freq[bid])
        plt.xlim(-180, 180)
        plt.ylim(0, 5)
        plt.xlabel('Phase ')
        plt.ylabel('Intensity')
        fig3.savefig('beam_%.3f_GHz_.png' %freq[bid]) 
        print 'Done!'"""
