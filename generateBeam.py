#!/usr/bin/env python

# Program to generate a 2D beam and the corresponding 1D line 
# of sight profiles at varying frequency.

import beamModel as bm
import numpy as np
import argparse
import time, os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import scipy.signal as sci_sig
import scipy.stats as stats
from matplotlib.ticker import MultipleLocator

#==============================================================================================================================================
#                                          IMPORTANT FUNCTION:
#==============================================================================================================================================
def generateBeam(P, alpha, beta, freq, heights, npatch, snr, do_ab, iseed, fanBeam=None, hollowCone=None):
    """Function to plot the patches for a given rotation period.
    
       A rgs:
       -----
       P          : rotational period (seconds)
       alpha      : inclination angle (degrees)
       beta       : impact parameter (degrees)
       heights    : emission heights (in km)
       centerx    : the patch center projection on the x-axis 
       centery    : the patch center projection on the y-axis
       snr        : signal to noise ratio       
       iseed      : seed for the random number generator
       do_ab      : option to include abberration effects
       fanBeam    : option to use fan beam model 
       hollowCone : option to use hollow cone model

       Returns:
       --------
       prof    : an array containing 1D profile
       Z       : a 2D array of the beam 
    
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
    patchwidths = bm.patch_width(P, heights, hollowCone)
#   An arbitrary peak of the profile:
    peakAmp = 1.
#   Get the line of sight:
    xlos, ylos, thetalos = bm.los(alpha, beta, res)
#   Get the centre of the emission patches on the xy-plane
    centerx, centery = bm.patch_center(P, heights, npatch, iseed, fanBeam, hollowCone)
#   Get the ofset due to abberation:
    ab_xofset, ab_yofset = bm.aberration(heights, P, alpha)
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
parser.add_argument('-iseed', metavar="<iseed>", type=int, default=None, help='integer seed for a pseudo-random number generator (default = 4)')
parser.add_argument('-snr', metavar="<snr>", type=float, default=None, help='signal to noise ratio (default = None)')
parser.add_argument('-dmFile', metavar="<psrcat file>", default='psrcatdm.dat', type=str, help='A file containing PSRCAT dm values.')
parser.add_argument('-dm', metavar="<dm>", type=int, help='A dm to use for scattering.')
parser.add_argument('--outfile', action="store_true", help="Write delta dm to file.")
parser.add_argument('--do_ab', action="store_true", help='include aberration ofset (default = None)')
parser.add_argument('--doFan', action="store_true", help='Fan beam - default: patchy beam')
parser.add_argument('--scatter', action="store_true", help='include scattering (default = None)')
parser.add_argument('--writeprofile', action="store_true", help='Option to write out profile array into a file profile.txt.')
parser.add_argument('--doHC', action="store_true", help='Hollow Cone beam - default: patchy beam')
parser.add_argument('--getPlot', action="store_true", help='Option plot and save the beam / profiles')
parser.add_argument('--showrfm', action="store_true", help='Option to produce rfm .gif image')
parser.add_argument('--diagnostic', action="store_true", help='Option to show diagnostic plots')
args = parser.parse_args()
P = args.p
ncomp = args.nc
npatch = args.npatch
#iseed = args.iseed
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
hollowCone = args.doHC

if not args.iseed:
    # Using current time as seed to avoid repeating the same random number generation
    iseed = int(time.time())
else:
    iseed = args.iseed

#====================================================================================================================================================
#                                                        MAIN BODY: 
#====================================================================================================================================================
#=====================
# initialize params:
#=====================
beam = []
prof = []
w10 = []
res = 1e3 # resolution
t_res = P/res # time-reso;ution
phase = np.linspace(-180, 180, num=res) # rotation phase in degrees
max_freq = (nch - 1) * chbw + min_freq # maximum frequency
freq = np.linspace(min_freq, max_freq, nch) #channel frequencies in GHz!!!

#=======================================
#     1. Find the emission height:
#=======================================
H = bm.emission_height(P, ncomp, iseed, hmin, hmax, fanBeam, hollowCone)
#========================================
#     2. Get profile at each frequency:
#========================================
for i in np.arange(len(freq)):
    heights = bm.height_f(H, freq[i]) # frequency dependent H
    pr, Z = generateBeam(P, alpha, beta, freq[i], heights, npatch, snr, do_ab, iseed, fanBeam, hollowCone)
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
    if not args.dm:
        rand_dm = bm.getadm(args.dmFile, iseed, nbins=500, n=1) # random dm value from a dist. of known psr dm
    else:
        rand_dm = args.dm
    #rand_dm = bm.getadm(args.dmFile, iseed, nbins=500, n=1) # random dm value from a dist. of known psr dm
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
    rms = bm.noise_rms(snr)                            # Determine the noise rms
    profile = bm.add_noise(sc_prof, rms, res)          # add noise to each profile
    for p in profile:
       SN.append(bm.signal_to_noise(np.max(p), rms))   # snr for each of the profiles(to use later below)

# Write out the profile into a file (scattered, with added noise, if specified)
if args.writeprofile:
    np.savetxt('profile_file.txt', np.array(profile))

#==================================================================
#      5. Fit a DM Curve:
#==================================================================
# Increase the resolution
highres_phase = np.linspace(-180,180,1000*res)
resampled = np.zeros((int(nch),int(1000*res)))
for nfr in range(len(freq)):
    resampled[nfr] = sci_sig.resample(profile[nfr], int(1000*res))


# delta dm search only for profiles with snr above threshold
#----------------- FIRST ITERATION --------------------------------
# Find a region that contain the best DM
if all(i > 10 for i in SN):
    #average_profile = []
    peaks_of_average = []
    phase_bin0 = bm.find_phase_bin(resampled[nch - 1])
    phase_bin1 = bm.find_phase_bin(resampled[0])
    dm_range = bm.find_delta_dm(P, resampled, highres_phase, phase_bin0, phase_bin1, freq[nch - 1], freq[0], nch)
#    print "First iteration: dm_range[0] = %.5f, dm range[-1] = %.5f " %(dm_range[0], dm_range[-1])
#    print "dm step:", (dm_range[1]-dm_range[0])
    for dm_id in range(len(dm_range)):
        shifted_profiles = []
        if args.diagnostic:
            fig = plt.figure(figsize=(10,5))
            #st = fig_dm.suptitle("Dm trial, delta DM = %.5f" %(dm_range[dm_id]), fontsize="x-large")
            plt.title('DM trial, delta DM = %.5f' %dm_range[dm_id])
            plt.xlabel('phase (degrees)')
            plt.ylabel('Intensity')
            plt.xlim(-75,75)
            plt.grid()
        for freq_id in range(nch):
            bin_shift = bm.delay(freq[nch - 1], freq[freq_id], dm_range[dm_id], t_res/1000.) # Res increased by 1000 more bins
            shifted_profiles.append(np.roll(resampled[freq_id], bin_shift))
            #plt.subplot(1,2,1)
            plt.plot(highres_phase, shifted_profiles[freq_id] + 2*freq_id)
        #average_profile.append(bm.avg_prof(shifted_profiles))
        #peaks_of_average.append(bm.find_peak(average_profile[dm_id]))
        average = bm.avg_prof(shifted_profiles)
        peaks_of_average.append(bm.find_peak(average))
#        if args.diagnostic:
#            plt.subplot(1,2,2)
#            plt.plot(highres_phase, average_profile[dm_id])
#            plt.ylim(np.min(profile), np.max(profile))
#            plt.xlim(-100, 100)
#            plt.grid()
        if args.diagnostic:
            if args.doHC:
                fig.savefig('Dm_trial_HC_%d_seed_%d_DM_%.5f_1.png' %(dm_id, int(iseed),dm_range[dm_id]))
                fig.clear()
                plt.close(fig)
            else:
                fig.savefig('Dm_trial_KJ07_%d_seed_%d_DM_%.5f_1.png' %(dm_id, int(iseed),dm_range[dm_id]))
                fig.clear()
                plt.close(fig)
    if args.getPlot:
        # Create a snr vs dm plot for visualization
        snrfig1 = plt.figure()
        plt.plot(dm_range, peaks_of_average, '.')
        plt.title('Best dm trial')
        plt.xlabel('delta dm (pc cm^-3)')
        plt.ylabel('SNR')
        if args.doHC:
            snrfig1.savefig('SNR_DM_HC_seed_%f_1.png' %(iseed))
            snrfig1.clear()
            plt.close(snrfig1)
        else:
            snrfig1.savefig('SNR_DM_KJ07_seed_%f_1.png' %(iseed))
            snrfig1.clear()
            plt.close(snrfig1)

#----------------- SECOND ITERATION ----------------------------
# Search around this region for a best dm
dm_bin = bm.find_phase_bin(peaks_of_average) # find the bin where the SNR-DM curve peak (best DM from the 1st iteration)
#binshift = bm.delay(freq[-1], freq[0], 1, t_res/1000)
#dm_step = 1/float(binshift)
# print dm_range[dm_bin-1], dm_range[dm_bin],  dm_range[dm_bin+1], dm_bin
# dm_search = np.arange(dm_range[dm_bin] - 100*dm_step,dm_range[dm_bin] + 100*dm_step ,dm_step)
dm_search = np.linspace(dm_range[dm_bin - 1], dm_range[dm_bin + 1], 200) # search for the range to try
dm_range = dm_search # set a new range
#dm_range = np.linspace(-.025, .025, 200)
#print "Second iterationd: dm_range[0] = %.5f, dm range[-1] = %.5f " %(dm_range[0], dm_range[-1])
#print "Step:", (dm_range[1]-dm_range[0])

if all(i > 10 for i in SN):
    #average_profile = []
    peaks_of_average = []
#    phase_bin0 = bm.find_phase_bin(resampled[nch - 1])
#    phase_bin1 = bm.find_phase_bin(resampled[0])
    # dm_range = bm.find_delta_dm(P, resampled, highres_phase, phase_bin0, phase_bin1, freq[nch - 1], freq[0], nch)
    for dm_id in range(len(dm_range)):
        shifted_profiles = []
#        if args.diagnostic:
#            fig_dm2 = plt.figure(figsize=(10,5))
#            #st = fig_dm.suptitle("Dm trial, delta DM = %.5f" %(dm_range[dm_id]), fontsize="x-large")
#            plt.title('DM trial, delta DM = %.5f' %dm_range[dm_id])
#            plt.xlabel('phase (degrees)')
#            plt.ylabel('Intensity')
#            plt.xlim(-180,180)
#            plt.grid()
        for freq_id in range(nch):
            bin_shift = bm.delay(freq[nch - 1], freq[freq_id], dm_range[dm_id], t_res/1000.) # Res increased by 1000 more bins
            shifted_profiles.append(np.roll(resampled[freq_id], bin_shift))
#            plt.subplot(1,2,1)
#            plt.plot(highres_phase, shifted_profiles[freq_id] + 2*freq_id)
        #average_profile.append(bm.avg_prof(shifted_profiles))
        #peaks_of_average.append(bm.find_peak(average_profile[dm_id]))
        average = bm.avg_prof(shifted_profiles)
        peaks_of_average.append(bm.find_peak(average))
#        if args.diagnostic:
#            plt.subplot(1,2,2)
#            plt.plot(highres_phase, average_profile[dm_id])
#            plt.ylim(np.min(profile), np.max(profile))
#            plt.xlim(-100, 100)
#            plt.grid()
#        if args.doHC:
#            fig_dm2.savefig('Dm_trial_HC_seed_%f_DM_%.5f_1.png' %(iseed,dm_range[dm_id]))
#        else:
#            fig_dm2.savefig('Dm_trial_KJ07_seed_%f_DM_%.5f_1.png' %(iseed,dm_range[dm_id]))
    if args.getPlot:
        # Create a snr vs dm plot for visualization
        snrfig2 = plt.figure()
        plt.plot(dm_range, peaks_of_average, '.')
        plt.title('Best dm trial')
        plt.xlabel('delta dm (pc cm^-3)')
        plt.ylabel('SNR')
        if args.doHC:
            snrfig2.savefig('SNR_DM_HC_seed_%f_2.png' %(iseed))
            snrfig2.clear()
            plt.close(snrfig2)
        else:
            snrfig2.savefig('SNR_DM_KJ07_seed_%f_2.png' %(iseed))
            snrfig2.clear()
            plt.close(snrfig2)


    # Find the best dm (dm that maximises SNR)
    for i in range(len(peaks_of_average)):
        if peaks_of_average[i] == np.max(peaks_of_average):
            best_dm = dm_range[i]
#            print "Best dm = ", best_dm
    
    # Write out important parameters into a file    
    if args.outfile:
        if args.doHC:
            model = 'Hollow_cone'
        else:
            model = 'Patchy_beam'
        pulsarParams = np.asarray([model, min_freq, chbw, nch, P, alpha, beta, \
                                  iseed, rand_dm, best_dm])
        f = open('dm_dat.txt', 'a')
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
        ax2.plot(phase, profile[k] + k, label='frequency = %0.2f GHz' %freq[k], color='black')
	plt.title('Profile evolution with frequency' )
	plt.xlabel('Phase (degrees)')
        plt.ylabel('Frequency (30 - 80 MHz)')
        plt.xlim(-50,50)
        plt.tick_params(axis='y', which='both', left='off', top='off', labelleft='off')
        plt.grid('on', which='both')
        #============================================
        #    2D emission region:
        #============================================
        xlos, ylos, thetalos = bm.los(alpha, beta, res)
        fig3 = plt.figure(figsize=(10,5))
        ax31 = fig3.add_subplot(1,2,1)
        plt.plot(xlos, ylos, '+r')
        plt.imshow(beam[k].T, extent=[-180, 180, 180, -180])
        plt.xlabel('X (degrees)')
        plt.title('Beam')
        plt.ylabel('Y (degrees)')
        # find zoomed extent for plot
        nonZero = np.where(beam[0]!= 0.0)    # Set the maximum zoom to beam at lowest frequency
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
        ax32 = fig3.add_subplot(1,2,2)
        plt.plot(phase, profile[0], color='blue')
        plt.title('Profile at %.3f GHz' % freq[0])
        plt.ylim(np.min(profile), np.max(profile))
        plt.xlim(-180, 180)
        plt.xlabel('Phase ')
        plt.ylabel('Intensity')
        if args.doHC:
            suffix = 'HC'
        else:
            suffix = 'KJ07'
        fig2.savefig('sequence_%s_%d.png' %(suffix, int(iseed)))
        fig3.savefig('beam_%.3f_GHz_%s_%d_.png' %(freq[0], suffix, int(iseed)))
        #fig2.clear()
        #fig3.clear()
        plt.close(fig2)
        plt.close(fig3)

    if args.showrfm:
        # show rfm on 2d beam plots
        fig4 = plt.figure(figsize=(10,5))
        ax41 = fig4.add_subplot(1,2,1)
        for bid in range(nch):
            xlos, ylos, thetalos = bm.los(alpha, beta, res)
            fig4 = plt.figure(figsize=(10,5))
            ax41 = fig4.add_subplot(1,2,1)
            plt.plot(xlos, ylos, '+r')
            plt.imshow(beam[bid].T, extent=[-180, 180, 180, -180])
            plt.xlabel('X (degrees)')
            plt.title('Radio pulsar beam')
            plt.ylabel('Y (degrees)')
            # find zoomed extent for plot
            nonZero = np.where(beam[0]!= 0.0)     # set the maximum zoom to beam at lowest frequency
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
            # 1D profile plot
            ax42 = fig4.add_subplot(1,2,2)
            plt.plot(phase, profile[bid])
            plt.title('Profile at %.3f GHz' % freq[bid])
            plt.xlim(-180, 180)
            plt.ylim(np.min(profile), np.max(profile))
            plt.xlabel('Phase ')
            plt.ylabel('Intensity')
            if args.doHC:
                suffix = 'HC'
            else:
                suffix = 'KJ07'
            fig4.savefig('beam_%.3f_GHz_%s_%d_.png' %(freq[bid], suffix, int(iseed))) 
            fig4.clear()
            plt.close(fig4)
# Clear up 
cwd = os.getcwd()
folder = '/' + str(iseed)
new_directory = cwd + folder
os.mkdir(new_directory)
files_to_move = '*%s*.png' %(str(iseed))
os.system('mv %s %s' %(files_to_move, new_directory))
