#!/usr/bin/env python
# Program to generate a 2D beam and the corresponding 1D line
# of sight profiles at varying frequency.
import beamModel as bm
import numpy as np
import argparse
import time
import sys
from scipy.signal import resample  
from scipy.optimize import curve_fit
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
# ==============================================================================================================================================
#                                          IMPORTANT FUNCTION:
# ==============================================================================================================================================
def generateBeam(P, alpha, beta, freq, heights, npatch, snr, do_ab, iseed, fanBeam=None, hollowCone=None):
    """Function to plot the patches for a given rotation period.
    
       A rgs:
       -----
       P          : rotational period (seconds)
       alpha      : inclination angle (degrees)
       beta       : impact parameter (degrees)
       heights    : emission heights (in km)
       snr        : signal to noise ratio       
       iseed      : seed for the random number generator
       do_ab      : option to include abberration effects
       fanBeam    : option to use fan beam model 
       hollowCone : option to use hollow cone model
       iseed      : seed for random number generator

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
#   find the highest height
    maxheight = np.max(heights)
#   find opening angle for maximum height
    opa_max = bm.rho(P, maxheight)
#   find the observed pulse/beam width using Gil formula (eq 3.26) [assumes a fully illuminated beam!]
    sin2W4 = (np.sin(np.deg2rad(opa_max))**2-np.sin(np.deg2rad(beta/2.0))**2)/ (np.sin(np.deg2rad(alpha))*np.sin(np.deg2rad(alpha+beta)))
    W = 4. * np.rad2deg(np.arcsin(np.sqrt(abs(sin2W4))))
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
            # if distance of current box from the patch center is
            # larger than 3 times the patchwidth, I do not want any
            # power contributed
            # first compute a grid of distance from the center of the patch, in units of the width
            distance = (np.sqrt((X - pc[0])**2 + (Y - pc[1])**2))/sigmax
            distance[np.where(distance > 3.0)] = 0.0
            distance[np.where(distance != 0.0)] = peakAmp
            if not do_ab:
                Z += distance * np.exp(-((X - pc[0])**2 / (2 * sigmax**2) + (Y - pc[1])**2 / (2 * sigmay**2)))
            else:
                Z += distance * np.exp(-((X - pc[0] - ab_xofset[cid])**2 / (2 * sigmax**2) + (Y - pc[1] - ab_yofset[cid])**2 / (2 * sigmay**2)))
#   1D profile from 2D patch, closest to the line of sight (select nearest neighbors):
#    ZxIdx = np.array((xlos-xmin)/dx, dtype=int) # x index
#    ZyIdx = np.array((ylos-ymin)/dy, dtype=int) # y index
#    prof = Z[ZxIdx, ZyIdx]
    prof = np.zeros(int(res))
# Find the appropriate range to fill the profile
# First converst our expected width W into bins
    halfbinsW = int(W/360. * res/2.)
    profilerange = range(int(res/2)-halfbinsW, int(res/2)+halfbinsW, 1)
    for i in profilerange:
        for cid, comp in enumerate(heights):
            #       widths for circular patches:        
            sigmax = patchwidths[cid]
            sigmay = patchwidths[cid]

#       center of the patch:
            patchCenterX = centerx[cid]
            patchCenterY = centery[cid]
            for pc in zip(patchCenterX, patchCenterY):
                distance = (np.sqrt((xlos[i] - pc[0])**2 + (ylos[i] - pc[1])**2))/sigmax
                if distance > 3.0:
                    distance  = 0.0
                else:
                    distance = peakAmp
                    if not do_ab:
                        prof[i] += distance * np.exp(-((xlos[i] - pc[0])**2 / (2 * sigmax**2) + (ylos[i] - pc[1])**2 / (2 * sigmay**2)))
                    else:
                        prof[i] += distance * np.exp(-((xlos[i] - pc[0] - ab_xofset[cid])**2 / (2 * sigmax**2) + (ylos[i] - pc[1] - ab_yofset[cid])**2 / (2 * sigmay**2)))
    
    return prof, Z, W

# ===============================================================================================================================================
#                                       MAIN:
# ===============================================================================================================================================

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
parser.add_argument('-dm', metavar="<dm>", type=float, help='A dm to use for scattering.')
parser.add_argument('-freqrange', metavar="<range>", nargs='+', help='Option to select range of frequency (used only for labeling)')
parser.add_argument('--outfile', action="store_true", help="Write delta dm to file.")
parser.add_argument('--do_ab', action="store_true", help='include aberration ofset (default = None)')
parser.add_argument('--doFan', action="store_true", help='Fan beam - default: patchy beam')
parser.add_argument('--scatter', action="store_true", help='include scattering (default = None)')
parser.add_argument('--disperse', action="store_true", help='include dispersion effects (default = None)')
parser.add_argument('--writeprofile', action="store_true", help='Option to write out profile array into a file profile.txt.')
parser.add_argument('--doHC', action="store_true", help='Hollow Cone beam - default: patchy beam')
parser.add_argument('--getPlot', action="store_true", help='Option plot and save the beam / profiles')
parser.add_argument('--showrfm', action="store_true", help='Option to produce rfm .gif image')
parser.add_argument('--diagnostic', action="store_true", help='Option to show diagnostic plots')
parser.add_argument('--random_beta', action="store_true", help='Option to chose a random beta')
parser.add_argument('--random_alpha', action="store_true", help='Option to chose a random alpha')
parser.add_argument('--template_matching', action="store_true", help='Option to use template matching method for delta dm search. S/N maximization method default.')
args = parser.parse_args()
P = args.p
ncomp = args.nc
npatch = args.npatch
# iseed = args.iseed
hmin = args.hmin
hmax = args.hmax
#alpha = args.alpha
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

# ====================================================================================================================================================
#                                                        MAIN BODY: 
# ====================================================================================================================================================
# =====================
# initialize params:
# =====================
beam = []
prof = []
res = 1e3 # resolution
t_res = P/res # time-reso;ution
phase = np.linspace(-180, 180, num=res) # rotation phase in degrees
max_freq = (nch - 1) * chbw + min_freq # maximum frequency
freq = np.linspace(min_freq, max_freq, nch) # channel frequencies in GHz!!!

#=======================================
#     1. Find the emission height:
#=======================================
H = bm.emission_height(P, ncomp, iseed, hmin, hmax, fanBeam, hollowCone)
emission_heights = []
for i in np.arange(len(freq)):
#   frequency dependent H 
    emission_heights.append(bm.height_f(H, freq[i]))

# ========================================
#      2. Get profile at each frequency:
# ========================================
# ------- Get the impact parameter -------
if args.random_beta:
    # Largest opening angle at lowest frequency:
    h = np.max(emission_heights)
    opening_angle = bm.rho(P, h)
    np.random.seed(iseed)
    beta = np.random.uniform(-np.max(opening_angle), np.max(opening_angle))
else:
    beta = args.beta
if args.random_alpha:
    alpha = np.rad2deg(np.arccos(np.random.uniform()))
else:
    alpha = args.alpha

# -------- Get the beam & profile --------
for i in range(nch):
    pr, Z, W = generateBeam(P, alpha, beta, freq[i], emission_heights[i], npatch, snr, do_ab, iseed, fanBeam, hollowCone)
#    w10.append(bm.find_width(pr)) # Width at 10% of the peak 
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
        rand_dm = bm.getadm(args.dmFile, iseed, nbins=1000, n=1) # random dm value from a dist. of known psr dm
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

#==========================================
#    Disperse the signal
#==========================================
if args.disperse:
    dispersed_profile = []
    if args.dm:
        for fid in range(len(freq)):
            # Delay profiles using high freq profile as reference
            bin_delay = bm.delay(freq[-1], freq[fid], args.dm, t_res)
            print bin_delay
            dispersed_profile.append(np.roll(sc_prof[fid], bin_delay))
    else:
        print 'Give a DM to disperse the profiles!'
        sys.exit(1)
else:
    dispersed_profile = sc_prof


#====================================================================
#           Determine the profile widths of the noiseless profile
#====================================================================
peaks = np.max(sc_prof, axis=1)
w10_noiseless = np.zeros_like(peaks)
for i, peak in enumerate(peaks):
    w10_phase_initial = phase[np.where(sc_prof[i] > peak/10)][0]
    w10_phase_final = phase[np.where(sc_prof[i] > peak/10)][-1]
    w10_phase = w10_phase_final - w10_phase_initial
    w10_noiseless[i] = w10_phase/360. * P # w10 in seconds
#===========================================
#     4. Add noise:
#===========================================i


# Find amplitudes of the profiles.
#peaks = []
#for j in np.arange(len(prof)):
#    peaks.append(bm.find_peak(sc_prof[j]))

SN = []
if snr == None:
    profile = dispersed_profile
else:
    rms = bm.noise_rms(snr)                            # Determine the noise rms
    profile = bm.add_noise(dispersed_profile, rms, res)          # add noise to each profile
    for p in profile:
       SN.append(bm.signal_to_noise(np.max(p), rms))   # snr for each of the profiles(to use later below)

# Write out the profile into a file (scattered, with added noise, if specified)
if args.writeprofile:
    np.savetxt('profile_file_%i.txt' %iseed, np.array(profile))
    np.savetxt('channels_%i.txt' %iseed, freq)
    sys.exit()
#==================================================================
#      5. Fit a DM Curve:
#==================================================================
# Increase the resolution
if not args.template_matching:
    highres_phase = np.linspace(-180,180,1000*res)
    resampled = np.zeros((int(nch),int(1000*res)))
    for nfr in range(len(freq)):
        resampled[nfr] = resample(profile[nfr], int(1000*res))

    # delta dm search only for profiles with snr above threshold
    #----------------- FIRST ITERATION --------------------------------
    # Find a region that contain the best DM
    #if all(i > 3 for i in SN):
    if np.min(np.amax(profile, axis=1)) > 0.5:
        #average_profile = []
        peaks_of_average = []
        phase_bin0 = bm.find_phase_bin(resampled[nch - 1])
        phase_bin1 = bm.find_phase_bin(resampled[0])
    #   Find the time delay between max and min frequency profiles
        delta_t = bm.find_delta_t(W, P)
    #   Scattering broadens the pulse profile, introducing a possible shift in the peak
    #   relative to the unscattered profile.
    #   If scattering included = The range of DM to tries shift by tau from DM used for scatteering
        if args.scatter:
            # Get the delay caused by scattering:
            delta_tau = bm.sc_time(freq[0], rand_dm, iseed)
            # Total time delay:
            delta_t = delta_t + delta_tau
            # Relate to dispersion measure:
            D = 4.148808 * 1e3 # +/- 3e-6 MHz^2 pc^-1 cm^3 s
            dm_scatter = delta_t/ (D * ((freq[-1] * 1e3)**(-2) - (freq[0] * 1e3)**(-2)))
            dm_range = np.linspace(-0.5*dm_scatter, dm_scatter*0.5 , num=20)
        else:
            D = 4.148808 * 1e3 # +/- 3e-6 MHz^2 pc^-1 cm^3 s
            dm_unscatter = delta_t / (D * ((freq[-1]*1e3)**(-2) - (freq[0]*1e3)**(-2))) 
            dm_range = np.linspace(-0.5*dm_unscatter, dm_unscatter*0.5 , num=20)
        for dm_id in range(len(dm_range)):
            shifted_profiles = []
            if args.diagnostic:
                fig = plt.figure(figsize=(10,5))
                plt.title('DM trial, delta DM = %.5f' %dm_range[dm_id], fontsize=18)
                plt.xlabel('Phase (degrees)', fontsize=18)
                plt.ylabel('Frequencies (%.1f - %.1f MHz) ' %(freq[0]*1e3, freq[-1]*1e3), fontsize=18)
                plt.tick_params(axis='y', which='both', left='off', top='off', labelleft='off')
                plt.xlim(-180, 180)
                plt.xticks(fontsize = 10)
                plt.yticks(fontsize = 10)
                plt.grid()

            for freq_id in range(nch):
                bin_shift = bm.delay(freq[nch - 1], freq[freq_id], dm_range[dm_id], t_res/1000.) # Res increased by 1000 more bins
                shifted_profiles.append(np.roll(resampled[freq_id], bin_shift))
                #plt.subplot(1,2,1)
                plt.plot(highres_phase, shifted_profiles[freq_id] + 2*freq_id, color='grey')
            #average_profile.append(bm.avg_prof(shifted_profiles))
            #peaks_of_average.append(bm.find_peak(average_profile[dm_id]))
            average = bm.avg_prof(shifted_profiles)
            peaks_of_average.append(bm.find_peak(average))
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


#------------   Find region that we will search for best dm -----------
#       find the bin where the SNR-DM curve peak (best DM from the 1st iteration)
        dm_bin = bm.find_phase_bin(peaks_of_average) 
#       binshift = bm.delay(freq[-1], freq[0], 1, t_res/1000)
#       dm_step = 1/float(binshift)
#       dm_search = np.arange(dm_range[dm_bin] - 100*dm_step,dm_range[dm_bin] + 100*dm_step ,dm_step)
#       Find the dm that maximises SNR from 1st search    
        best_dm1 = dm_range[dm_bin]
#       A finer search around this dm
        if args.scatter:
            new_dm_search = np.linspace(-best_dm1, best_dm1, 200)
        else:
            new_dm_search = np.linspace(best_dm1-0.5*best_dm1, best_dm1+0.5*best_dm1, 200)
#       dm_search = np.linspace(dm_range[dm_bin - 1], dm_range[min(dm_bin + 1,19)], 200) # search for the range to try
        dm_range_2 = new_dm_search # set a new range
#       dm_range = np.linspace(-.025, .025, 200)
#       Plot the region we search for best dm
        if args.getPlot:
            # Create a snr vs dm plot for visualization
            snrfig1 = plt.figure()
            plt.plot(dm_range, peaks_of_average, '.')
            plt.axvline(x=best_dm1, color='r', ymax=np.max(peaks_of_average))
#        plt.fill_betweenx(peaks_of_average, x1=-best_dm1, x2=best_dm1, color='grey', alpha='0.5')
            plt.title('Best dm trial', fontsize=18)
            plt.xticks(fontsize=18)
            plt.yticks(fontsize=18)
            plt.xlabel('Delta dm (pc cm^-3)', fontsize=18)
            plt.ylabel('Average profile peak', fontsize=18)
            if args.doHC:
                snrfig1.savefig('SNR_DM_HC_seed_%d_1.png' %(iseed))
                snrfig1.clear()
                plt.close(snrfig1)
            else:
                snrfig1.savefig('SNR_DM_KJ07_seed_%d_1.png' %(iseed))
                snrfig1.clear()
                plt.close(snrfig1)

#----------------- SECOND ITERATION ----------------------------
#       Search for a best dm
        peaks_of_average = []
        for dm_id in range(len(dm_range_2)):
            shifted_profiles = []
            for freq_id in range(nch):
                bin_shift = bm.delay(freq[nch - 1], freq[freq_id], dm_range_2[dm_id], t_res/1000.) # Res increased by 1000 more bins
                shifted_profiles.append(np.roll(resampled[freq_id], bin_shift))
            average = bm.avg_prof(shifted_profiles)
            peaks_of_average.append(bm.find_peak(average))
        
#       Find the best dm (dm that maximises SNR)
        for i in range(len(peaks_of_average)):
            if peaks_of_average[i] == np.max(peaks_of_average):
                best_dm = dm_range_2[i]
        if args.getPlot:
#       Create a snr vs dm plot for visualization
            snrfig2 = plt.figure()
            plt.plot(dm_range_2, peaks_of_average, '.')
            if best_dm:
                plt.axvline(x=best_dm, color='r', ymax=np.max(peaks_of_average))
                #plt.annotate('annotate', xy=(best_dm, peaks_of_average[i]), xytext=(np.min(peaks_of_average), 0.01),\
                #arrowprops=dict(facecolor='black', shrink=0.05))
            plt.title('Best dm trial', fontsize=18)
            plt.xlabel('Delta dm (pc cm^-3)', fontsize=18)
            plt.ylabel('Average profile peak', fontsize=18)
            plt.yticks(fontsize=18)
            plt.xticks(fontsize=18)
            if args.doHC:
                snrfig2.savefig('SNR_DM_HC_seed_%d_bestDM_%.3f_2.png' %(iseed, best_dm))
                snrfig2.clear()
                plt.close(snrfig2)
            else:
                snrfig2.savefig('SNR_DM_KJ07_seed_%d_bestDm_%.3f_2.png' %(iseed, best_dm))
                snrfig2.clear()
                plt.close(snrfig2)
        
        # Write out important parameters into a file    
        if args.outfile:
            if args.doHC:
                model = 'Hollow_cone'
            else:
                model = 'Patchy_beam'
            pulsarParams = np.asarray([model, min_freq, chbw, nch, P, alpha, beta, \
                                      iseed, int(rand_dm), best_dm])
            f = open('dm_dat.txt', 'a')
            f.write(' '.join([str(item) for item in pulsarParams]) + ' \n')
else:
    template = profile[-1]
    lag_time = bm.cross_correlate(profile, template, period=P)
    if snr:
        sigma_lag = w10_noiseless/float(snr)
    else:
        sigma_lag = w10_noiseless
    dm_guess = (lag_time[-1] - lag_time[0]) / (4.148808 * 1e3 * ((freq[-1] * 1e3) ** (-2) - (freq[0] * 1e3)  ** (-2)))
    C_guess = 0.0 # Initial constant to shift the time delay across the x-axis
    popt, pcov = curve_fit(bm.dispersive_delay, freq * 1e3, lag_time, p0=[dm_guess, C_guess], sigma=sigma_lag)
    residual = lag_time - bm.dispersive_delay(freq * 1e3, popt[0], popt[1])
    if args.diagnostic:
        fig, ax = plt.subplots(1, 2, figsize=(16, 5))
        ax[0].plot(lag_time, freq, '.', label='Delta t')
        ax[0].plot(bm.dispersive_delay(freq * 1e3, popt[0], popt[1]), freq, '--', label='fit')
        ax[0].set_ylabel('Frequency (GHz)', fontsize=14)
        ax[0].set_xlabel(r'$\delta_t$ (s)', fontsize=14)
        ax[0].legend(fontsize=14)
        ax[1].plot(residual, '*', label='residual')
        #ax[1].errorbar(lag_time, freq, xerr=sigma_lag, ls='--', label=r'$\sigma_{TOA}$')
        ax[1].legend()
        ax[1].set_ylabel('Frequency (GHz)', fontsize=14)
        ax[1].set_xlabel(r'$\delta_t$ (s)', fontsize=14)

        fig.savefig('template_matching_%i.png' %iseed)
    # Write out important parameters into a file    
    if args.outfile:
        if args.doHC:
            model = 'Hollow_cone'
        else:
            model = 'Patchy_beam'
        pulsarParams = np.asarray([model, min_freq, chbw, nch, P, alpha, beta, \
        iseed, int(rand_dm), popt[0]])
        f = open('dm_dat_template_matching.txt', 'a')
        f.write(' '.join([str(item) for item in pulsarParams]) + ' \n')
        print lag_time, sigma_lag
#===========================================================================================================================
#                     PRODUCE PLOTS:
#===========================================================================================================================
if args.getPlot:
#===========================================
# SET A ZERO BASELINE AND PLOT THE PROFILE:
#===========================================
    fig2, ax2 = plt.subplots()
    for k in np.arange(len(profile)):
        ax2.plot(phase, profile[k] + k, label='frequency = %0.2f GHz' %freq[k], color='grey')
	plt.title('Profile evolution with frequency', fontsize=18)
	plt.xlabel('Phase (degrees)', fontsize=18)
        if not args.freqrange:
            plt.ylabel('Frequency (%.1f - %.1f MHz)' %(freq[0]*1e3, freq[-1]*1e3), fontsize=18)
        else:
            if len(args.freqrange) > 1:
                plt.ylabel('Frequency (%.1f - %.1f MHz)' %(args.freqrange[0], args.freqrange[-1]), fontsize=18)
            else:
                plt.ylabel('Frequency (%s MHz)' %(args.freqrange[0]), fontsize=18)
        plt.xticks(fontsize=18) 
        plt.yticks(fontsize=18) 
        plt.xlim(-180, 180)
        plt.tick_params(axis='y', which='both', left='off', top='off', labelleft='off')
#        plt.grid('on', which='both')
        #============================================
        #    2D emission region:
        #============================================
        xlos, ylos, thetalos = bm.los(alpha, beta, res)
        fig3 = plt.figure(figsize=(14,6))
        ax31 = fig3.add_subplot(1,2,1)
        plt.plot(xlos, ylos, '-', lw=2)
        plt.imshow(beam[k].T, extent=[-180, 180, 180, -180], cmap='gist_heat')
        plt.xlabel('X (degrees)', fontsize=18)
        #plt.title('Beam', fontsize=18)
        plt.ylabel('Y (degrees)', fontsize=18)
        plt.tick_params(axis='both', which='major', labelsize=18)
        # find zoomed extent fior plot
        nonZero = np.where(beam[0]!= 0.0)    # Set the maximum zoom to beam at lowest frequency
        nzx_min = np.min(np.amin(nonZero,0))
        nzy_min = np.min(np.amin(nonZero,1))
        nzx_max = np.max(np.amax(nonZero,0))
        nzy_max = np.max(np.amax(nonZero,1))
        x1 = phase[nzx_min]
        x2 = phase[nzx_max]
        y1 = phase[nzy_min]
        y2 = phase[nzy_max]
        plt.xlim(-20,20)
        plt.ylim(-20,20)
#        plt.xlim(x1,x2)
#        plt.ylim(y1,y2)
        cb = plt.colorbar()
        #cb.set_label('Intensity', labelpad=-40, y=0.45)
        cb.set_label('Intensity', size=18)
        cb.ax.tick_params(labelsize=18) 
        # 1D plot
        ax32 = fig3.add_subplot(1,2,2)
        #plt.tight_layout()
        plt.plot(phase, profile[0], color='grey', label='Scattered')
        # For comparision, overplot a non scattered profile with added noise
        if args.scatter:
            if args.snr:
                prof_noise = bm.add_noise(prof[0], rms, res)
                plt.plot(phase, prof_noise, alpha=0.5, label='Intrinsic')
            else:
                plt.plot(phase, prof[0], alpha=0.5, label='Intrinsic')
            plt.legend()
        plt.title('Profile at %.3f GHz' % freq[0], fontsize=18)
        plt.ylim(0,10)
        plt.xlim(-180, 180)
        plt.xlabel('Phase ', fontsize=18)
        plt.tight_layout()
        #plt.ylabel('Intensity', fontsize=18)
        plt.xticks(fontsize=16)
        plt.yticks(fontsize=16)
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
            plt.tight_layout()
            plt.tick_params(axis='both', which='major', labelsize=18)
            plt.plot(xlos, ylos, '+r')
            plt.imshow(beam[bid].T, extent=[-180, 180, 180, -180])
            plt.xlabel('X (degrees)', fontsize=18)
            plt.title('Radio pulsar beam', fontsize=18)
            plt.ylabel('Y (degrees)', fontsize=18)
            #plt.xticks(fontsize=18)
            #plt.yticks(fontsize=18)
            #plt.tick_params(axis='y', which='both', left='off', top='off', labelleft='off')
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
            cb = plt.colorbar()
            #cb.set_label('Intensity', labelpad=-40, y=0.45)
            #cb.set_label('Intensity', size=18)
            cb.ax.tick_params(labelsize=18)
            # 1D profile plot
            ax42 = fig4.add_subplot(1,2,2)
            plt.plot(phase, profile[bid], label='Scattered')
            # Overplot non scattered profile (noise added) for comparision
            if args.scatter:
                if args.snr:
                    prof_noise = bm.add_noise(prof[bid], rms, res)
                    plt.plot(phase, prof_noise, color='grey', alpha=0.5, label='Intrinsic')
                else:
                    plt.plot(phase, prof[bid], color='grey', alpha=0.5, label='Intrinsic')
                plt.legend()
            plt.title('Profile at %.3f GHz' % freq[bid], fontsize=18)
            plt.xlim(-180, 180)
            plt.ylim(np.min(profile), np.max(profile))
            plt.xlabel('Phase ', fontsize=18)
            plt.ylabel('Intensity', fontsize=18)
            plt.xticks(fontsize=18)
            plt.yticks(fontsize=18)
            if args.doHC:
                suffix = 'HC'
            else:
                suffix = 'KJ07'
            fig4.savefig('beam_%.3f_GHz_%s_%d_.png' %(freq[bid], suffix, int(iseed))) 
            fig4.clear()
            plt.close(fig4)
# Clean up 
#cwd = os.getcwd()
#folder = '/' + str(iseed)
#new_directory = cwd + folder
#os.mkdir(new_directory)
#files_to_move = '*%s*.png' %(str(iseed))
#os.system('mv %s %s' %(files_to_move, new_directory))
