#!/usr/bin/env python

import numpy as np
import pickle
import matplotlib.pyplot as plt
import argparse
import os, sys, time


#======================================================
#    DEFINE FUNCTIONS NEEDED:
#======================================================
#====================================================================================================================================================
#                                                       SCATTERING TIME:
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


#====================================================================================================================================================
#                                                       DM CURVE:
#====================================================================================================================================================
def find_phase_bin(prof):
    """Function to find a bin where the profile peaks.
       
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
    if phase_bin0==phase_bin1:
        delta_phase = phase[10]
    else:
        delta_phase = phase_at_peak1 - phase_at_peak0
    # Convert the phase to time and find a corresponding delta_dm:
    delta_t = delta_phase/360. * P
    D = 4.148808 * 1e3 # +/- 3e-6 MHz^2 pc^-1 cm^3 s
    dm = delta_t / (D * ((freq_ref * 1e3)**(-2) - (freq * 1e3)**(-2))) 
    delta_dm = np.linspace(-2*dm, 0, num=20)# try only 20 for now

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
#                                                       ADD NOISE:
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
    noisy_prof = prof + noise

    return noisy_prof


#====================================================================================================================================================================
#                                              INITIALIZE PARAMETERS FROM THE COMMAND LINE:
#====================================================================================================================================================================
parser = argparse.ArgumentParser(description='Find excess DM due to profile evolution with frequency. Uses the files produced when running gen_profile.py.', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('-iseed', metavar="<iseed>", type=int, default='4', help='integer seed for a pseudo-random number generator.')
parser.add_argument('-snr', metavar="<snr>", type=float, default=10000, help='signal to noise ratio.')
parser.add_argument('-dm', metavar="<dm>", type=float, default=1, help='dispersion measure in cm^-3 pc.')
parser.add_argument('-scatter', metavar="<0/1>", default=None, help='include scattering.')
parser.add_argument('-out', metavar="<outFile>", type=str, default="outFile", help='output file name.')
parser.add_argument('-outfile', metavar="<output file>", help="Write to file.")
args = parser.parse_args()
dm = args.dm
iseed = args.iseed
scr = args.scatter
snr = args.snr
outFile = args.out
pulsarParamsFile = args.outfile
pickleFile = pulsarParamsFile+'_p'
pickleFile2 = pulsarParamsFile+'_p2'
#======================================================
#    1. Load the files containing the profile and beam:
#======================================================

#prof_file = open('prof_data.txt', 'rb')
prof_file = open(pickleFile, 'rb')
#profDict = pickle.load(prof_file)
#prof = profDict['prof']
#phase = profDict['phase']
#beam = profDict['beam']
prof = pickle.load(prof_file)
phase = pickle.load(prof_file)
beam = pickle.load(prof_file)
prof_file.close()

params_file = open(pickleFile2, 'rb')
freq = pickle.load(params_file)
P = pickle.load(params_file)
iseed = pickle.load(params_file)
res = pickle.load(params_file)
t_res = pickle.load(params_file)
params_file.close()

nch = len(prof)

#======================================================
#     2. Scatter the line of sight profile: 
#======================================================
train = []
bf = []
tau = sc_time(freq, dm, iseed)
if scr == None:
    sc_prof = prof # profile without scattering

else:
    sc_prof = []
    for pid in range(nch):
        train.append(pulsetrain(3, res, prof[pid]))

    for fid in range(nch):
        tau = sc_time(freq[fid], dm, iseed)
        bf = broadening(tau, P, res)
        sc_train = scatter(train[fid], bf)
        sc_prof.append(extractpulse(sc_train, 2, res))

#======================================================
#     3. Add noise to the scattered profile:
#======================================================
peaks = []
for j in range(nch):
    peaks.append(find_peak(sc_prof[j]))

if snr == None:
    profile = sc_prof
else:
    rms = noise_rms(snr, np.max(peaks))
    profile = add_noise(sc_prof, rms, iseed, res)

#======================================================
#      4. Fit a DM Curve:
#======================================================

average_profile = []
peaks_of_average = []
phase_bin0 = find_phase_bin(profile[nch - 1])
phase_bin1 = find_phase_bin(profile[0])
dm_range = find_delta_dm(P, profile, phase, phase_bin0, phase_bin1, freq[nch - 1], freq[0], nch) 

# Create a directory to save the files.
#date = time.strftime("%Y-%m-%d")
date = time.ctime(time.time())
#newDir = os.mkdir(date)
#prevDir = os.getcwd()
#os.chdir(str(prevDir)/str(newDir))
for dm_id in range(len(dm_range)):
    shifted_profile = []
    #fig = plt.figure()
    #plt.grid()
    #plt.xlim(-180, 180)
    #plt.title("shifted profiles with respect to dm = " + str(dm_range[dm_id]))
    #plt.xlabel("phase (degrees)")
    #plt.ylabel("Intensity")
    #print "Shifting the profiles with dm = " + str(dm_range[dm_id])
    for freq_id in range(nch):
        #print "frequency " + str(freq[freq_id])       
        bin_shift = delay(freq[nch - 1], freq[freq_id], dm_range[dm_id], t_res)
        shifted_profile.append(np.roll(profile[freq_id], bin_shift))
        plt.plot(phase, shifted_profile[freq_id])
    average_profile.append(avg_prof(shifted_profile))
    peaks_of_average.append(find_peak(average_profile[dm_id]))
 #   fig.savefig(outFile+'_'+str(time.time())+'.png')

# Make an animation of the files to see the dm_shift.
#os.system("convert -delay 50 -loop 0 outFile*.png animation.gif")

#print "I------------------------------------------------------------I"
#print "I                 EXCESS DM                                  I"
#print "I------------------------------------------------------------I"
#print "Dm range in bins = " + str(dm_range/P * 360)
#print "Profile at minimum frequency will be shifted w.r.t profile at maximum frequency in " + str(delay(freq[nch - 1], freq[0], dm_range[0], t_res)/P * 360) + " bins"
#plt.figure()
#plt.grid()
#plt.xlim(-180, 180)
#plt.title("average profile")
#plt.xlabel("phase (degrees)")
#lt.ylabel("Intensity")
#peaks_of_average = []
#for i in np.arange(len(average_profile)):
#    plt.plot(phase, average_profile[i])
#    peaks_of_average.append(find_peak(average_profile[i]))

for i in range(len(peaks_of_average)):
    if peaks_of_average[i] == np.max(peaks_of_average):
       best_dm = dm_range[i]
f = open(pulsarParamsFile, 'a')
f.write(str(best_dm) + '\n')
f.close()
       #print "Best dm = " +  str(best_dm) + " pc cm^-3 " 
       #print "Best average profile at index " + str(i)
       #print "Highest peak of average profile " + str(peaks_of_average[i])    
       #print "\n"
       #print "I-----------------------------------------------------------------I"
       #plt.plot(phase, average_profile[i])
       #plt.show()
#shifted_with_best_dm = shifted_profile[450:500]
#print shifted_with_best_dm
#plt.figure()
#plt.xlim(-180,180)
#plt.title("Shifted profile with excess dm = %.5f pc cm^-3" %best_dm)
#for i in np.arange(len(shifted_with_best_dm)):
#    plt.plot(phase, shifted_with_best_dm)


#======================================================
#      5. Plot and visualise the data:
#======================================================
#    2D emission region:

#xlos, ylos, thetalos = los(alpha, beta, res)
#plt.figure(figsize=(10,5))
#plt.subplot(1, 2, 1)
#plt.plot(xlos, ylos, '--r')
#plt.imshow(beam[0], extent=[-np.amax(beam[0]),np.amax(beam[0]),-np.amax(beam[0]),np.amax(beam[0])]), cmap=cm.gray)
#plt.title('Patchy emission')
#plt.xlabel('X (degrees)')
#plt.ylabel('Y (degress)')
#plt.colorbar()
#plt.subplot(1, 2, 2)
#plt.plot(phase, averageP[0]) # average profile using first DM?
#plt.legend()
#plt.xlim(-180, 180)
#plt.title('Average pulse profile')

# Sequence of pulse profiles:

#plt.figure()
#plt.grid()
#plt.title('A sequence of %i pulse profile' %len(sc_prof))
#plt.xlabel('phase (degrees)')
#plt.xlim(-180, 180)
#plt.ylabel('profile number')

#for pid in np.arange(len(sc_prof)):
#    plt.plot(phase, sc_prof[pid])

# 2D Emission:

#plt.figure()
#plt.title('Beam')
#plt.xlabel('phase bin number')
#plt.ylabel('channels')
#plt.imshow(prof, aspect='auto', origin='lower')

#=========================================================
#   
#=========================================================
#plt.show()
