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
    print "bins of peaks: ", phase_bin0, phase_bin1
    D = 4.148808 * 1e3 # +/- 3e-6 MHz^2 pc^-1 cm^3 s
    range_dm = 0.01 # dm units
    if phase_bin0==phase_bin1:
        dm = 0
    #        delta_phase = phase[10] - phase[0]
    else:
        delta_phase = phase_at_peak1 - phase_at_peak0
        # Convert the phase to time and find a corresponding delta_dm:
        delta_t = delta_phase/360. * P
        dm = delta_t / (D * ((freq_ref * 1e3)**(-2) - (freq * 1e3)**(-2))) 
        print "find_delta_dm: ", dm, delta_t  
    delta_dm = np.linspace(-dm - range_dm, -dm + range_dm, num=20)# try only 20 for now
    print delta_dm
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
fileName =  args.outfile
beamFile = fileName+'_profiles'
paramsFile = fileName+'_params'
#======================================================
#    1. Load the files containing the profile and beam:
#======================================================

#prof_file = open('prof_data.txt', 'rb')
prof_file = open(beamFile, 'rb')
#profDict = pickle.load(prof_file)
#prof = profDict['prof']
#phase = profDict['phase']
#beam = profDict['beam']
prof = pickle.load(prof_file)
phase = pickle.load(prof_file)
beam = pickle.load(prof_file)
prof_file.close()

params_file = open(paramsFile, 'rb')
freq = pickle.load(params_file)
P = pickle.load(params_file)
iseed = pickle.load(params_file)
res = pickle.load(params_file)
t_res = pickle.load(params_file)
params_file.close()

nch = len(prof)


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
#newDir = os.mkdir(date)
#prevDir = os.getcwd()
#os.chdir(str(prevDir)/str(newDir))
for dm_id in range(len(dm_range)):
    shifted_profile = []
    '''
    fig = plt.figure()
    plt.grid()
    plt.xlim(-180, 180)
    plt.title("shifted profiles with respect to dm = " + str(dm_range[dm_id]))
    plt.xlabel("phase (degrees)")
    plt.ylabel("Intensity")
    '''
    #print "Shifting the profiles with dm = " + str(dm_range[dm_id])
    for freq_id in range(nch-1):
        #print "frequency " + str(freq[freq_id])       
        bin_shift = delay(freq[nch - 1], freq[freq_id], dm_range[dm_id], t_res)
        #print "frequency and shift: ", freq[freq_id], bin_shift
        shifted_profile.append(np.roll(profile[freq_id], bin_shift))
        #plt.plot(phase, shifted_profile[freq_id])
    average_profile.append(avg_prof(shifted_profile))
    peaks_of_average.append(find_peak(average_profile[dm_id]))

for i in range(len(peaks_of_average)):
    if peaks_of_average[i] == np.max(peaks_of_average):
       best_dm = dm_range[i]
f = open(pulsarParamsFile, 'a')
f.write(str(best_dm) + '\n')
f.close()

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
