# -*- coding: utf-8 -*-
"""
Created on Tue Sep 26 12:25:20 2023

@author: Hendrik
"""

from __future__ import print_function 
import numpy as np
import matplotlib.pyplot as plt
import tifffile as tf   # External download
from scipy.fft import fft, ifft
from nptdms import TdmsFile
from scipy.stats import norm



' ~ Import Data as TDMS ~ '

def read_1d_tdms(file_path):
    with TdmsFile.open(file_path) as tdms_file:
        group = tdms_file.groups()[0]
        channel = group.channels()[0]
        data = channel[:]
        return data


file_path = r"C:\80k Data\SNR\26-04-2024\Bscan_5000_OD1_LinDisp_GalvoOff.tdms"
Bscan_data = read_1d_tdms(file_path)
Complex_Bscan = Bscan_data.reshape(4096, 5000, order='F')

# plt.imshow(np.abs(Complex_Bscan), cmap='gray', aspect='auto')   # Raw B-scan (Absolute Value)



' ~ Data Processing ~ '

# Performing Fourier transform on spectrum
N_size = 2**14   # Zero-padding scale
spectrum = Complex_Bscan[:,0]
spectrumFFT = np.fft.fft(spectrum, n=N_size)   # Peforms a Fourier Transform on Spectrum
ascan = spectrumFFT[0:int(N_size/2)]   # Takes first A-scan from B-scan
#depth = np.linspace(0, 2047, 2048)   # Depth = 4096/2



' ~ Phase Stability ~ '

bscanFFT = np.fft.fft(Complex_Bscan, axis=0)
bscanPhase = np.angle(bscanFFT)
bscanUWPhase = np.unwrap(bscanPhase)


# plt.figure()
# plt.plot(np.abs(np.fft.fft(bscanUWPhase[114,:])))
bscanUWPhaseFFT = np.fft.fft(bscanUWPhase[114,:])

plt.figure()
plt.plot(np.abs(bscanUWPhaseFFT))
plt.xlim([0,100])
plt.ylim([-20,150])
plt.title('Fourier Transform of Phase at Low Frequencies')
plt.xlabel('Frequency (Hz)')
plt.ylabel('Amplitude (Linear)')
plt.grid('True')


# plt.figure()
# plt.plot(bscanUWPhaseFFT)
bscanUWPhaseFFT[1:75] = 0   # Filtering out very low frequency noise.
bscanUWPhaseFFT[4910:5000] = 0   # Filtering out low frequency from mirrored signal.
# plt.figure()
# plt.plot(np.abs(bscanUWPhaseFFT))
bscanUWPhaseiFFT = np.fft.ifft(bscanUWPhaseFFT)
bscanUWPhaseiFFT[:50] = 2.5475   # Removing sinc function ringing along edges from previous filtering.
bscanUWPhaseiFFT[4950:] = 2.5475   # ""


bscanSubsetNoFilt = np.std(bscanUWPhase[114,:], axis=0)   # Phase stability measurement (Unfiltered).
print('Phase Stability (Unfiltered):', round(bscanSubsetNoFilt*1000,1), 'mrad')

bscanSubsetFilt = np.std(bscanUWPhaseiFFT)   # Phase stability measurement (Filtered).
print('Phase Stability (Filtered):', round(bscanSubsetFilt*1000000,1), 'urad')



# ' ~ Phase Stability (Old) ~ '

# bscanFFT = np.fft.fft(Complex_Bscan, axis=0)
# bscanPhase = np.angle(bscanFFT)
# bscanUWPhase = np.unwrap(bscanPhase)
# bscanDiff = np.diff(bscanUWPhase, axis=1)
# bscanDiffSubset = np.mean(bscanDiff[49:51,:], axis=0)
# #bscanDiff1D = bscanDiff.flatten()

# 'Histogram'
# y, x, _ = plt.hist(bscanDiffSubset, 100)
# #plt.xlim(-0.0025,0.0025)
# x_lower = -0.035
# plt.axvline(x_lower, c='r')
# x_upper = 0.035
# plt.axvline(x_upper, c='r')
# plt.title('Histogram of Phase Difference in Adjacent A-scans (Sub-set) - FWHM: 70 mrad')
# plt.ylabel('Amount in Bins')
# plt.xlabel('Phase Difference (Rad)')
# print('Phase Stability FWHM (Sub-set): ', round((x_upper-x_lower)*1000, 3), 'mRad')

# 'Gaussian Fit'
# mu, sigma = norm.fit(bscanDiffSubset)
# best_fit_line = norm.pdf(x, mu, sigma)
# plt.plot(x, best_fit_line*8)
# halfMax = (best_fit_line*8).max()/2
# plt.axhline(halfMax, c='r')
# plt.figure()

# y_full, x_full, _ = plt.hist(bscanDiff1D, 50)
# halfMax_full = y_full.max()/2
# plt.axhline(halfMax_full, c='r')
# x_lower_full = -0.75
# plt.axvline(x_lower_full, c='r')
# x_upper_full = 0.75
# plt.axvline(x_upper_full, c='r')
# plt.title('Histogram of Phase Difference in Adjacent A-scans')
# plt.ylabel('Amount in Bins')
# plt.xlabel('Phase Difference (Rad)')
# plt.figure()
# print('Phase Stability FWHM: ', round((x_upper_full-x_lower_full)*1000, 3), 'mRad')
# plt.imshow(20*np.log(np.abs(bscanFFT[0:2048])), cmap='gray', aspect='auto')   # Processed B-scan



' ~ Plotting ~ '

plt.figure()
plt.plot(bscanUWPhase[114,:])
plt.plot(bscanUWPhaseiFFT)
plt.legend(['Unfiltered','Filtered'])
plt.title('Phase Stability - Low Frequencies Filtered')
plt.xlabel('A-scan Number')
plt.ylabel('Phase (rad)')
plt.grid('True')

plt.figure()
plt.hist(bscanUWPhaseiFFT, 50)   # Gaussian distribution that agrees with the theory.
plt.title('Histogram of Filtered Unwrapped Phase')
plt.xlabel('Bins (50)')
plt.ylabel('Amount in Bins')

# halfMax = bscanUWPhaseiFFT.max()/2
# plt.axhline(halfMax, c='r')
# x_lower = 2.55
# plt.axvline(x_lower, c='r')
# x_upper = 10.05
# plt.axvline(x_upper, c='r')

# plt.figure()
# PSFascan = np.abs(ascan[160:180])
# plt.plot(PSFascan)
# #plt.xlim([40,60])
# plt.title('Isolated Signal of Mirror - FWHM: 4.2 um in air')
# plt.xlabel('Pixels (px)')
# plt.ylabel('Amplitude (Linear)')
# halfMax = PSFascan.max()/2
# plt.axhline(halfMax, c='r')
# x_lower = 5.5
# plt.axvline(x_lower, c='r')
# x_upper = 13
# plt.axvline(x_upper, c='r')
# print('Axial Resolution:', np.round(((x_upper-x_lower)*0.56), 1), 'microns in air,', np.round((((x_upper-x_lower)*0.56)/1.3), 1), 'microns in tissue.')
