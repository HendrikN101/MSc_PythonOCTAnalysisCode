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


file_path = r"C:\80k Data\SNR\10-07-2024\5000Bscan_NoDispComp_BGsub.tdms"
#file_path = r"C:\80k Data\SNR\01-05-2024\Bscan_5000_LinDisp_exp30_BGsub_OD2.tdms"
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


' ~ Signal-to-Noise Ratio ~ '

# Define signal and noise regions
signal_start_index = 75   # Index where the signal region starts on the x-axis (for automated plotting)
signal_region = ascan[signal_start_index:250]   # Signal region
noise_region = ascan[600:int(N_size/3.5)]   # Noise region

# Calculate mean amplitude in the noise region and get max amplitude from signal region
max_signal = np.max(np.abs(signal_region))
index_max_signal = np.argmax(np.abs(signal_region))
stdDev_noise = np.std(noise_region)
#stdDev_error = np.std(noise_region, ddof=1) / np.sqrt(np.size(noise_region))

# Calculate SNR in dB
R_samp = 1   # Reflectance of sample, assume 1 for mirror.
T_filt = 0.01441   # Transmittance of density filter - Change to correct value, look up optical density-to-transmittance conversion.
#snr_db_upper = (20 * np.log10(max_signal / (stdDev_noise + stdDev_error))) - (10 * np.log10(R_samp * (T_filt)**2))
snr_db_avg = (20 * np.log10(ascan / stdDev_noise)) - (10 * np.log10(R_samp * (T_filt)**2))
#snr_db_lower = (20 * np.log10(max_signal / (stdDev_noise - stdDev_error))) - (10 * np.log10(R_samp * (T_filt)**2))

#print("Signal-to-Noise Ratio (SNR):", round(snr_db_avg, 3), "dB")



' ~ Plotting ~ '

plt.figure(figsize=(10, 6))
plt.plot(snr_db_avg, label='OCT Data')
plt.plot(index_max_signal + signal_start_index, snr_db_avg[169], marker='x', color='red', label='Maximum Signal')
plt.axvspan(75, 250, color='r', alpha=0.3, label='Signal Region')
plt.axvspan(600, int(N_size/3.5), color='g', alpha=0.3, label='Noise Region')
plt.xlabel('Depth (px)')
plt.ylabel('Amplitude (dB)')
plt.legend()
plt.grid(True)
plt.title('OCT A-scan Signal - SNR: 94.5 dB')
plt.show()

plt.figure()
PSFascan = np.abs(ascan[163:177])   # Make sure only exact signal peak with.
plt.plot(PSFascan)
#plt.xlim([40,60])
plt.title('Isolated Signal of Mirror - FWHM: 4.2 um in air')
plt.xlabel('Pixels (px)')
plt.ylabel('Amplitude (Linear)')

print('Axial Resolution (StDev):', round(np.std(PSFascan)/1e5, 1), 'um in air.')

halfMax = PSFascan.max()/2
plt.axhline(halfMax, c='r')
x_lower = 2.55
plt.axvline(x_lower, c='r')
x_upper = 10.05
plt.axvline(x_upper, c='r')

print('Axial Resolution (FWHM):', np.round(((x_upper-x_lower)*0.56), 1), 'microns in air,')

# plt.figure()
# plt.plot(np.abs(spectrum))