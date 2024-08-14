# -*- coding: utf-8 -*-
"""
Created on Wed Feb 14 10:46:14 2024

@author: Hendrik
"""

from nptdms import TdmsFile
import matplotlib.pyplot as plt
import numpy as np

def read_1d_tdms(file_path):
    with TdmsFile.open(file_path) as tdms_file:
        # Assuming only one group and one channel in the TDMS file
        group = tdms_file.groups()[0]
        channel = group.channels()[0]
        
        # Extracting data from the channel
        data = channel[:]
        return data

# Example usage:
file_path = r"C:\80k Data\TDMS Data\Spectrums\FP_Spec4Lin.tdms"   # Fabry-Perot Interferometer spectra to avoid any dispersion.
data = read_1d_tdms(file_path)
print(data)


plt.plot(data)
plt.title('850 nm Spectrum')
plt.ylabel('Amplitude (12-bit)')
plt.xlabel('Line Camera Width (Pixels)')
plt.grid('True')


fourier_data = np.fft.fft(data)   # Fourier transform of the data
fourier_data[0:45] = 0   # Setting everything but the signal to zero
fourier_data[80:] = 0
inv_fourier_data = np.fft.ifft(fourier_data)   # Inverse Fourier transform of the modified data


phase = np.angle(inv_fourier_data)   # Extracting the phase of the modified data
non_lin = np.unwrap(phase)   # Unwrapping the phase to get the non-linear data [Extract to .txt]
plt.figure()
plt.plot(non_lin, label='Non-Linear')

' ~ Writes the Non-Linear curve to a .txt File [Uncomment] ~ '

# non_lin_txt = open(r'C:\80k Data\Lin_Disp_Files\From Python\non_lin_80k.txt', 'w')   # Putting the non-linear data into a .txt file

# for i in range(0, len(non_lin)):
#     non_lin_txt.write(str(non_lin[i]) + '\t')

# non_lin_txt.close()


m = (np.max(non_lin)-np.min(non_lin))/4095
axis = np.linspace(0,4095,4096)
linear = m * axis + np.min(non_lin)   # Forming the linear fit [Extract to .txt]
plt.plot(linear, label='Linear')
plt.legend()
plt.grid('True')

' ~ Writes the Linear curve to a .txt File [Uncomment] ~ '

# lin_txt = open(r'C:\80k Data\Lin_Disp_Files\From Python\lin_80k.txt', 'w')   # Putting the linear data into a .txt file

# for i in range(0, len(linear)):
#     lin_txt.write(str(linear[i]) + '\t') #  + '\n'

# lin_txt.close()

pre_spectrum = np.fft.fft(data)
interpolate = np.interp(linear, non_lin, data)   # Interpolating the linear, non-linear and original data
plt.figure()
plt.plot(np.log(np.abs(pre_spectrum)), label='Original A-scan')
plt.plot(np.log(np.abs(np.fft.fft(interpolate))), label='Linearised A-scan')
plt.xlim(0,2048)
plt.grid('True')
plt.legend()