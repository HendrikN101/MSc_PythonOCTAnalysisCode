# -*- coding: utf-8 -*-
"""
Created on Mon May 13 11:31:25 2024

@author: hnie464



This Python script is for analysing individual A-scans.

Note: Zero-padding is implemented.


"""

' ~ Importing Required Modules ~ '

import numpy as np
import matplotlib.pyplot as plt
from nptdms import TdmsFile

' ~ Function to import TDMS Files ~ '

def read_1d_tdms(file_path):
    with TdmsFile.open(file_path) as tdms_file:
        group = tdms_file.groups()[0]
        channel = group.channels()[0]
        data = channel[:]
        return data



# ' ~ Importing TDMS Spectral Data ~ '

# file_path1 = r"C:\80k Data\Projects\AxialConversion\RefArm\RefArm_Mirror_0.tdms"
# file_path2 = r"C:\80k Data\Projects\AxialConversion\RefArm\RefArm_Mirror_10um.tdms"
# file_path3 = r"C:\80k Data\Projects\AxialConversion\RefArm\RefArm_Mirror_20um.tdms"
# file_path4 = r"C:\80k Data\Projects\AxialConversion\RefArm\RefArm_Mirror_30um.tdms"
# file_path5 = r"C:\80k Data\Projects\AxialConversion\RefArm\RefArm_Mirror_40um.tdms"
# file_path6 = r"C:\80k Data\Projects\AxialConversion\RefArm\RefArm_Mirror_50um.tdms"

# Spectrum_data1 = read_1d_tdms(file_path1)*np.hanning(4096)
# Spectrum_data2 = read_1d_tdms(file_path2)
# Spectrum_data3 = read_1d_tdms(file_path3)
# Spectrum_data4 = read_1d_tdms(file_path4)
# Spectrum_data5 = read_1d_tdms(file_path5)
# Spectrum_data6 = read_1d_tdms(file_path6)

zp = 2**14
# FFTdata1 = np.abs(np.fft.fft(Spectrum_data1, zp))
# FFTdata2 = np.abs(np.fft.fft(Spectrum_data2, zp))
# FFTdata3 = np.abs(np.fft.fft(Spectrum_data3, zp))
# FFTdata4 = np.abs(np.fft.fft(Spectrum_data4, zp))
# FFTdata5 = np.abs(np.fft.fft(Spectrum_data5, zp))
# FFTdata6 = np.abs(np.fft.fft(Spectrum_data6, zp))
x_axis = np.linspace(0,4095,zp)

# plt.figure(1)
# plt.plot(x_axis, FFTdata1)
# plt.plot(x_axis, FFTdata2)
# plt.plot(x_axis, FFTdata3)
# plt.plot(x_axis, FFTdata4)
# plt.plot(x_axis, FFTdata5)
# plt.plot(x_axis, FFTdata6)
# plt.xlim([200,250])
# plt.ylim([0,2.5e5])
# plt.legend(['0 um', '10 um', '20 um', '30 um', '40 um', '50 um'])
# plt.title('Comparison of PSFs at Different Positions (Reference Arm, 0-50 um)')
# plt.xlabel('Pixels')
# plt.ylabel('Intensity (Linear)')

# max1 = np.argmax(FFTdata1[80:int(zp/2)])+80
# max2 = np.argmax(FFTdata6[80:int(zp/2)])+80
# px_diff = max2-max1
# conversion = 50/(px_diff/4)
# print('Pixel-to_Micrometre Conversion (RefArm, 0-50um):', round(conversion, 2), 'um/px')



# ' ~ Importing TDMS Spectral Data ~ '

# file_path1 = r"C:\80k Data\Projects\AxialConversion\SampleArm\SampArm_Mirror_0.tdms"
# file_path2 = r"C:\80k Data\Projects\AxialConversion\SampleArm\SampArm_Mirror_10um.tdms"
# file_path3 = r"C:\80k Data\Projects\AxialConversion\SampleArm\SampArm_Mirror_20um.tdms"
# file_path4 = r"C:\80k Data\Projects\AxialConversion\SampleArm\SampArm_Mirror_30um.tdms"
# file_path5 = r"C:\80k Data\Projects\AxialConversion\SampleArm\SampArm_Mirror_40um.tdms"
# file_path6 = r"C:\80k Data\Projects\AxialConversion\SampleArm\SampArm_Mirror_50um.tdms"

# Spectrum_data1 = read_1d_tdms(file_path1)
# Spectrum_data2 = read_1d_tdms(file_path2)
# Spectrum_data3 = read_1d_tdms(file_path3)
# Spectrum_data4 = read_1d_tdms(file_path4)
# Spectrum_data5 = read_1d_tdms(file_path5)
# Spectrum_data6 = read_1d_tdms(file_path6)

# FFTdata1 = np.abs(np.fft.fft(Spectrum_data1, zp))
# FFTdata2 = np.abs(np.fft.fft(Spectrum_data2, zp))
# FFTdata3 = np.abs(np.fft.fft(Spectrum_data3, zp))
# FFTdata4 = np.abs(np.fft.fft(Spectrum_data4, zp))
# FFTdata5 = np.abs(np.fft.fft(Spectrum_data5, zp))
# FFTdata6 = np.abs(np.fft.fft(Spectrum_data6, zp))

# plt.figure(2)
# plt.plot(x_axis, FFTdata1)
# plt.plot(x_axis, FFTdata2)
# plt.plot(x_axis, FFTdata3)
# plt.plot(x_axis, FFTdata4)
# plt.plot(x_axis, FFTdata5)
# plt.plot(x_axis, FFTdata6)
# plt.xlim([200,250])
# plt.ylim([0,2.75e5])
# plt.legend(['0 um', '10 um', '20 um', '30 um', '40 um', '50 um'])
# plt.title('Comparison of PSFs at Different Positions (Sample Arm, 0-50 um)')
# plt.xlabel('Pixels')
# plt.ylabel('Intensity (Linear)')

# max1 = np.argmax(FFTdata1[80:int(zp/2)])+80
# max2 = np.argmax(FFTdata6[80:int(zp/2)])+80
# px_diff = max2-max1
# conversion = 50/(px_diff/4)
# print('Pixel-to_Micrometre Conversion (SampArm, 0-50um):', round(conversion, 2), 'um/px')



' ~ Importing TDMS Spectral Data ~ '

file_path1 = r"C:\80k Data\Projects\AxialConversion\Trial2\RefArm_100um\RefArm_Mirror_0um.tdms"
file_path2 = r"C:\80k Data\Projects\AxialConversion\Trial2\RefArm_100um\RefArm_Mirror_100um.tdms"
file_path3 = r"C:\80k Data\Projects\AxialConversion\Trial2\RefArm_100um\RefArm_Mirror_200um.tdms"
file_path4 = r"C:\80k Data\Projects\AxialConversion\Trial2\RefArm_100um\RefArm_Mirror_300um.tdms"
file_path5 = r"C:\80k Data\Projects\AxialConversion\Trial2\RefArm_100um\RefArm_Mirror_400um.tdms"
file_path6 = r"C:\80k Data\Projects\AxialConversion\Trial2\RefArm_100um\RefArm_Mirror_500um.tdms"

Spectrum_data1 = read_1d_tdms(file_path1)
Spectrum_data2 = read_1d_tdms(file_path2)
Spectrum_data3 = read_1d_tdms(file_path3)
Spectrum_data4 = read_1d_tdms(file_path4)
Spectrum_data5 = read_1d_tdms(file_path5)
Spectrum_data6 = read_1d_tdms(file_path6)

FFTdata1 = np.abs(np.fft.fft(Spectrum_data1, zp))
FFTdata2 = np.abs(np.fft.fft(Spectrum_data2, zp))
FFTdata3 = np.abs(np.fft.fft(Spectrum_data3, zp))
FFTdata4 = np.abs(np.fft.fft(Spectrum_data4, zp))
FFTdata5 = np.abs(np.fft.fft(Spectrum_data5, zp))
FFTdata6 = np.abs(np.fft.fft(Spectrum_data6, zp))

FFTList = [FFTdata1, FFTdata2, FFTdata3, FFTdata4, FFTdata5, FFTdata6]
    

signal_start_index = 75   # Index where the signal region starts on the x-axis (for automated plotting)
signal_region = FFTdata1[signal_start_index:250]   # Signal region
noise_region = FFTdata1[600:int(zp/3.5)]   # Noise region

# Calculate mean amplitude in the noise region and get max amplitude from signal region
max_signal = np.max(np.abs(signal_region))
index_max_signal = np.argmax(np.abs(signal_region))
stdDev_noise = np.std(noise_region)

R_samp = 1   # Reflectance of sample, assume 1 for mirror.
T_filt = 0.01441   # Transmittance of density filter - Change to correct value, look up optical density-to-transmittance conversion.
#snr_db_upper = (20 * np.log10(max_signal / (stdDev_noise + stdDev_error))) - (10 * np.log10(R_samp * (T_filt)**2))

snr_db_avg = (20 * np.log10(FFTList / stdDev_noise)) - (10 * np.log10(R_samp * (T_filt)**2))


'Plotting'

# plt.figure(3)
# plt.plot(x_axis, snr_db_avg[0,:])
# plt.plot(x_axis, snr_db_avg[1,:])
# plt.plot(x_axis, snr_db_avg[2,:])
# plt.plot(x_axis, snr_db_avg[3,:])
# plt.plot(x_axis, snr_db_avg[4,:])
# plt.plot(x_axis, snr_db_avg[5,:])
# plt.xlim([10,350])
# plt.ylim([0,110])
# plt.legend(['0 um', '100 um', '200 um', '300 um', '400 um', '500 um'])
# plt.title('Comparison of PSFs at Different Positions (Reference Arm, 0-500 um)')
# plt.xlabel('Pixels')
# plt.ylabel('Intensity (dB)')


plt.figure(3, figsize=(7,4))
plt.plot(x_axis, np.abs(FFTdata1))
plt.plot(x_axis, np.abs(FFTdata2))
plt.plot(x_axis, np.abs(FFTdata3))
plt.plot(x_axis, np.abs(FFTdata4))
plt.plot(x_axis, np.abs(FFTdata5))
plt.plot(x_axis, np.abs(FFTdata6))
plt.xlim([10,350])
plt.ylim([0,0.8e6])
plt.legend(['0 um', '100 um', '200 um', '300 um', '400 um', '500 um'])
plt.title('Comparison of PSFs at Different Positions (Reference Arm, 0-500 um)')
plt.xlabel('Pixels')
plt.ylabel('Intensity (Linear)')
plt.grid('True')


max1 = np.argmax(FFTdata1[10:int(zp/2)])+10
max2 = np.argmax(FFTdata6[10:int(zp/2)])+10
px_diff = max2-max1
conversion = 500/(px_diff/4)   # Dividing by 4 to account for zero-padding.
print('Pixel-to_Micrometre Conversion (RefArm, 0-500um):', round(conversion, 2), 'um/px')




print('Note: Micrometre has an uncertainty of 10 um, so reference arm 0-500 um value is more accurate.')