# -*- coding: utf-8 -*-
"""
Created on Tue Mar  5 11:38:08 2024

@author: Hendrik



This Python script is for analysing individual B-scans.

Note: Zero-padding is implemented.


"""

' ~ Importing Required Modules ~ '

from nptdms import TdmsFile
import matplotlib.pyplot as plt
import numpy as np



' ~ Importing Linearisation & Dispersion Compensation Files ~ '

lin_txt = open(r'C:\80k Data\Lin_Disp_Files\From Python\lin_80k.txt', 'r')
lin_array = np.array(lin_txt.read().split())
lin_txt.close()
non_lin_txt = open(r'C:\80k Data\Lin_Disp_Files\From Python\non_lin_80k.txt', 'r')
non_lin_array = np.array(non_lin_txt.read().split())
non_lin_txt.close()
disp_comp_txt = open(r'C:\80k Data\Lin_Disp_Files\From Python\disp_comp_80k.txt', 'r')
disp_comp_array = np.array(disp_comp_txt.read().split())
disp_comp_txt.close()

def convert_strings_to_floats(input_array):
    output_array = []
    for element in input_array:
        converted_float = float(element)
        output_array.append(converted_float)
    return output_array

lin_array = convert_strings_to_floats(lin_array)
non_lin_array = convert_strings_to_floats(non_lin_array)
disp_comp_array = np.array(convert_strings_to_floats(disp_comp_array))



' ~ Function to import TDMS Files ~ '

def read_1d_tdms(file_path):
    with TdmsFile.open(file_path) as tdms_file:
        group = tdms_file.groups()[0]
        channel = group.channels()[0]
        data = channel[:]
        return data



' ~ Importing TDMS B-scan Data ~ '

file_path = r"C:\80k Data\TDMS Data\Bscans\Bscan_Mirror_Test.tdms"

Bscan_data = read_1d_tdms(file_path)
Bscan = Bscan_data.reshape((4096, 512), order='F')


' ~ Dispersion Compensation Loop ~ '

for n in range(0,512):
    corrected_spec = Bscan[:,n]*np.exp(1j*(disp_comp_array))
    Bscan[:,n] = corrected_spec


# plt.figure(frameon=False, figsize=(6, 6))
# plt.imshow(np.abs(Bscan), cmap='gray', aspect='auto')   # Raw 2D spectral image.
# plt.title('B-scan (Raw Spectra)', c='white')
# plt.axis('off')


# ' ~ Background Subtraction ~ '

# Background = np.mean(Bscan, axis = 1)
# Bscan = Bscan - Background[:,None]


' ~ B-scan Intensity Image ~ '

N_size = 2**14   # Zero-padding scale
FFT_Bscan = np.fft.fft(Bscan, axis=0)
FFT_Bscan = FFT_Bscan[0:2048, :]

plt.figure(figsize=(8, 6))
plt.imshow(np.abs(np.log(FFT_Bscan)), cmap='gray', aspect='auto', vmin=6, vmax=12)   # Fourier transformed 2D image.
plt.title('B-scan (No BG Sub.)')
plt.xlabel('A-scan Number (px)')
plt.ylabel('Depth (px)')
plt.axis('on')

# plt.figure()
# PSF = np.abs(FFT_Bscan[190:210,0])
# plt.plot(PSF)
# plt.title('Mirror PSF')
# plt.xlabel('Depth (px)')
# plt.ylabel('Amplitude (Linear)')
# halfMax = PSF.max()/2
# plt.axhline(halfMax, c='r')
# x_lower = 5.4
# plt.axvline(x_lower, c='r')
# x_upper = 13.8
# plt.axvline(x_upper, c='r')
# print('Axial Resolution:', round((x_upper - x_lower)*0.56, 2), 'um in air')



# ' Linearisation of the B-scan '

# Lin_Bscan_List = []
# for i in range(len(Bscan[0,:])):
#     interpolate = np.interp(lin_array, non_lin_array, Bscan[:,i])
#     Lin_Bscan_List.append(interpolate)
    
# Lin_Bscan = np.array(Lin_Bscan_List)


# ' ~ Dispersion Compensation Loop ~ '

# for i in range(len(Bscan[0,:])):
#     corrected_spec = Lin_Bscan[i,:]*np.exp(1j*(non_lin_array-lin_array))
#     Lin_Bscan[i,:] = corrected_spec

# FFT_Lin_Bscan = np.fft.fft(Lin_Bscan, axis=1)
# FFT_Lin_Bscan = FFT_Lin_Bscan[:, 0:2048]
# FFT_Lin_Bscan = np.transpose(FFT_Lin_Bscan, axes=(1,0))

# plt.figure(frameon=False, figsize=(8, 6))
# plt.imshow(np.abs(np.log(FFT_Lin_Bscan)), cmap='gray', aspect='auto')   # Fourier transformed 2D image.
# plt.title('B-scan (Compensated)', c='white')
# plt.axis('off')

# angleBscan = np.angle(FFT_Bscan[:,0])
# UWbscan = np.unwrap(angleBscan)
# plt.figure()
# plt.plot(UWbscan)   # Check to see if linear (mirror).

# plt.plot(FFT_Bscan[:,0]) # PSF should not disperse with depth if linearisation worked (mirror).
