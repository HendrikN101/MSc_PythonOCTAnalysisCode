# -*- coding: utf-8 -*-
"""
Created on Tue Mar  5 11:38:08 2024

@author: Hendrik



This Python script is for analysing full C-scan volumes.

(!) WARNING: Requires a lot of system memory (RAM) to operate [4.2 GB per 512x512 C-scan].

Note: Does NOT implement zero-padding to avoid crashing PC.

"""

' ~ Importing Required Modules ~ '

from nptdms import TdmsFile
import matplotlib.pyplot as plt
import numpy as np
import tifffile as tf



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

lin_array = np.array(convert_strings_to_floats(lin_array))
non_lin_array = np.array(convert_strings_to_floats(non_lin_array))
disp_comp_array = np.array(convert_strings_to_floats(disp_comp_array))



' ~ Importing TDMS Volume Scan into a Python Array ~ '

file_path = r"C:\80k Data\TDMS Data\Cscans\USAF Lat Res\USAF_C_512x512_90d2.tdms"

# Set C-scan Size
Width = 512
Depth = 512
Length = 4096

Array = np.zeros([Depth,Length,Width])   # [Depth, Length, Width]
Array = Array.astype(complex)

with TdmsFile.read(file_path) as tdms_file:   # TdmsFile.read writes straight to memory!
    for i in range(0,Depth):
        for n in range(0,Width):
              group = tdms_file["Buffer " + str('{0:03}'.format(i))]
              all_group_channels = group.channels()
              Array[i,:,n] = all_group_channels[n][:]


' ~ Linearisation Loop ~ '

for n in range(0,Width):
    for i in range(0,Depth):
        interpolate = np.interp(lin_array, non_lin_array, Array[i,:,n])
        Array[i,:,n] = interpolate


' ~ Dispersion Compensation Loop ~ '

for n in range(0,Width):
    for i in range(0,Depth):
        corrected_spec = Array[i,:,n]*np.exp(1j*(disp_comp_array))
        Array[i,:,n] = corrected_spec


' ~ Background Subtraction ~ '

ArraySample = Array[int(np.ceil(np.shape(Array)[0]/2)),:,:]
Background = np.mean(ArraySample,axis = 1)
Background = np.repeat(Background[:,None],(np.shape(ArraySample)[1]),axis = 1)
Array = Array - Background[None,:,:]


' ~ Fourier Transform Loop ~ '

for i in range(0,Depth):
    Bscan = Array[i,:,:]
    FFT_Bscan = np.fft.fft(Bscan, axis=0)
    Array[i,:,:] = FFT_Bscan
    
Array = Array[:,0:2048,:]


# ' ~ A-scan Averaging Fourier Transform Loop ~ '

# for i in range(0,Depth):
#     Bscan = Array[i,:,:]
#     FFT_Bscan = np.fft.fft(Bscan, axis=0)
#     Array[i,:,:] = FFT_Bscan
    
#     image_height = Array[:,0].size
#     image_width = Array[0,:].size
#     target_columns = 256
#     columns_per_group = image_width // target_columns
#     Ascan_averaged = Array[:target_columns*columns_per_group, :].reshape(image_height, target_columns, columns_per_group)
#     Ascan_averaged = np.mean(Ascan_averaged, axis=2)   # A-scan averaged phase difference.
    
# Array = Ascan_averaged[:,0:2048,:]


# ' ~ Fourier Transform Loop ~ '

# ArrayFT = np.zeros([Depth,int(Length/2),Width])   # [Depth, Length, Width]

# for i in range(0,Depth):
#     Bscan = Array[i,:,:]
#     FFT_Bscan = np.fft.fft(Bscan, axis=0)
#     FFT_Bscan = FFT_Bscan[0:2048, :]
#     ArrayFT[i,:,:] = FFT_Bscan



' ~ En-face Image ~ '
arr = np.log(np.median(np.abs((Array[:,390:410,:]+1)), axis=1))
arr = np.roll(arr,27,axis = 0)   # Fix horizontal displacement [ __ ]
arr = np.roll(arr,-444,axis = 1)  # Fix vertical displacement [ | ]
plt.figure(figsize=(6, 6))
plt.imshow(arr, cmap='gray', aspect='equal', vmin=6, vmax=12)   # Fourier transformed 2D image.
plt.title('En-Face of 1951 USAF Target - F230FC-850', c='black')
plt.xlabel('A-scan # (1mm) [px]')
plt.ylabel('B-scan # (1mm) [px]')
#plt.axis('off')


' ~ B-scan Intensity Image ~ '

Bscan = Array[100,:,:]

plt.figure(frameon=False, figsize=(8, 6))
plt.imshow(np.log(np.abs(Bscan)), cmap='gray', aspect='auto', vmin=5, vmax=11)   # Fourier transformed 2D image.
plt.title('B-scan (Original)', c='black')
#plt.axis('off')


' ~ Saving to a Tiff ~ '

# tf.imwrite('C:\80k Data\TIFF Data\Droplet.tif', Array)


