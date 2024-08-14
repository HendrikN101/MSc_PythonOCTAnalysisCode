# -*- coding: utf-8 -*-
"""
Created on Tue Mar  5 11:38:08 2024

@author: Hendrik



This Python script is for analysing M-scan volumes.

(!) WARNING: May require a decent amount of system memory (RAM) to operate.

Note: Zero-padding is implemented.


"""

' ~ Importing Required Modules ~ '

from nptdms import TdmsFile
import matplotlib.pyplot as plt
import numpy as np
import tifffile as tf
import scipy.io



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


Disp_mat = scipy.io.loadmat(r'C:\80k Data\Lin_Disp_Files\From Python\Disp.mat')
Disp_data = np.array(Disp_mat['Disp'])


' ~ Importing TDMS Volume Scan into a Python Array ~ '

file_path = r"C:\80k Data\Projects\DopplerFlow\Trial4\DopplerDoubleTube1024_1p5mm_5deg_800uLmin_700um.tdms"

# Set C-scan Size
Width = 1024
Depth = 4
Length = 4096

Array = np.zeros([Depth,Length,Width])   # [Depth, Length, Width]

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
        corrected_spec = Array[i,:,n]*np.exp(-1j*(np.squeeze(Disp_data)))   # Darven's disp comp.
        #corrected_spec = Array[i,:,n]*np.exp(1j*(disp_comp_array)) # [i,:,n]
        Array[i,:,n] = corrected_spec


' ~ Background Subtraction ~ '
ArraySample = Array[int(np.ceil(np.shape(Array)[0]/2)),:,:]
Background = np.mean(ArraySample,axis = 1)
Background = np.repeat(Background[:,None],(np.shape(ArraySample)[1]),axis = 1)
Array = Array - Background[None,:,:]


' ~ Fourier Transform Loop ~ '

ZeroPad = 2**14
ArrayZeroPad = np.zeros([Depth,Length*4,Width]) # Length*4
for i in range(0,Depth):
    Bscan = Array[i,:,:]
    FFT_Bscan = np.fft.fft(Bscan, n=ZeroPad, axis=0) #, n=ZeroPad, 
    ArrayZeroPad[i,:,:] = FFT_Bscan
    
Array = ArrayZeroPad[:,0:int(ZeroPad/2),:] # 0:int(ZeroPad/2)
Array = np.roll(Array,-780,axis = 2)  # Fix vertical displacement [ | ]


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



# ' ~ En-face Image ~ '
# arr = np.log(np.median(np.abs((Array[:,190:220,:]+1)), axis=1))
# arr = np.roll(arr,15,axis = 0)   # Fix horizontal displacement [ __ ]
# #arr = np.roll(arr,-250,axis = 1)  # Fix vertical displacement [ | ]
# plt.figure(figsize=(6, 6))
# plt.imshow(arr, cmap='gray', aspect='equal', vmin=7, vmax=11)   # Fourier transformed 2D image.
# plt.title('En-Face', c='black')
# plt.xlabel('A-scan # (2mm) [px]')
# plt.ylabel('B-scan # (2mm) [px]')
# #plt.axis('off')


' ~ B-scan Intensity Image ~ '

Bscan = Array[3,:,:]
#Bscan = np.roll(Bscan,-660,axis = 1)  # Fix vertical displacement [ | ]

plt.figure(frameon=False, figsize=(8, 6))
plt.imshow(np.log(np.abs(Bscan)), cmap='gray', aspect='auto', vmin=5, vmax=11)   # Fourier transformed 2D image.
plt.title('B-scan (Original)', c='black')
#plt.axis('off')


' ~ Saving to a Tiff ~ '

# tf.imwrite('C:\80k Data\TIFF Data\Droplet.tif', Array)


