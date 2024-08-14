# -*- coding: utf-8 -*-
"""
Created on Thu Jun 20 10:26:31 2024

@author: Hendrik
"""

from nptdms import TdmsFile
import matplotlib.pyplot as plt
import numpy as np
import tifffile as tf
import scipy.io
#import cv2
import colorsys
import matplotlib.colors as colors



' ~ Importing Linearisation & Dispersion Files ~ '


lin_txt = open(r'C:\80k Data\Lin_Disp_Files\From Python\lin_80k.txt', 'r')
lin_array = np.array(lin_txt.read().split())
lin_txt.close()
non_lin_txt = open(r'C:\80k Data\Lin_Disp_Files\From Python\non_lin_80k.txt', 'r')
non_lin_array = np.array(non_lin_txt.read().split())
non_lin_txt.close()
disp_comp_txt = open(r'C:\80k Data\Lin_Disp_Files\From Python\disp_comp_80k.txt', 'r')
disp_comp_array = np.array(disp_comp_txt.read().split())
disp_comp_txt.close()


Disp_mat = scipy.io.loadmat(r'C:\80k Data\Lin_Disp_Files\From Python\Disp.mat')
Disp_data = np.array(Disp_mat['Disp'])


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

file_path = r"C:\80k Data\Projects\Agarose\Agarose\1%Agarose\Agarose2_SquareWave_6mm1cycl_10Vpp.tdms"

# Set C-scan Size
Width = 512
Depth = 512
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
        #corrected_spec = Array[i,:,n]*np.exp(1j*(disp_comp_array))   # Original disp comp.
        Array[i,:,n] = corrected_spec


' ~ Mean Background Subtraction ~ '

ArraySample = Array[int(np.ceil(np.shape(Array)[0]/2)),:,:]
Background = np.mean(ArraySample,axis = 1)
Background = np.repeat(Background[:,None],(np.shape(ArraySample)[1]),axis = 1)
Array = Array - Background[None,:,:]


' ~ Fourier Transform Loop ~ '

#ZeroPad = 2**14
#ArrayZeroPad = np.zeros([Depth,Length,Width])
for i in range(0,Depth):
    Bscan = Array[i,:,:]
    FFT_Bscan = np.fft.fft(Bscan, axis=0)
    Array[i,:,:] = FFT_Bscan
    
Array = Array[:,0:2048,:]



' ~ Phase Analysis ~ '

PhaseArray = np.zeros([Depth,int(Length/2),Width-1])
for i in range(0,Depth):
    Phase = np.angle(Array[i,:,:])   # Extract phase information.
    PhaseUW = np.unwrap(Phase)   # Unwrapped phase.
    PhaseDiff = np.diff(PhaseUW, axis=1)   # Phase difference of adjacent A-scans in the same B-scan.
    PhaseArray[i,:,:] = PhaseDiff   # Append to new array.



' ~ Plotting ~ '

plt.figure(1, figsize=(6, 8))
Preview_Bscan = np.log(np.abs(Array[:,:,100]))   # Plotting B-scan in the time axis.
plt.imshow(np.transpose(Preview_Bscan), cmap='gray', vmin=6, vmax=12)
plt.title('B-scan (t = 0)')
plt.xlabel('Width (px) - 3 mm')
plt.ylabel('Length (px) - 4 mm')

plt.figure(2, figsize=(6, 8))
Preview_Phase = PhaseArray[:,:,100]   # Plotting Phase B-scan in the time axis.
plt.imshow(np.transpose(Preview_Phase), cmap='seismic', norm=colors.CenteredNorm())
plt.title('Phase of B-scan (t = 0)')
plt.xlabel('Width (px) - 3 mm')
plt.ylabel('Length (px) - 4 mm')


plt.figure(3, figsize=(6, 8))
Preview_Bscan2 = np.log(np.abs(Array[:,:,150]))   # Plotting B-scan in the time axis.
plt.imshow(np.transpose(Preview_Bscan2), cmap='gray', vmin=6, vmax=12)
plt.title('B-scan (t = 1.25 ms)')
plt.xlabel('Width (px) - 3 mm')
plt.ylabel('Length (px) - 4 mm')

plt.figure(4, figsize=(6, 8))
Preview_Phase2 = PhaseArray[:,:,4]   # Plotting Phase B-scan in the time axis.
plt.imshow(np.transpose(Preview_Phase2), cmap='seismic', norm=colors.CenteredNorm())
plt.title('Phase of B-scan (t = 1.25 ms)')
plt.xlabel('Width (px) - 3 mm')
plt.ylabel('Length (px) - 4 mm')


surface = 480
plt.figure(5, figsize=(8, 8))
Preview_Phase_surface = PhaseArray[:,surface,:]   # Plotting Phase B-scan in the time axis.
plt.imshow(np.flipud(np.transpose(Preview_Phase_surface)), cmap='seismic', norm=colors.CenteredNorm())
plt.title('Phase at surface level (t = 1.25 ms)')
plt.xlabel('Width (px) - 3 mm')
plt.ylabel('Length (px) - 4 mm')



# ' ~ Building HSV Image & Converting to RGB~ '

# Vmin = 6
# Vmax = 12
# Int_Image_Norm = (Preview_Bscan - Vmin)/(Vmin + Vmax)   # Normalising Intensity Values between 0-1. TEST: np.uint8(255* ((Preview_Bscan - Vmin)/(Vmin + Vmax)))

# Sat_Image = np.ones((2048, 512))*0.8   # Creating a saturation channel of 0.8 values. (Not used in analysis).

# Pmin = 6   # Change these to phase values.
# Pmax = 12
# Phase_Image_Norm = ((Preview_Phase - Pmin)/(Pmin + Pmax))   # Normalising Phase Values between 0-1.

# HSV_Image = np.append(np.append(Phase_Image_Norm, Sat_Image), Int_Image_Norm)   # Building HSV Image.

# RBG_Image = colorsys.hsv_to_rgb(HSV_Image[0], HSV_Image[1], HSV_Image[2])   # Converting HSV to RGB.

# RGB_Image = np.uint8(255*RBG_Image)   # Convert RGB image to 8-bit to save memory.
