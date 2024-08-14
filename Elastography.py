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
import scipy
from sklearn.linear_model import LinearRegression



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

file_path = r"C:\80k Data\Projects\OCE\AgaroseOCE_10Vpp_3x0mm_512x10.tdms"

# Set C-scan Size
Width = 512
Depth = 10
Length = 4096

wavelength = 842.5
ref_in = 1.32

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


' ~ Fourier Transform ~ '

Array = Array.astype(complex)
Array = np.fft.fft(Array, axis=1)
Array = Array[:,0:2048,:]
Array = np.roll(Array, -357, axis=2)


' ~ Phase Analysis ~ '

ArrayPhi = np.angle(Array)
phaseDif = np.diff(ArrayPhi, axis = 0)

PhaseDiff0 = scipy.ndimage.median_filter(phaseDif[0,:,:], 5)
PhaseDiff1 = scipy.ndimage.median_filter(phaseDif[1,:,:], 5)

'Averging every second phase image'
phaseDiff_avg = phaseDif[0:9:2]
phaseDiff_avg = np.mean(phaseDiff_avg, axis=0)
PhaseDiff_avg = scipy.ndimage.median_filter(phaseDiff_avg, 5)

ImgMask = np.ones_like(Array[0,:,:])
ImgMask = ImgMask.astype(float)
ImgMask = ImgMask*(scipy.ndimage.median_filter(np.log(np.abs(Array[0,:,:]))>6, 5))
#ImgMask = ImgMask*(scipy.ndimage.convolve(np.log(np.abs(Array[0,:,:]))>3, np.ones((7,7))/(7*7)))
ImgMask= ImgMask.astype(int)
#plt.imshow(ImgMask)


' ~ Plotting ~ '

plt.figure(1, figsize=(8, 6))
Preview_Bscan = np.log(np.abs(Array[0,:,:]))   # Plotting B-scan in the time axis.
plt.imshow(Preview_Bscan, cmap='gray', aspect='auto', vmin=5, vmax=12)
plt.title('Agarose - Intensity Image')
plt.xlabel('Width (px)')
plt.ylabel('Length (px)')
plt.ylim([1024, 0])

# plt.figure(2, figsize=(10, 6))
# Preview_Phase = PhaseDiff0*ImgMask[:,:]   # Plotting Phase B-scan in the time axis.
# Preview_Phase = (wavelength*Preview_Phase)/(4*np.pi*ref_in)   # Converting to displacement.
# plt.imshow(Preview_Phase, cmap='seismic', aspect='auto', norm=colors.CenteredNorm())
# plt.title('Agarose - Phase Difference (Pos.)')
# plt.xlabel('Width (px)')
# plt.ylabel('Length (px)')
# plt.ylim([1024, 0])
# plt.colorbar()


# # plt.figure(3, figsize=(8, 6))
# # Preview_Bscan2 = np.log(np.abs(Array[1,:,:]))   # Plotting B-scan in the time axis.
# # plt.imshow((Preview_Bscan2), cmap='gray', aspect='auto', vmin=5, vmax=12)
# # plt.title('B-scan (2)')
# # plt.xlabel('Width (px)')
# # plt.ylabel('Length (px)')
# # plt.ylim([1024, 0])

# plt.figure(3, figsize=(10, 6))
# Preview_Phase2 = PhaseDiff1*ImgMask[:,:]   # Plotting Phase B-scan in the time axis.
# Preview_Phase2 = (wavelength*Preview_Phase2)/(4*np.pi*ref_in)   # Converting to displacement.
# plt.imshow(Preview_Phase2, cmap='seismic', aspect='auto', norm=colors.CenteredNorm())
# plt.title('Piezo Surface - Phase Difference (Neg.)')
# plt.xlabel('Width (px)')
# plt.ylabel('Length (px)')
# plt.ylim([1024, 0])
# plt.colorbar()

plt.figure(4, figsize=(10, 6))
Preview_Phase_Avg = PhaseDiff_avg*ImgMask[:,:]   # Plotting Phase B-scan in the time axis.
Preview_Phase_Avg = (wavelength*Preview_Phase_Avg)/(4*np.pi*ref_in)   # Converting to displacement.
plt.imshow(Preview_Phase_Avg, cmap='seismic', aspect='auto', vmin=-150, vmax=150)
plt.title('Agarose - Phase Difference (Avg.)')
plt.xlabel('Width (px)')
plt.ylabel('Length (px)')
plt.ylim([1024, 0])
plt.colorbar()



avgPlot = []
for i in range(155,232):
    avgPlot.append(Preview_Phase_Avg[400:500, i])

avgPlot = np.mean(avgPlot, axis=0)
x = np.reshape(np.linspace(400, 499, 100), (-1,1))

model = LinearRegression().fit(x, avgPlot)
r_sq = model.score(x, avgPlot)

plt.figure()
plt.plot(x, avgPlot)
plt.plot(x, model.predict(x), color='r')
plt.title('Agarose Displacement Gradient (400-500px, Avg.)')
plt.xlabel('Depth (px)')
plt.ylabel('Displacement (nm)')
plt.legend(['Raw Gradient', 'Linear Fit (R2 = 0.885)'])
plt.grid('True')


Strain = avgPlot/(2.24*(np.linspace(1,100, 100))*1e3)
avgStrain = np.mean(Strain[20:100])

LimStrain = Strain[20:100]
LimX = x[20:100]

model2 = LinearRegression().fit(LimX, LimStrain)
r_sq2 = model.score(LimX, LimStrain)

plt.figure(figsize=(7,4))
plt.plot(x, Strain)
plt.plot(LimX, model2.predict(LimX), '--', color='r')
plt.title('Strain vs Depth (Agarose)')
plt.xlabel('Depth [px]')
plt.ylabel('Strain')
plt.legend(['Strain [400-500px]','Linear Fit [420-500px]'])
plt.grid('True')

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
