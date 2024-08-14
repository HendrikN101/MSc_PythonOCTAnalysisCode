# -*- coding: utf-8 -*-
"""
Created on Tue May  7 14:10:27 2024

@author: hnie464
"""

import numpy as np
import matplotlib.pyplot as plt
#import tifffile as tf   # External download
#from scipy.fft import fft, ifft
from scipy import signal
#from nptdms import TdmsFile   # External Download
#from PIL import Image
#from skimage import data, img_as_float, color, exposure
from skimage.restoration import unwrap_phase
from scipy.signal import convolve
#import matplotlib as mpl
from matplotlib.colors import ListedColormap
import matplotlib.colors as colors
from scipy.signal import butter,filtfilt
from nptdms import TdmsFile
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


' Function to import TDMS Files '

def read_1d_tdms(file_path):
    with TdmsFile.open(file_path) as tdms_file:
        group = tdms_file.groups()[0]
        channel = group.channels()[0]
        data = channel[:]
        return data



' Setting Variables [Change to match your data] '

Width = 1024   # Number of A-scans per B-scan.
Length = 4096   # Number of pixels per A-scan.
Depth = 4
lambda0 = 842.5   # Central wavelength in nm.
sampFreq = 40000   # Samapling frequency of system in Hz (i.e., number of A-scans per second).
surf = 917


' Importing TDMS B-scan Data '

file_path = r"C:\80k Data\TDMS Data\Mscans\Piezodisp05.08.24\102.6Vpp_1.5kHz_37pulses_1024Ascans4Bscans.tdms"

# Bscan_data = read_1d_tdms(file_path)
# Bscan = Bscan_data.reshape((length, width), order='F')


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
        corrected_spec = Array[i,:,n]*np.exp(-1j*(np.squeeze(Disp_data)))   # Darven's disp comp.
        #corrected_spec = Array[i,:,n]*np.exp(1j*(disp_comp_array)) # [i,:,n]
        Array[i,:,n] = corrected_spec


' ~ Background Subtraction ~ '
ArraySample = Array[int(np.ceil(np.shape(Array)[0]/2)),:,:]
Background = np.mean(ArraySample,axis = 1)
Background = np.repeat(Background[:,None],(np.shape(ArraySample)[1]),axis = 1)
Array = Array - Background[None,:,:]


' ~ Fourier Transform Loop ~ '

N_size = 2**14
ArrayN_size = np.zeros([Depth,Length*4,Width]) # Length*4
ArrayN_size = ArrayN_size.astype(complex)
for i in range(0,Depth):
    Bscan = Array[i,:,:]
    FFT_Bscan = np.fft.fft(np.transpose(np.transpose(Bscan)*np.hanning(Length)), n=N_size, axis=0) #, n=N_size, 
    ArrayN_size[i,:,:] = FFT_Bscan
    
Array = ArrayN_size[:,0:int(N_size/2),:] # 0:int(N_size/2)
#Array = np.roll(Array,-780,axis = 2)  # Fix vertical displacement [ | ]


' Data Processing '

# N_size = 2**14  # Zero-padding scale
surface = int(surf/Length*N_size)   # Row in B-scan that corresponds to the surface of the sample.
# FFT_Bscan = np.fft.fft(np.transpose(np.transpose(Array[0,:,:])*np.hanning(Length)), n=N_size, axis=0)
# #FFT_Bscan = np.fft.fft(Bscan, n=N_size, axis=0)
# FFT_Bscan = FFT_Bscan[0:int(N_size/2), :]
# #FFT_Bscan = FFT_Bscan[0:int(length/2), :]

plt.figure(1, figsize=(8, 6))
plt.imshow(np.log(np.abs(Array[0,:,:])), cmap='gray', aspect='auto', vmin=6, vmax=12)   # Fourier transformed 2D image.
plt.title('B-scan')
plt.xlabel('A-scan Number (px)')
plt.ylabel('Depth (px)')
plt.axis('on')



' Phase Analysis '

bscanPhase1 = np.angle(Array)   # Extract phase information.

bscanPhaseUW1 = np.unwrap(bscanPhase1)   # Unwrapped phase.

bscanPhaseUW1 = np.mean(bscanPhaseUW1, axis=0)

meanSub = np.mean(np.abs(bscanPhaseUW1), axis=1)
frequencySignal = np.zeros([int(N_size/2), Width])
for i in range(0,int(N_size/2)):
        frequencySignal[i,:] =  (np.abs(bscanPhaseUW1[i,:]) - meanSub[i])

bscanPhaseUW = np.zeros([int(N_size/2), Width])
for i in range(0,Width):
    bscanPhaseUW[:,i] = frequencySignal[:,i]*np.hanning(int(N_size/2))


disp = (bscanPhaseUW[917,:])*lambda0/(4*np.pi) # Rescaling the y-axis amplitude to get surface displacement.

plt.figure()
plt.plot(disp)

FFT_UW_signal = np.fft.fft(bscanPhaseUW, axis=1)   # Fourier transform of unwrapped phase.
FFT_UW_signal = (FFT_UW_signal/Width)*lambda0/(4*np.pi) # Rescaling the y-axis amplitude to get surface displacement.


' Plotting '

# plt.figure()3
# plt.imshow(bscanPhase1)
# plt.title('B-scan Phase')
# plt.xlabel('Depth (px)')
# plt.ylabel('Phase (rad)')

# plt.figure()
# plt.imshow(bscanPhaseUW)
# plt.title('B-scan Unwrapped Phase')
# plt.xlabel('Depth (px)')
# plt.ylabel('Phase (rad)')

# plt.figure()
# plt.plot(np.abs(FFT_UW_signal))
# plt.title('B-scan FT Unwrapped Phase')
# plt.xlabel('Depth (px)')
# plt.ylabel('Amplitude (Linear)')

plt.figure(2)
finalSignal = np.abs(FFT_UW_signal[surf,0:int(Width/2)])  # Multiplying by 4 to compensate for 4x data points from zero-padding.
rescaledX = np.arange(0, int(Width/2))*(sampFreq/Width)
plt.plot(rescaledX/1e3, finalSignal)
plt.title('Piezo Surface Signal')
plt.ylabel('Displacement (nm)')
plt.xlabel('Frequency (kHz)')
plt.xlim([0.5, 10])
#plt.ylim([0, 10])
plt.grid('True')

# standDev = np.std(np.abs(FFT_UW_signal[surface,0:int(width/2)])[250:500])
# print('Phase Stability:', (standDev/(2*np.pi)*842.5e-9)*1000, 'mrad')