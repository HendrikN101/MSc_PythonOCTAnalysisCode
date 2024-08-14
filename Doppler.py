# -*- coding: utf-8 -*-
"""
Created on Wed Sep 20 13:55:37 2023

@author: Hendrik
"""

import numpy as np
import matplotlib.pyplot as plt
#import tifffile as tf   # External download
#from scipy.fft import fft, ifft
from scipy import signal
#from nptdms import TdmsFile   # External Download
#import utils as ut
#from PIL import Image
#from skimage import data, img_as_float, color, exposure
from skimage.restoration import unwrap_phase
from scipy.signal import convolve
#import matplotlib as mpl
from matplotlib.colors import ListedColormap
import matplotlib.colors as colors
from scipy.signal import butter,filtfilt,medfilt
from nptdms import TdmsFile
import scipy.io

plt.close('all')

' Importing Linearisation Files '

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

Disp_mat = scipy.io.loadmat(r'C:\80k Data\Lin_Disp_Files\From Python\Disp.mat')
Disp_data = np.array(Disp_mat['Disp'])


' Function to import TDMS Files '

def read_1d_tdms(file_path):
    with TdmsFile.open(file_path) as tdms_file:
        group = tdms_file.groups()[0]
        channel = group.channels()[0]
        data = channel[:]
        return data


' Importing TDMS M-scan Data '

file_path = r"C:\80k Data\Projects\DopplerFlow\Trial6\5Degrees\80uLmin_1.5mm_1024x10.tdms"

# Set C-scan Size
Width = 1024
Depth = 10
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
        corrected_spec = Array[i,:,n]*np.exp(-1j*(np.squeeze(Disp_data)))
        Array[i,:,n] = corrected_spec


' ~ Background Subtraction ~ '

ArraySample = Array[int(np.ceil(np.shape(Array)[0]/2)),:,:]
Background = np.mean(ArraySample,axis = 1)
Background = np.repeat(Background[:,None],(np.shape(ArraySample)[1]),axis = 1)
Array = Array - Background[None,:,:]


' ~ Fourier Transform Loop ~ '

# ZeroPad = 2**14
# ArrayZeroPad = np.zeros([Depth,Length*4,Width]) # Length*4
for i in range(0,Depth):
    Bscan = Array[i,:,:]
    FFT_Bscan = np.fft.fft(Bscan, axis=0) #, n=ZeroPad, 
    Array[i,:,:] = FFT_Bscan
    
Array = Array[:,0:int(Length/2),:] # 0:int(ZeroPad/2)

#Array = np.roll(Array,-425,axis = 2)  # Fix vertical displacement [ | ], 40 uL/min Double Flow

'5 Degrees Trial 7'
#Array = np.roll(Array,-579,axis = 2)  # Fix vertical displacement [ | ], 100 uL/min Double Flow
#Array = np.roll(Array,-258,axis = 2)  # Fix vertical displacement [ | ], 200 uL/min Double Flow
#Array = np.roll(Array,-258,axis = 2)  # Fix vertical displacement [ | ], 300 uL/min Double Flow
#Array = np.roll(Array,-715,axis = 2)  # Fix vertical displacement [ | ], 400 uL/min Double Flow
#Array = np.roll(Array,-701,axis = 2)  # Fix vertical displacement [ | ], 450 uL/min Double Flow

'5 Degrees Trial 6'
#Array = np.roll(Array,-734,axis = 2)  # Fix vertical displacement [ | ], 40 uL/min Double Flow
#Array = np.roll(Array,-643,axis = 2)  # Fix vertical displacement [ | ], 60 uL/min Double Flow
Array = np.roll(Array,-685,axis = 2)  # Fix vertical displacement [ | ], 80 uL/min Double Flow
#Array = np.roll(Array,-733,axis = 2)  # Fix vertical displacement [ | ], 100 uL/min Double Flow
#Array = np.roll(Array,-614,axis = 2)  # Fix vertical displacement [ | ], 120 uL/min Double Flow
#Array = np.roll(Array,-701,axis = 2)  # Fix vertical displacement [ | ], 140 uL/min Double Flow
#Array = np.roll(Array,-670,axis = 2)  # Fix vertical displacement [ | ], 200 uL/min Double Flow
#Array = np.roll(Array,-581,axis = 2)  # Fix vertical displacement [ | ], 250 uL/min Double Flow
#Array = np.roll(Array,-533,axis = 2)  # Fix vertical displacement [ | ], 275 uL/min Double Flow
#Array = np.roll(Array,-537,axis = 2)  # Fix vertical displacement [ | ], 300 uL/min Double Flow
#Array = np.roll(Array,-537,axis = 2)  # Fix vertical displacement [ | ], 300 uL/min Double Flow

'10 Degrees Trial 6'
#Array = np.roll(Array,-688,axis = 2)  # Fix vertical displacement [ | ], 40 uL/min Double Flow
#Array = np.roll(Array,-385,axis = 2)  # Fix vertical displacement [ | ], 60 uL/min Double Flow
#Array = np.roll(Array,-555,axis = 2)  # Fix vertical displacement [ | ], 80 uL/min Double Flow
#Array = np.roll(Array,-562,axis = 2)  # Fix vertical displacement [ | ], 100 uL/min Double Flow
#Array = np.roll(Array,-370,axis = 2)  # Fix vertical displacement [ | ], 120 uL/min Double Flow
#Array = np.roll(Array,-623,axis = 2)  # Fix vertical displacement [ | ], 140 uL/min Double Flow

'5 Degrees Trial 5'
#Array = np.roll(Array,-627,axis = 2)  # Fix vertical displacement [ | ], 40 uL/min Double Flow
#Array = np.roll(Array,-837,axis = 2)  # Fix vertical displacement [ | ], 60 uL/min Double Flow
#Array = np.roll(Array,-742,axis = 2)  # Fix vertical displacement [ | ], 80 uL/min Double Flow
#Array = np.roll(Array,-744,axis = 2)  # Fix vertical displacement [ | ], 100 uL/min Double Flow
#Array = np.roll(Array,-792,axis = 2)  # Fix vertical displacement [ | ], 200 uL/min Double Flow


'10 Degrees Trial 5'
#Array = np.roll(Array,-474,axis = 2)  # Fix vertical displacement [ | ], 20 uL/min Double Flow
#Array = np.roll(Array,-615,axis = 2)  # Fix vertical displacement [ | ], 40 uL/min Double Flow
#Array = np.roll(Array,-776,axis = 2)  # Fix vertical displacement [ | ], 60 uL/min Double Flow
#Array = np.roll(Array,-700,axis = 2)  # Fix vertical displacement [ | ], 80 uL/min Double Flow
#Array = np.roll(Array,-735,axis = 2)  # Fix vertical displacement [ | ], 100 uL/min Double Flow
#Array = np.roll(Array,-524,axis = 2)  # Fix vertical displacement [ | ], 150 uL/min Double Flow




Phi = 0.5467 #np.mean(phase_plot)   # Phase difference [Radians].
Lambda_0 = 842.5e-9   # Central wavelength [Metres].
Tau = 1/40000   # Time period - inverse A-scan rate [Seconds].
Ref_Index = 1.33   # For milk.
Theta = 5 * np.pi/180  # Angle between horizontal and flow direction [Radians].

Vel = (Phi * Lambda_0) / (4*np.pi*Tau*Ref_Index*np.cos((np.pi/2)-Theta))

print('Average Flow Velocity:', np.round(Vel, 5), 'm/s', ' ', '[', np.round(Vel*1e3, 2), 'mm/s ]')



' ~ Data Processing ~ '

depth = np.linspace(0, int(np.size(Array, axis=0)-1), int(np.size(Array, axis=0)))

Ascan_1 = Array[0,:,0]
Ascan_2 = Array[1,:,0]
Ascan_3 = Array[2,:,0]
Ascan_4 = Array[3,:,0]

Bscan_1 = Array[0,:,:]
Bscan_2 = Array[1,:,:]
Bscan_3 = Array[2,:,:]
Bscan_4 = Array[3,:,:]


' ~ Doppler Imaging ~ '

# Improvement: Complex conjugate before np.angle, remove np.diff after.

bscanPhase1 = np.angle(Bscan_1)   # Extract phase information.

bscanPhaseUW1 = np.unwrap(bscanPhase1)   # Unwrapped phase.

bscanPhaseDiff = np.diff(bscanPhaseUW1, axis=0)   # Phase difference of adjacent A-scans in the same B-scan.

image_height = bscanPhaseDiff[:,0].size
image_width = bscanPhaseDiff[0,:].size
target_columns = 100
columns_per_group = image_width // target_columns
Ascan_averaged = bscanPhaseDiff[:, :target_columns*columns_per_group].reshape(image_height, target_columns, columns_per_group)
Ascan_averaged = np.mean(Ascan_averaged, axis=2)   # A-scan averaged phase difference.

Bscan_List = [Bscan_1, Bscan_2, Bscan_3, Bscan_4]
Bscan_Processed_List = []
for i in range(0,Depth):
    # A = Array[i,:,:]
    # B = np.roll(A, -1, axis = 1)
    # N = A * np.conj(B)
    bscanPhase = np.angle(Array[i,:,:])   # Extract phase information.
    bscanPhaseUW = np.unwrap(bscanPhase)   # Unwrapped phase.
    bscanPhaseUW = bscanPhaseUW[:,12:1013]   # Reduces to 1001 so after diff will be 1000, which is then scaled down to 100 later.
    bscanPhaseDiff = np.diff(bscanPhaseUW, axis=1)   # Phase difference of adjacent A-scans in the same B-scan.
    image_height = bscanPhaseDiff[:,0].size
    image_width = bscanPhaseDiff[0,:].size
    columns_per_group = image_width // target_columns
    Ascan_averaged = bscanPhaseDiff[:, :target_columns*columns_per_group].reshape(image_height, target_columns, columns_per_group)
    Ascan_averaged = np.mean(Ascan_averaged, axis=2)   # A-scan averaged phase difference.
    Bscan_Processed_List.append(Ascan_averaged)
    
Bscan_Processed_Array = np.array(Bscan_Processed_List)
Bscan_Averaged_Array = np.mean(Bscan_Processed_Array, axis=0)




' ~ Zero Flow ~ '

' Importing TDMS M-scan Data '

file_path_0 = r"C:\80k Data\Projects\DopplerFlow\Trial6\5Degrees\ZeroFlow_1.5mm_1024x10.tdms"

Array_0 = np.zeros([Depth,Length,Width])   # [Depth, Length, Width]
Array_0 = Array_0.astype(complex)

with TdmsFile.read(file_path_0) as tdms_file:   # TdmsFile.read writes straight to memory!
    for i in range(0,Depth):
        for n in range(0,Width):
              group = tdms_file["Buffer " + str('{0:03}'.format(i))]
              all_group_channels = group.channels()
              Array_0[i,:,n] = all_group_channels[n][:]


' ~ Linearisation Loop ~ '

for n in range(0,Width):
    for i in range(0,Depth):
        interpolate = np.interp(lin_array, non_lin_array, Array_0[i,:,n])
        Array_0[i,:,n] = interpolate


' ~ Dispersion Compensation Loop ~ '

for n in range(0,Width):
    for i in range(0,Depth):
        corrected_spec = Array_0[i,:,n]*np.exp(-1j*(np.squeeze(Disp_data)))
        Array_0[i,:,n] = corrected_spec


' ~ Background Subtraction ~ '

ArraySample = Array_0[int(np.ceil(np.shape(Array_0)[0]/2)),:,:]
Background = np.mean(ArraySample,axis = 1)
Background = np.repeat(Background[:,None],(np.shape(ArraySample)[1]),axis = 1)
Array_0 = Array_0 - Background[None,:,:]



' ~ Fourier Transform Loop ~ '

# ZeroPad = 2**14
# ArrayZeroPad = np.zeros([Depth,Length*4,Width]) # Length*4
for i in range(0,Depth):
    Bscan = Array_0[i,:,:]
    FFT_Bscan = np.fft.fft(Bscan, axis=0) #, n=ZeroPad, 
    Array_0[i,:,:] = FFT_Bscan
    
Array_0 = Array_0[:,0:int(Length/2),:] # 0:int(ZeroPad/2)

#Array_0 = np.roll(Array_0,-584,axis = 2)  # Fix vertical displacement [ | ]

'Trial 7'
#Array_0 = np.roll(Array_0,-684,axis = 2)  # Fix vertical displacement [ | ]

'5 Degrees Trial 6'
Array_0 = np.roll(Array_0,-774,axis = 2)  # Fix vertical displacement [ | ]

'10 Degrees Trial 6'
#Array_0 = np.roll(Array_0,-805,axis = 2)  # Fix vertical displacement [ | ]

'5 Degrees Trial 5'
#Array_0 = np.roll(Array_0,-536,axis = 2)  # Fix vertical displacement [ | ]

'10 Degrees Trial 5'
#Array_0 = np.roll(Array_0,-586,axis = 2)  # Fix vertical displacement [ | ]


Ascan_1_0 = Array_0[0,:,0]
Ascan_2_0 = Array_0[1,:,0]
Ascan_3_0 = Array_0[2,:,0]
Ascan_4_0 = Array_0[3,:,0]

Bscan_1_0 = Array_0[0,:,:]
Bscan_2_0 = Array_0[1,:,:]
Bscan_3_0 = Array_0[2,:,:]
Bscan_4_0 = Array_0[3,:,:]

# plt.imshow(np.log(np.abs(Bscan_1_0)), cmap='gray', vmin=5, vmax=10)


Bscan_List_0 = [Bscan_1_0, Bscan_2_0, Bscan_3_0, Bscan_4_0]
Bscan_Processed_List_0 = []
for i in range(0,Depth):
    # A = Array_0[i,:,:]
    # B = np.roll(A, 1, axis = 1)
    # N = A * np.conj(B)
    bscanPhase = np.angle(Array_0[i,:,:])   # Extract phase information.
    bscanPhaseUW = np.unwrap(bscanPhase)   # Unwrapped phase.
    bscanPhaseUW = bscanPhaseUW[:,12:1013]   # Reduces to 1001 so after diff will be 1000, which is then scaled down to 100 later.
    bscanPhaseDiff = np.diff(bscanPhaseUW, axis=1)   # Phase difference of adjacent A-scans in the same B-scan.
    image_height = bscanPhaseDiff[:,0].size
    image_width = bscanPhaseDiff[0,:].size
    columns_per_group = image_width // target_columns
    Ascan_averaged = bscanPhaseDiff[:, :target_columns*columns_per_group].reshape(image_height, target_columns, columns_per_group)
    Ascan_averaged = np.mean(Ascan_averaged, axis=2)   # A-scan averaged phase difference.
    Bscan_Processed_List_0.append(Ascan_averaged)
    
Bscan_Processed_Array_0 = np.array(Bscan_Processed_List_0)
Bscan_Averaged_Array_0 = np.mean(Bscan_Processed_Array_0, axis=0)

Final_Phase_Image = Bscan_Averaged_Array - Bscan_Averaged_Array_0

Final2_Phase_Image = unwrap_phase(Final_Phase_Image)

Filtered_Final = medfilt(Final_Phase_Image,5)   # Final image filtered.



' ~ Flow Velocity [UPDATE VALUES] ~ '

phase_prof_x = np.mean(Filtered_Final[:,26:31], axis=1)   # Phase profile in the vertical direction through flow.
phase_prof_savgol_x = signal.savgol_filter(phase_prof_x, 150, 3)
#print('Maximum Phase Difference (Vertical):', round(np.max(phase_prof_x), 5), 'Radians')
print('Minimum Phase Difference (Vertical, Left Flow):', round(np.max(phase_prof_x), 5), 'Radians')


phase_prof_x2 = np.mean(Filtered_Final[:,72:77], axis=1)   # Phase profile in the vertical direction through flow. USED FOR FLOW PAIR.
phase_prof_savgol_x2 = signal.savgol_filter(phase_prof_x2, 150, 3)
print('Maximum Phase Difference (Vertical, Right Flow):', round(np.min(phase_prof_x2), 5), 'Radians')
#print('Minimum Phase Difference (Vertical):', round(np.min(phase_prof_x2), 5), 'Radians')


phase_prof_y = np.mean(Filtered_Final[445:495,:], axis=0)  # Phase profile in the horizontal direction through flow. 454:494
phase_prof_savgol_y = signal.savgol_filter(phase_prof_y, 15, 3)
print('Maximum Phase Difference (Horizontal):', round(np.max(phase_prof_y), 5), 'Radians')
print('Minimum Phase Difference (Horizontal):', round(np.min(phase_prof_y), 5), 'Radians')



phase_plot = Final_Phase_Image[250:280,45:50]   # Finding the mean phase difference ([:, x] = [Depth Range, A-scan Num])







' ~ Plotting ~ '

# # Plotting the data for an individual A-scan
# plt.figure(figsize=(10, 8))
# plt.subplot(4, 1, 1)
# plt.plot(depth, np.real(Ascan_1), label='Real Part')
# plt.plot(depth, np.imag(Ascan_1), label='Imaginary Part', color='darkorange')
# plt.xlabel('Depth (Pixels)')
# plt.ylabel('Amplitude')
# plt.title('Real & Imag. A-scan Data')
# plt.legend()

# plt.subplot(4, 1, 2)
# plt.plot(depth, np.abs(Ascan_1))
# plt.xlabel('Depth (Pixels)')
# plt.ylabel('Amplitude')
# plt.title('A-scan')

# plt.subplot(4, 1, 3)
# plt.plot(depth, phase, label='Phase', color='darkorange')
# plt.xlabel('Depth (Pixels)')
# plt.ylabel('Phase (Radians)')
# plt.title('Phase')

# plt.subplot(4, 1, 4)
# plt.plot(depth, unwrapped_phase, label='Unwrapped Phase', color='red')
# plt.xlabel('Depth (Pixels)')
# plt.ylabel('Phase (Radians)')
# plt.title('Unwrapped Phase')

# plt.tight_layout()
# plt.show()


# Plotting Structural OCT Image (B-scan)
plt.figure()
plt.imshow(np.log(np.abs(Bscan_1[0:1000,:])), cmap='gray', aspect='equal', vmin=5, vmax=10)
plt.title('OCT Image (Single B-scan)')
plt.xlabel('A-scan Number')
plt.ylabel('Depth (Pixels)')
#plt.axis('off')
plt.show()

# Plotting Structural OCT Image Zero Flow (B-scan)
plt.figure()
plt.imshow(np.log(np.abs(Array_0[1,0:1000,:])), cmap='gray', aspect='equal', vmin=5, vmax=10)
plt.title('OCT Image (Single B-scan - Zero Flow)')
plt.xlabel('A-scan Number')
plt.ylabel('Depth (Pixels)')
#plt.axis('off')
plt.show()

# Plotting Phase Image
plt.figure()
plt.imshow(bscanPhase1[0:1000,:], cmap='seismic', aspect='equal', norm=colors.CenteredNorm())
plt.title('Phase Image (Single B-scan)')
plt.xlabel('A-scan Number')
plt.ylabel('Depth (Pixels)')
plt.colorbar()
plt.show()

# Plotting Unwrapped Phase Image
plt.figure()
plt.imshow(bscanPhaseUW1[0:1000,:], cmap='seismic', aspect='equal', norm=colors.CenteredNorm())
plt.title('Unwrapped Phase Image (Single B-scan)')
plt.xlabel('A-scan Number')
plt.ylabel('Depth (Pixels)')
plt.colorbar()
plt.show()

# Plotting Phase Difference Image (Adjacent A-scans in the same B-scan)
plt.figure()
plt.imshow(bscanPhaseDiff[0:1000,:], cmap='seismic', aspect='equal', norm=colors.CenteredNorm())
plt.title('Phase Difference Image (Single B-scan)')
plt.xlabel('A-scan Number')
plt.ylabel('Depth (Pixels)')
plt.colorbar()
plt.show()

# Plotting Phase Difference Image w/ Averaged A-scans
plt.figure()
plt.imshow(Ascan_averaged[0:1000,:], cmap='seismic', aspect=1/10, norm=colors.CenteredNorm())
plt.title('Phase Difference Image (Single B-scan, A-scan Averaged)')
plt.xlabel('A-scan Number')
plt.ylabel('Depth (Pixels)')
plt.colorbar()
plt.show()

# Plotting Phase Difference Image w/ Averaged A-scans and B-scans
plt.figure()
plt.imshow(Bscan_Averaged_Array[0:1000,:], cmap='seismic', aspect=1/10, norm=colors.CenteredNorm())
plt.title('Phase Difference Image (A-scan & B-scan Averaged)')
plt.xlabel('A-scan Number')
plt.ylabel('Depth (Pixels)')
plt.colorbar()
plt.show()

# Plotting Phase Difference Image w/ Averaged A-scans and B-scans (ZERO FLOW)
plt.figure()
plt.imshow(Bscan_Averaged_Array_0[0:1000,:], cmap='seismic', aspect=1/10, norm=colors.CenteredNorm())
plt.title('Zero Flow Phase Difference Image (A-scan & B-scan Averaged)')
plt.xlabel('A-scan Number')
plt.ylabel('Depth (Pixels)')
plt.colorbar()
plt.show()

# Plotting Phase Difference Image w/ Averaged A-scans and B-scans (ZERO FLOW SUBTRACTED)
plt.figure()
plt.imshow(Final_Phase_Image[0:1000,:], cmap='seismic', aspect=1/10, norm=colors.CenteredNorm())
plt.title('Final Phase Difference Image (A-scan & B-scan Averaged)')
plt.xlabel('A-scan Number')
plt.ylabel('Depth (Pixels)')
plt.colorbar()
plt.show()

# Plotting Phase Difference Image w/ Averaged A-scans and B-scans (ZERO FLOW SUBTRACTED, FILTERED)
plt.figure()
plt.imshow(Filtered_Final[0:1000,:], cmap='seismic', aspect=1/10, norm=colors.CenteredNorm())
plt.title('Final Phase Difference Image (A-scan & B-scan Averaged - Filtered)')
plt.xlabel('A-scan Number')
plt.ylabel('Depth (Pixels)')
plt.colorbar()
plt.show()
 


# Velocity Profile in the vertical direction through flow.
plt.figure()
plt.plot(phase_prof_x[0:1000], 'r')
plt.plot(phase_prof_savgol_x[0:1000], 'b');plt.title('Vertical Phase Profile through Left Flow');plt.xlabel('Depth (Pixels)');plt.ylabel('Phase Difference (rad)')
plt.legend(['Raw Phase Difference','Filtered Phase Difference'])
plt.grid('True')
plt.show()

# Velocity Profile in the vertical direction through flow.
plt.figure()
plt.plot(phase_prof_x2[0:1000], 'r')
plt.plot(phase_prof_savgol_x2[0:1000], 'b');plt.title('Vertical Phase Profile through Right Flow');plt.xlabel('Depth (Pixels)');plt.ylabel('Phase Difference (rad)')
plt.legend(['Raw Phase Difference','Filtered Phase Difference'])
plt.grid('True')
plt.show()

# Velocity Profile in the horizontal direction through flow
plt.figure(figsize=(7,4))
plt.plot(phase_prof_y, 'r')
plt.plot(phase_prof_savgol_y, 'b');plt.title('Horizontal Phase Profile through Flow');plt.xlabel('Width (Pixels)');plt.ylabel('Phase Difference (rad)')
plt.legend(['Raw Phase Difference','Filtered Phase Difference'])
plt.grid('True')
plt.show()


# # Plotting Individual A-scan to see Phase Difference
# plt.plot(phase_plot)
# plt.title('[Individual A-scan] Phase Difference Plot (A-scan & B-scan Averaged)')
# plt.xlabel('Depth (Pixels)')
# plt.ylabel('Phase Difference (Radians)')
# plt.show()

# # Plotting Phase Difference Image w/ Averaged A-scans and B-scans (ZERO FLOW SUBTRACTED)
# plt.imshow(Final2_Phase_Image, cmap='seismic', aspect='auto', norm=colors.CenteredNorm())
# plt.title('Final Phase Difference Image (A-scan & B-scan Averaged)')
# plt.xlabel('A-scan Number')
# plt.ylabel('Depth (Pixels)')
# plt.colorbar()
# plt.show()

# # Better Doppler Plotting
# plt.figure();plt.imshow(g,aspect="auto",cmap='seismic',vmin=-2,vmax=2);plt.colorbar(label='Unwrapped Phase Difference (rad)')
# plt.xlabel('A-scan Number');plt.ylabel('Depth (px)');plt.title('Unwrapped Phase Difference of Milk Flowing in Opposite Directions');plt.show()


# def butter_lowpass_filter(data, cutoff, fs, order):
#     nyq = 2 * 100000
#     normal_cutoff = cutoff / nyq
#     b, a = butter(order, normal_cutoff, btype='low', analog=False)
#     y = filtfilt(b, a, data)
#     return y



