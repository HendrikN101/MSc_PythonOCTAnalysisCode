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


# ' ~ Importing Linearisation Files ~ '

# lin_txt = open(r'C:\Users\User\OneDrive\Documents\2023-2024\Uni\Masters\80k System\Spectrum\lin_80k.txt', 'r')
# lin_array = np.array(lin_txt.read().split())
# lin_txt.close()
# non_lin_txt = open(r'C:\Users\User\OneDrive\Documents\2023-2024\Uni\Masters\80k System\Spectrum\non_lin_80k.txt', 'r')
# non_lin_array = np.array(non_lin_txt.read().split())
# non_lin_txt.close()

# def convert_strings_to_floats(input_array):
#     output_array = []
#     for element in input_array:
#         converted_float = float(element)
#         output_array.append(converted_float)
#     return output_array

# lin_array = convert_strings_to_floats(lin_array)
# non_lin_array = convert_strings_to_floats(non_lin_array)


' ~ Function to import TDMS Files ~ '

def read_1d_tdms(file_path):
    with TdmsFile.open(file_path) as tdms_file:
        group = tdms_file.groups()[0]
        channel = group.channels()[0]
        data = channel[:]
        return data



' ~ Importing TDMS Spectral Data ~ '

file_path = r"C:\80k Data\TDMS Data\Spectrums\Spec4Disp_OD2_LinApplied_NEW.tdms"

Ascan_data = read_1d_tdms(file_path)


' ~ Background Subtraction ~ '

Background = np.mean(Ascan_data, axis = 0)
Ascan_data = Ascan_data - Background


plt.figure(1)
plt.plot(Ascan_data)

FFTdata = np.fft.fftshift(np.abs(np.fft.fft(Ascan_data, n=2**14)))
x_axis = np.linspace(0,4095,2**14)

plt.figure(2)
plt.plot(x_axis, FFTdata)