# -*- coding: utf-8 -*-
"""
Created on Mon Apr  8 11:13:54 2024

@author: hnie464
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
file_path = r"C:\80k Data\Projects\Agarose\Agarose\1.5%Agarose\Agarose1.5%_1kHz_to_20kHz_piezo\06.08.24\BscanDispComp.tdms"   # Spectra of mirror surface with linearisation applied.
data = read_1d_tdms(file_path)
print(data)


# plt.plot(data)
# plt.title('850 nm Spectrum')
# plt.ylabel('Amplitude (12-bit)')
# plt.xlabel('Line Camera Width (Pixels)')
# plt.grid('True')

data = data[0:4096]

window_data = np.hanning(4096)*data
fourier_data = np.fft.fft(window_data)   # Fourier transform of the data
plt.figure()
plt.plot(np.abs(fourier_data[10:300]))
fourier_data[0:104] = 0   # Setting everything but the signal to zero
fourier_data[135:] = 0
inv_fourier_data = np.fft.ifft(fourier_data)   # Inverse Fourier transform of the modified data


phase = np.angle(inv_fourier_data)   # Extracting the phase of the modified data
non_lin = np.unwrap(phase)   # Unwrapping the phase to get the non-linear data [Extract to .txt]
plt.figure()
plt.plot(non_lin, label='Non-Linear')


m = (np.max(non_lin)-np.min(non_lin))/4095
axis = np.linspace(0,4095,4096)
linear = m * axis + np.min(non_lin)   # Forming the linear fit [Extract to .txt]
plt.plot(linear, label='Linear')
plt.legend()
plt.grid('True')

# Polyfit non-linear curve
Xindex = np.linspace(0, 4095, 4096)
p = np.poly1d(np.polyfit(Xindex[225:3690], non_lin[225:3690], 7))  #[225:3690]
plt.plot(Xindex, p(Xindex))
# non_lin = p(Xindex)

plt.figure()
plt.plot(p(Xindex)-non_lin)

corrected_spec = data*np.exp(1j*(non_lin-linear))
disp_comp_data = 1*(non_lin-linear)   # Data required for compensation. [May have to multiply by -1 to get correct orientation]
ft_corrected_spec = np.fft.fft(corrected_spec)

plt.figure()
plt.plot(np.abs(np.fft.fft(data)), label='Before')
plt.plot(np.abs(ft_corrected_spec), label='After')
plt.legend()
plt.grid('True')


' ~ Writes a New Dispersion Compensation File [Uncomment] ~ '

disp_comp_txt = open(r'C:\80k Data\Lin_Disp_Files\From Python\disp_comp_80k.txt', 'w')   # Putting the dispersion compensation data into a .txt file

for i in range(0, len(disp_comp_data)):
    disp_comp_txt.write(str(disp_comp_data[i]) + '\t') #  + '\n'

disp_comp_txt.close()