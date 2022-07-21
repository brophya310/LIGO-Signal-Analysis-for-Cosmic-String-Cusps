#CBC vs CUSP injection testing
#CBC Injection Recovery Program

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
import readligo as rl
import template 
import FindCuspTemplate 

#m1 = input(print("What is the first mass?"))
#m2 =input(print("What is the second mass?"))
#Deff = input(print("What is the effective distance?"))

#Read in an already downloaded file, set dt and sampling frequency
#Read Ligo takes the arrays for strain and time from the hdf5 file
strain, strain_times, dq = rl.loaddata('H-H1_LOSC_4_V1-816578560-4096.hdf5')
if len(strain) < 2:
    print("Failed to read file.")
    exit()

dt = strain_times[1] - strain_times[0]
fs = 1.0 / dt


#Creating an Injection Segment from GPS Time
#given the gps time of the injection
gpstime_event = 816579666   
start = gpstime_event - 32
stop = gpstime_event + 32


#seg = rl.getsegs(start, stop, 'H1')[0]
#for (begin, end) in segs:
seg_strains, meta, dq = rl.getstrain(start, stop, 'H1')
seg_duration = len(seg_strains)/fs
seg_times = np.arange(0, len(seg_strains))*dt + start
print("Last segment time: ", seg_times[-1], " GPS end: " , stop)

print ("The injection segment is {0} s long".format(seg_duration))

#create template, frequency domain inspiral template
temp, temp_freq = template.createTemplate(fs, seg_duration, 10, 10)
ctemp, ctemp_freq = FindCuspTemplate.create_a_template(fs, seg_duration, 25, 220)
#set template to zero at low frequencies
temp[temp_freq < 25] = 0
ctemp[ctemp_freq < 25] = 0

#Show template value vs frequency
plt.figure()
plt.plot(temp_freq, abs(temp))
plt.plot(ctemp_freq, abs(ctemp))
plt.axis([10, 1000, 1e-22, 1e-19])
plt.xlabel("Frequency (Hz)")
plt.ylabel("Template value (Strain/Hz")
plt.grid()
#plt.show()

#Show IFFT of template and plot as time series
"""
t_temp = abs(temp)
t = np.arange(len(t_temp))
invtemp = np.fft.ifft(temp)
plt.figure()
plt.plot(t, invtemp.real)
plt.axis([0,1, 1e-22, 1e-19])
plt.xlabel("Time (s)")
plt.ylabel("Template value (Strain)")
plt.grid
plt.show()
"""
         
#fft the data and add window
window = np.blackman(seg_strains.size)
data_fft = np.fft.rfft(seg_strains*window)

#get noise segment for pxx
noise_strain, meta1, dq1 = rl.getstrain(start-512, gpstime_event-1, 'H1')

#setting NFFT equal to the size of the data segment
Pxx, psd_freq = mlab.psd(noise_strain, Fs=fs, NFFT=len(seg_strains))

#computing matched filter output, multiply data by the template and weight by the PSD
integrand = data_fft*np.ma.conjugate(temp)/Pxx
cusp_integrand = data_fft*np.ma.conjugate(ctemp)/Pxx
#
num_zeros = len(seg_strains) - len(data_fft)
#print("The number of zeros to pad:", num_zeros)
padded_int = np.append( integrand, np.zeros(num_zeros) )
cusp_padded_int = np.append( cusp_integrand, np.zeros(num_zeros) )
z = 4*np.fft.ifft(padded_int)
cz = 4*np.fft.ifft(cusp_padded_int)


###Compute the Normalization
kernal =  (np.abs(temp))**2 / Pxx
df = psd_freq[1] - psd_freq[0]
sig_sqr = 4*kernal.sum()*df
sigma = np.sqrt(sig_sqr)

###Compute the Normalization
cusp_kernal =  (np.abs(ctemp))**2 / Pxx
cusp_df = psd_freq[1] - psd_freq[0]
cusp_sig_sqr = 4*cusp_kernal.sum()*df
cusp_sigma = np.sqrt(cusp_sig_sqr)


expected_SNR = sigma / 10

####Compute the SNR
inv_win = (1.0 / window)
inv_win[:20*4096] = 0
inv_win[-20*4096:] = 0
rho = abs(z) / sigma * inv_win
crho = abs(cz) / cusp_sigma * inv_win


####Output the results
# -- Plot rho as a function of time
#print(seg_times.size)
#print(rho.size)
plt.figure()
plt.plot(seg_times[::8]-seg_times[0], rho[::8])
plt.plot(seg_times[::8]-seg_times[0], crho[::8])
plt.xlabel("Seconds since GPS {0:.0f}".format(seg_times[0]) )
plt.ylabel("SNR")
plt.show()

#Compute and display results
snr = rho.max()
snrx = np.average(rho)
snrr = snr / snrx
found_time = seg_times[ np.where(rho == snr) ]

cusp_snr = crho.max()
cusp_snrx = np.average(crho)
cusp_snrr = cusp_snr / cusp_snrx
cusp_found_time = seg_times[ np.where(rho == cusp_snr) ]


print("The expected SNR is ", expected_SNR)
print ("Recovered SNR is ", snr)
print ("The SNRR is " , snrr)
print ("Recovered time GPS {0:.1f}".format( found_time[0] ))

print("The expected Cusp SNR is ", expected_SNR)
print ("Recovered  Cusp SNR is ", snr)
print ("The Cusp SNRR is " , snrr)
print ("Recovered time GPS for Cusp{0:.1f}".format( found_time[0] ))
