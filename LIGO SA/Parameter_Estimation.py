#Parameter estimation loops
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
import readligo as rl
import template
import math
from matplotlib import cm           # color maps:
# Required for projection='3d'
from mpl_toolkits.mplot3d import Axes3D
# To customize the z axis
from matplotlib.ticker import LinearLocator, FormatStrFormatter

#LOAD Data, Setup event time and fs
strain, strain_times, dq = rl.loaddata("L-L1_GWOSC_4KHZ_R1-1128676853-4096.hdf5")
#if len(strain) < 2:
    #print("Failed to read file.")
    #exit()

dt = strain_times[1] - strain_times[0]
fs = 1.0 / dt


gpstime_event =  1128678900
start = gpstime_event - 10
stop = gpstime_event + 10

#universal necessary information for templates
seg_strains, meta, dq = rl.getstrain(start, stop, 'L1')
seg_duration = len(seg_strains)/fs
seg_times = np.arange(0, len(seg_strains))*dt + start

#put data in frequency domain and window data
window = np.blackman(seg_strains.size)
windowed_strain = seg_strains*window
data_fft = np.fft.rfft(windowed_strain)

#get noise segment for pxx
noise_strain, meta1, dq1 = rl.getstrain(start-512, gpstime_event-1, 'L1')
#noise_strain, meta1, dq1 = rl.getstrain( gpstime_event+2, 'H1')
#setting NFFT equal to the size of the data segment
Pxx, psd_freq = mlab.psd(noise_strain, Fs=fs, NFFT=len(seg_strains))


#GWOSC has reported 23 and 13 as the masses for this event
#this loop will calculate snr for the template for the range of possible m1 and m2

SNR_Max = 0
m1_optimal = 1
m2_optimal = 1
arr=[]

for m1 in range(1,31):
    for m2 in range(1,31):
        temp, temp_freq = template.createTemplate(fs, seg_duration, m1, m2)
        temp[temp_freq < 25] = 0
        #computing matched filter output, multiply data by the template and weight by the PSD
        integrand = data_fft*np.ma.conjugate(temp)/Pxx
        #
        num_zeros = len(seg_strains) - len(data_fft)
        #print("The number of zeros to pad:", num_zeros)
        padded_int = np.append( integrand, np.zeros(num_zeros) )
        z = 4*np.fft.ifft(padded_int)
        ###Compute the Normalization
        kernal =  (np.abs(temp))**2 / Pxx
        df = psd_freq[1] - psd_freq[0]
        sig_sqr = 4*kernal.sum()*df
        sigma = np.sqrt(sig_sqr)

        ####Compute the SNR
        inv_win = (1.0 / window)
        inv_win[:20*4096] = 0
        inv_win[-20*4096:] = 0
        rho = abs(z) / sigma * inv_win
        snr = rho.max()

        arr.append(snr)
        

        if(snr>SNR_Max):
            SNR_Max = snr
            m1_optimal = m1
            m2_optimal= m2

print("The template values for M1 and M2 that produce the highest SNR for the data are: M1 = ", m1_optimal , " and M2= ", m2_optimal)





#######CREATING COLOR MAP FOR SNR VS MASSES#######################





Xlow = 1
Xhigh = 30
Yhigh = 30
Ylow = 1
Zlow = 0
Zhigh = SNR_Max

X = np.arange(Xlow, Xhigh, 1)
Y = np.arange(Ylow,Yhigh, 1)
X, Y = np.meshgrid(X, Y)

Z = arr


# Make data.


fig, ax = plt.subplots(subplot_kw={"projection": "3d"})

# Plot the surface.
surf = ax.plot_surface(X, Y, np.array(Z), cmap=cm.coolwarm,
                       linewidth=0, antialiased=False)

# Customize the z axis.
ax.zaxis.set_major_locator(LinearLocator(2))
ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))

# Add a color bar which maps values to colors.
fig.colorbar(surf, shrink=0.5, aspect=5)

plt.show()







