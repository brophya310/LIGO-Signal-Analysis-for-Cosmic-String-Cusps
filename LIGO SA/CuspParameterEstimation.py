#Cusp Parameter Estimation
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
import readligo as rl
import math
from matplotlib import cm           # color maps:
# Required for projection='3d'
from mpl_toolkits.mplot3d import Axes3D
# To customize the z axis
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import FindCuspTemplate
#np.seterr(divide='ignore', invalid='ignore')



#Read in an already downloaded file, set dt and sampling frequency
#Read Ligo takes the arrays for strain and time from the hdf5 file
strain, strain_times, dq = rl.loaddata("H-H1_LOSC_4_V1-816332800-4096 (2).hdf5")
#if len(strain) < 2:
    #print("Failed to read file.")
    #exit()

dt = strain_times[1] - strain_times[0]
fs = 1.0 / dt


#Creating an Injection Segment from GPS Time
#given the gps time of the injection
gpstime_event =  816335770
start = gpstime_event - 32
stop = gpstime_event + 32

#for (begin, end) in segs:
seg_strains, meta, dq = rl.getstrain(start, stop, 'H1')
seg_duration = len(seg_strains)/fs
seg_times = np.arange(0, len(seg_strains))*dt + start

#put data in frequency domain and window data
window = np.blackman(seg_strains.size)
windowed_strain = seg_strains*window
data_fft = np.fft.rfft(windowed_strain)

#get noise segment for pxx
noise_strain, meta1, dq1 = rl.getstrain(start-512, gpstime_event-1, 'H1')
#noise_strain, meta1, dq1 = rl.getstrain( gpstime_event+2, 'H1')
#setting NFFT equal to the size of the data segment
Pxx, psd_freq = mlab.psd(noise_strain, Fs=fs, NFFT=len(seg_strains))




#GWOSC has reported 23 and 13 as the masses for this event
#this loop will calculate snr for the template for the range of possible m1 and m2

def getSNR(F1, F2):
        temp, temp_freq = FindCuspTemplate.create_a_template(fs, seg_duration, F1, F2)
        #computing matched filter output, multiply data by the template and weight by the PSD
        integrand = data_fft*np.ma.conjugate(temp)/Pxx
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
        return snr
    


Xlow = 50
Xhigh = 500
Ylow = 550
Yhigh = 1000

X = np.arange(Xlow, Xhigh,25)
Y = np.arange(Ylow,Yhigh,25)
X, Y = np.meshgrid(X, Y)

#Z = getSNR(X,Y)#Does not work becuase template can not take array

Z = np.zeros((X.shape))
for i in range(50,X.shape[0],25):
        Z[i-1][j-1] = getSNR(freq1, freq2)
        
print(Z)
"""
     
fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
# Plot the surface.
surf = ax.plot_surface(X, Y, Z, cmap=cm.inferno,
                       linewidth=0, antialiased=False)
# Customize the z axis.
#ax.set_zlim(0, Z.max())
ax.zaxis.set_major_locator(LinearLocator(10))
ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))
plt.xlabel("Low Frequency Cutoff")
plt.ylabel("High Frequency Cutoff")
plt.title("Parameter Estimation for Cosmic String Cusp Template")

# Add a color bar which maps values to colors.
fig.colorbar(surf, shrink=0.5, aspect=5)
plt.show()

"""
