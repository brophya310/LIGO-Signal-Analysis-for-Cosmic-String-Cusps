#Template Sounds

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
import readligo as rl
import template 
import FindCuspTemplate
from scipy.io import wavfile
from playsound import playsound

strain, strain_times, dq = rl.loaddata("L-L1_GWOSC_4KHZ_R1-1128676853-4096.hdf5")
#if len(strain) < 2:
    #print("Failed to read file.")
    #exit()

dt = strain_times[1] - strain_times[0]
fs = 1.0 / dt


#Creating an Injection Segment from GPS Time
#given the gps time of the injection
gpstime_event =  1128678900
start = gpstime_event - 10
stop = gpstime_event + 10


#seg = rl.getsegs(start, stop, 'H1')[0]
#for (begin, end) in segs:
seg_strains, meta, dq = rl.getstrain(start, stop, 'L1')
seg_duration = len(seg_strains)/fs
seg_times = np.arange(0, len(seg_strains))*dt + start

#create template, frequency domain
temp, temp_freq = template.createTemplate(fs, seg_duration, 23, 13)
ctemp, ctemp_freq = FindCuspTemplate.create_a_template(fs, seg_duration, 25, 220)

#set template to zero at low frequencies
temp[temp_freq < 25] = 0
ctemp[ctemp_freq < 25] = 0

#template as time series
t_temp = abs(temp)
t = np.arange(len(t_temp))*dt
inv_cusp_temp = np.fft.ifft(ctemp)
invtemp = np.fft.ifft(temp)

#roll to the center
middleIndex = (len(inv_cusp_temp) - 1)/2
middleIndex2 = (len(invtemp)-1)/2
inv_cusp_tempT = np.roll(inv_cusp_temp.real, -int(middleIndex))
invtempT = np.roll(invtemp.real, -int(middleIndex2))

CBC_sound = []
CBC_sound = invtempT[73764:90156]
CUSP_sound = []
CUSP_sound = inv_cusp_tempT[73764:90156]

#Making Sounds for Both Templates:

def write_wavfile(filename,fs,data):
    d = np.int16(data/np.max(np.abs(data)) * 32767 * 0.9)
    wavfile.write(filename,int(fs), d)

#deltat_sound = 2.
#indxd = np.where((time >= tevent-deltat_sound) & (time < tevent+deltat_sound))

def reqshift(data,fshift=100,sample_rate=4096):
    """Frequency shift the signal by constant
    """
    x = np.fft.rfft(data)
    T = len(data)/float(sample_rate)
    df = 1.0/T
    nbins = int(fshift/df)
    # print T,df,nbins,x.real.shape
    y = np.roll(x.real,nbins) + 1j*np.roll(x.imag,nbins)
    y[0:nbins]=0.
    z = np.fft.irfft(y)
    return z


fs = 4096
fshift = 400.
speedup = 1.
fss = int(float(fs)*float(speedup))

# shift frequency of the data
CBC_shifted = reqshift(CBC_sound,fshift=fshift,sample_rate=fs)
CUSP_shifted = reqshift(CUSP_sound,fshift=fshift,sample_rate=fs)

# write the files:
write_wavfile("CBC_shifted.wav",int(fs), CBC_shifted)
write_wavfile("CUSP_shifted.wav",int(fs), CUSP_shifted)

#playsound("CBC_shifted.wav")
playsound("CUSP_shifted.wav")


