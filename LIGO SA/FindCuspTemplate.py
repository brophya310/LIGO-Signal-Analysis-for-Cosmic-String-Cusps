# Generate a template for detector strain
# as a function of frequency using
#eq.  (1) and eq (3)
# in Gravitational Wave Bursts from Cosmic Super Strings: Quantitative analysis and constraints
# by Aidan Brophy 2020

import readligo as rl
import numpy as np
import math
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt

#Constants

#G = 6.67E-11
#mu = 1 #mass per unit length of string
#r = 1 #Mpc
#L = 1 #size of the feature that produces the cusp
A = 1.00E-18 #scaling

#theta function for high and low frequency cutoff
def theta(x):
   if x<=0.0:
       return 0.0
   return 1.0


###########################################
#Define all functions which generate template
###########################################


#Calculate strain using frequency domain
def h(freqs,fl,fh):
    h_f = A*(abs((freqs))**((-4)/(3)))
    #(f-central frequency)
    #*theta((f-fl))*theta((fh-f))
    for i in np.arange(len(freqs)):
      f = freqs[i]
      if fl >= f or  f >= fh:
         h_f[i] = 0
    return h_f

def create_a_template(fs,dataChunk, fl, fh):
    
    dt=dataChunk/2
    # Create array of frequencies
    nyquist = fs/2
    f_i = 1./(2*dt)
    frequency = np.arange(0,nyquist+1./(2*dt),1./(2*dt))
    frequency[0] = 1./(2*dt)

    #use frequencies to find strain:
    strain = h(frequency,fl,fh)

    return strain,frequency

###########################
#Testing template
###########################



