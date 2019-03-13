# -*- coding: utf-8 -*-
"""
Created on Thu Dec  6 15:30:20 2018

@author: Martin Shahi


Calculating the order parameter of the system using Hilbert transform.
Order parameter needs to be calculated for every neuron for every time value.
"""

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.signal import hilbert

df = pd.read_csv("FHN.csv")
N = 9 #number of neurons

phases = []
system_data = []
rho_t = []
psi_t = []
real = []
imaginary = []

def orderParameter(data):
    #Calculate Hilbert transform of time series to create analytic signal,
         #then find the instantaneous phase of time series for each t, then make sure
         #phase is within 0->2pi
     hilbert_data = hilbert(data)
     return [np.unwrap(np.angle(hilbert_data)),hilbert_data]
     
         
for i in range(N):
    
    phases.append(orderParameter(df['v%d'%(i+1)])[0])

for i in range(len(phases[0])):
    temp = 0
    for j in range(N):
        temp += np.exp(1j*phases[j][i])
    temp = temp/N
    system_data.append(temp)
for i in system_data:
    rho_t.append(np.abs(i))
    psi_t.append(np.angle(i))
    real.append(i.real)
    imaginary.append(i.imag)
    
    
plt.scatter(real,imaginary)
plt.show()
