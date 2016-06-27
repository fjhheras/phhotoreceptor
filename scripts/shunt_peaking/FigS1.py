#!/usr/bin/python3
from __future__ import division

from numpy import *
from matplotlib import pyplot as plt
from phhotoreceptor.DepolarisePhotoreceptor import DepolarisePhotoreceptor
import phhotoreceptor.Experiment as Experiment
import FlyFactory as FlyFactory


__author__ = 'Francisco J. H. Heras'


T=400 #ms
dt=0.1 #ms
time_array = arange(0,T+dt,dt)
I = zeros_like(time_array)

option = 1 # If 1 -> S1b, else -> S1a

if (option ==1):
    N_rep = 8
    a = -0.15
    b = 0.3
else:
    N_rep = 6
    a = -0.5
    b = 0.2

V_membrane = -60 #mV

photoreceptor = FlyFactory.CalliphoraR16(channel_choice = "Weckstrom")

for ii in range(N_rep) :
    for i, t in enumerate(time_array):
        if 10 <= t <= 210: I[i] = 1e-3*(a+b*ii)  # nA->uA  --- Fig S1b
    DepolarisePhotoreceptor.WithLight(photoreceptor, V = V_membrane)

    V_array, g_Ch = Experiment.inject_current(photoreceptor,I,dt)
    plt.ylabel_set = False
    ax = plt.subplot(2,1,1)
    ax.plot(time_array, V_array,color='black') #mV
    ax.set_title('Hodgkin-Huxley voltage (top) and conductances (bottom)')
    plt.ylabel('Potential (mV)')
    ax.set_xticklabels([])
    ax = plt.subplot(2,1,2)
    ax.plot(time_array, g_Ch[0]*1e6,color='blue') #Fast conductance, nS
    ax.plot(time_array, g_Ch[1]*1e6,color='red') #Slow conductance, nS
    plt.xlabel('Time (msec)')
    plt.ylabel('Conductances (nS)')

plt.show()