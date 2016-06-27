#!/usr/bin/python3
from __future__ import division

from numpy import *
from matplotlib import pyplot as plt

from phhotoreceptor.DepolarisePhotoreceptor import DepolarisePhotoreceptor
import phhotoreceptor.Experiment as Experiment
import phhotoreceptor.Linearise as Linearise
import FlyFactory as FlyFactory

__author__ = 'Francisco J. H. Heras'

T=200 #ms
dt=0.05 #ms
time_array = arange(0,T+dt,dt)
I = zeros_like(time_array)

fig1 = plt.figure(1)
ax = fig1.add_subplot(211)
ax_t = fig1.add_subplot(212)

V_membrane = -60 #mV
photoreceptor = FlyFactory.CalliphoraR16(channel_choice = "Weckstrom")
DepolarisePhotoreceptor.WithLight(photoreceptor, V = V_membrane)

fig1.text(0.06, 0.5, 'Membrane potential deflection (mV)', ha='center', va='center', rotation='vertical')

for ii in range(4) :
    for i, t in enumerate(time_array):
        if 20 <= t <= 120: I[i] = 1e-3*(-0.075+0.05*ii)  # nA->uA
    DepolarisePhotoreceptor.WithLight(photoreceptor, V = V_membrane)

    V_array, g_Ch = Experiment.inject_current(photoreceptor,I,dt)
    ax.plot(time_array, V_array-V_membrane,color='black') #mV
    ax.set_xticklabels([])

    V_array = Linearise.inject_current(photoreceptor,I,dt)
    ax.plot(time_array, V_array,'k--') #mV
    ax.set_xticklabels([])

I = zeros_like(time_array)
V_membrane = -40 #mV
photoreceptor = FlyFactory.CalliphoraR16(channel_choice = "Weckstrom")
DepolarisePhotoreceptor.WithLight(photoreceptor, V = V_membrane)

for ii in [0,2,3,4,5,7] :
    for i, t in enumerate(time_array):
        if 20 <= t <= 120: I[i] = 1e-3*(-0.7+0.2*ii)  # nA->uA
    DepolarisePhotoreceptor.WithLight(photoreceptor, V = V_membrane)

    V_array, g_Ch = Experiment.inject_current(photoreceptor,I,dt)
    ax_t.plot(time_array, V_array-V_membrane,color='black') #mV
    ax_t.set_xlabel('Time (ms)')

    V_array = Linearise.inject_current(photoreceptor,I,dt)
    ax_t.plot(time_array, V_array,'k--') #mV
    ax_t.set_xticklabels([])

plt.show()