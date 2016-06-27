#!/usr/bin/python3
from __future__ import division

from numpy import *
from matplotlib import pyplot as plt

import Drone as Drone
from phhotoreceptor.DepolarisePhotoreceptor import DepolarisePhotoreceptor
import phhotoreceptor.Experiment as Experiment
import phhotoreceptor.Linearise as Linearise

__author__ = 'Francisco J. H. Heras'

T=300 #ms
dt=0.05 #ms
time_array = arange(0,T+dt,dt)
I = zeros_like(time_array)

fig1 = plt.figure(1)
ax = fig1.add_subplot(211)
ax_t = fig1.add_subplot(212)

fig2 = plt.figure(2)
ax_curr = fig2.add_subplot(211)

V_membrane = -38 #mV
photoreceptor = Drone.Vallet92()
DepolarisePhotoreceptor.WithLight(photoreceptor, V = V_membrane)

fig1.text(0.06, 0.5, 'Membrane potential deflection (mV)', ha='center', va='center', rotation='vertical')

for ii in range(5) :
    for i, t in enumerate(time_array):
        if 10 <= t <= 110: I[i] = 1e-3*(-0.02+0.01*ii)  # nA->uA
    DepolarisePhotoreceptor.WithLight(photoreceptor, V = V_membrane)

    V_array, g_Ch = Experiment.inject_current(photoreceptor,I,dt)
    ax_curr.plot(time_array, I, 'k')
    ax.plot(time_array, V_array-V_membrane,color='black') #mV
    ax.set_xticklabels([])

    V_array = Linearise.inject_current(photoreceptor,I,dt)
    ax.plot(time_array, V_array,'k--') #mV
    ax.set_xticklabels([])

V_membrane = -38 #mV
photoreceptor = Drone.Vallet92(k_h = 1)
DepolarisePhotoreceptor.WithLight(photoreceptor, V = V_membrane)

for ii in range(5) :
    for i, t in enumerate(time_array):
        if 10 <= t <= 110: I[i] = 1e-3*(-0.02+0.01*ii)  # nA->uA
    DepolarisePhotoreceptor.WithLight(photoreceptor, V = V_membrane)
    V_array, g_Ch = Experiment.inject_current(photoreceptor,I,dt)
    ax_t.plot(time_array, V_array-V_membrane,color='black') #mV
    V_array = Linearise.inject_current(photoreceptor,I,dt)

ax_t.plot(time_array, V_array, 'k--')  # mV
ax_t.set_xticklabels([])
ax_t.set_xlabel('Time (ms)')

ax_curr.set_xlabel('Time (ms)')
ax_curr.set_ylabel('Current (uA)')

plt.show()