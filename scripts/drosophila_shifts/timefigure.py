#!/usr/bin/python3
from __future__ import division

from numpy import *
from matplotlib.pyplot import show,figure
import FlyFactory
from phhotoreceptor.DepolarisePhotoreceptor import DepolarisePhotoreceptor
import phhotoreceptor.Experiment as Experiment

__author__ = 'Francisco J. H. Heras'

option = 1 # 1:Compare with full conductance freeze, 2: Compare with inactivation freeze

V_membrane_ = [-68,-59,-41] #mV
plot_window = array([-.5,3.5])
drosophila = FlyFactory.DrosophilaR16()
fig1 = figure(1)

T=200 #ms
dt=0.5 #ms
time_array = arange(0,T+dt,dt)
I = zeros_like(time_array)

for ii,V_membrane in enumerate(V_membrane_) :
    for i, t in enumerate(time_array):
        if 10 <= t <= 160: I[i] = 1e-3*(0.01)  # nA->uA
    DepolarisePhotoreceptor.WithLight(drosophila, V = V_membrane)

    V_array, g_Ch = Experiment.inject_current(drosophila,I,dt)
    ax = fig1.add_subplot(3,1,ii+1)
    ax.plot(time_array, V_array,color='black') #mV

    #Experiment.unfreeze_conductances(drosophila)

    DepolarisePhotoreceptor.WithLight(drosophila, V = V_membrane) #To make sure that all channels are back at rest

    Experiment.freeze_inactivations(drosophila)
    V_array, g_Ch = Experiment.inject_current(drosophila,I,dt)
    Experiment.unfreeze_inactivations(drosophila)
    ax.plot(time_array, V_array,'k:') #mV

    DepolarisePhotoreceptor.WithLight(drosophila, V = V_membrane) #To make sure that all channels are back at rest

    Experiment.freeze_conductances(drosophila)
    V_array, g_Ch = Experiment.inject_current(drosophila,I,dt)
    Experiment.unfreeze_conductances(drosophila)
    ax.plot(time_array, V_array,'k--') #mV

    ax.set_ylim(plot_window + V_membrane)

show()
