#!/usr/bin/python3

from pylab import *
from numpy import *
import FlyFactory as FlyFactory
from phhotoreceptor.DepolarisePhotoreceptor import DepolarisePhotoreceptor
import phhotoreceptor.Experiment as Experiment

__author__ = 'Francisco J. H. Heras'

photoreceptor  = FlyFactory.CalliphoraR16(channel_choice = "Weckstrom")
V = -60
tau_array = [0.1,3,8]
f_plot = logspace(0,2.5)

T=100 #ms
dt=0.05 #ms
time_array = arange(0,T+dt,dt)
I = zeros_like(time_array)

fig1 = figure(1)
ax1 = fig1.add_subplot(122)
ax2 = fig1.add_subplot(121)


DepolarisePhotoreceptor.WithLight(photoreceptor,V)
tau_fast_original = photoreceptor.body.voltage_channels[0].m_time(V)

for i, t in enumerate(time_array):
    if 0.01 <= t <= 100: I[i] = 1e-3*(0.1)  # nA->uA

for tau_fast in tau_array:
    DepolarisePhotoreceptor.WithLight(photoreceptor,V)
    photoreceptor.body.voltage_channels[0].time_multiplier = tau_fast/tau_fast_original
    Z = photoreceptor.body.impedance(f_plot)
    V_array, g_Ch = Experiment.inject_current(photoreceptor,I,dt)
    ax1.loglog(f_plot,abs(Z)/1000,'k')
    ax1.set_ylim(ymin = 4, ymax = 100)
    ax1.set_xlabel("Frequency (Hz)")
    ax1.set_ylabel("Impedance (MOhm)")
    ax2.plot(time_array,V_array,'k')
    ax2.set_xlabel("Time (ms)")
    ax2.set_ylabel("Membrane voltage (mV)")

show()
