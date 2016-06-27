#!/usr/bin/python3

from pylab import *
from numpy import *
import FlyFactory as FlyFactory
from phhotoreceptor.DepolarisePhotoreceptor import DepolarisePhotoreceptor
from GBWPutils import GBWP

HH  = FlyFactory.CalliphoraR16(channel_choice = "Weckstrom")
V = -60 #steady state voltage: -60 -> Fig 1e,  -40 -> Fig 1f
n = 100 #number of divisions

####### BODY STARTS HERE

tau_fast = linspace(.10,6,n) #
tau_slow = linspace(.10,50,n) #
Tau_fast, Tau_slow = meshgrid(tau_fast,tau_slow)

fig1 = figure(1)
ax_gbwp = fig1.add_subplot(111)

DepolarisePhotoreceptor.WithLight(HH,V)

GBWP1_a = zeros_like(Tau_slow)

passive_gbwp = 1/2/pi/HH.body.C
tau_fast_original = HH.body.voltage_channels[0].m_time(V)
tau_slow_original = HH.body.voltage_channels[1].m_time(V)
print("Original time constant for FDR is ", tau_fast_original, "ms, and SDR's is ", tau_slow_original)

it = nditer(Tau_slow, flags=['multi_index'])
max_GBWP1 = 0
while not it.finished:
    HH.body.voltage_channels[0].time_multiplier = Tau_fast[it.multi_index]/tau_fast_original
    HH.body.voltage_channels[1].time_multiplier = Tau_slow[it.multi_index]/tau_slow_original
    GBWP1_a[it.multi_index]= GBWP(HH.body.impedance)/passive_gbwp
    if (GBWP(HH.body.impedance)/passive_gbwp > max_GBWP1):
        max_tau_fast = Tau_fast[it.multi_index]
        max_tau_slow = Tau_slow[it.multi_index]
        max_GBWP1 = GBWP(HH.body.impedance)/passive_gbwp
    it.iternext()

CS = ax_gbwp.contour(Tau_fast,Tau_slow,GBWP1_a)
clabel(CS, inline=1, fontsize=10)
ax_gbwp.plot(tau_fast_original, tau_slow_original, '+')
ax_gbwp.set_xlabel("Time constant fast channel (ms)")
ax_gbwp.set_ylabel("Time constant slow channel (ms)")
print("Maximum GBWP1 obtained when FDR and SDR a. time constant are respectively ", max_tau_fast, " ms and ", max_tau_slow, " ms")

show()
