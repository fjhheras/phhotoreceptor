#!/usr/bin/python3
from pylab import *
from numpy import *
import Drone as Drone
from phhotoreceptor.DepolarisePhotoreceptor import DepolarisePhotoreceptor
from GBWPutils import GBWP,Is_Band_Pass

HH  = Drone.Vallet92() #Modified drone photoreceptor

####### BODY STARTS HERE

n = 100
tau_m = linspace(.01,3,n) #.01,3
tau_h = linspace(.01,15,n) #.01,15
extent=(tau_m[0],tau_m[-1],tau_h[0],tau_h[-1])
Tau_m, Tau_h = meshgrid(tau_m,tau_h)

fig1 = figure(1) # Plot comparing against RC
ax_gbwp1 = fig1.add_subplot(111)
fig2 = figure(2)
ax_low_pass = fig2.add_subplot(111)

V = -38
DepolarisePhotoreceptor.WithLight(HH,V)

GBWP1_a = zeros_like(Tau_h)

Band_pass_a = zeros_like(Tau_h)
passive_gbwp = 1/2/pi/HH.body.C
tau_h_original = HH.body.voltage_channels[0].h_time(V)
tau_m_original = HH.body.voltage_channels[0].m_time(V)

print("Original time constant for activation is ", tau_m_original, "ms, and inactivation is ", tau_h_original)

it = nditer(Tau_h, flags=['multi_index'])
max_GBWP1 = 0
while not it.finished:
    HH.body.voltage_channels[0].h_time_multiplier = Tau_h[it.multi_index]/tau_h_original
    HH.body.voltage_channels[0].m_time_multiplier = Tau_m[it.multi_index]/tau_m_original
    GBWP1_a[it.multi_index]= GBWP(HH.body.impedance)/passive_gbwp
    Band_pass_a[it.multi_index] = Is_Band_Pass(HH.body.impedance)
    if (GBWP(HH.body.impedance)/passive_gbwp > max_GBWP1):
        max_tau_m = Tau_m[it.multi_index]
        max_tau_h = Tau_h[it.multi_index]
        max_GBWP1 = GBWP(HH.body.impedance)/passive_gbwp
    it.iternext()

ax_low_pass.imshow(Band_pass_a, interpolation='bicubic', cmap = get_cmap('gray'), extent = extent, origin = 'Lower', aspect='auto')


CS = ax_gbwp1.contour(Tau_m,Tau_h,GBWP1_a)
clabel(CS, inline=1, fontsize=10)

ax_gbwp1.set_xlabel("Time constant of activation (ms)")
ax_gbwp1.set_ylabel("Time constant of inactivation (ms)")

show()
