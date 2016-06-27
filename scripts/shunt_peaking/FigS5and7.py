#!/usr/bin/python3

from pylab import *
from numpy import *
import FlyFactory as FlyFactory
from phhotoreceptor.DepolarisePhotoreceptor import DepolarisePhotoreceptor
import phhotoreceptor.Experiment as Experiment
from GBWPutils import GBWP, Is_Band_Pass, GDD, GDV

__author__ = 'Francisco J. H. Heras'

####### BODY STARTS HERE

HH  = FlyFactory.CalliphoraR16(channel_choice = "Weckstrom")
V = -40 # Membrane voltage in mV
f_down = 1 #Lower frequency where the group delay is calculated, Hz

n = 50
tau_fast = linspace(0.1,8,n) #ms

### FIGURES

fig1 = figure(1) # Plot of GBWP vs FDR activation time constant
ax_gbwp1 = fig1.add_subplot(111)

fig2 = figure(2)
ax_gdv = fig2.add_subplot(111)

### INITIALISE

GBWP_a = zeros_like(tau_fast)
GBWPf_a = zeros_like(tau_fast)
GDV_a = zeros_like(tau_fast)
GDVf_a = zeros_like(tau_fast)

### RUN

DepolarisePhotoreceptor.WithLight(HH,V)
passive_gbwp = 1/2/pi/HH.body.C
tau_fast_original = HH.body.voltage_channels[0].m_time(V)
print("Passive GBWP is ", passive_gbwp, " MOhm Hz")
print("Original time constant for FDR is ", tau_fast_original, "ms")

low_pass_found = False
for i,tau in enumerate(tau_fast):
    HH.body.voltage_channels[0].time_multiplier = tau_fast[i]/tau_fast_original
    GBWP_a[i]= GBWP(HH.body.impedance)/passive_gbwp
    GDV_a[i] = GDV(HH.body.impedance, f_down=f_down)
    if (Is_Band_Pass(HH.body.impedance) and low_pass_found==False):
        print ("First low pass when tau = ", tau)
        low_pass_found = True

Experiment.freeze_conductances(HH,index=1) #freeze SDR

low_pass_found = False
for i,tau in enumerate(tau_fast):
    HH.body.voltage_channels[0].time_multiplier = tau_fast[i]/tau_fast_original
    GBWPf_a[i]= GBWP(HH.body.impedance)/passive_gbwp
    GDVf_a[i] = GDV(HH.body.impedance, f_down=f_down)
    if (Is_Band_Pass(HH.body.impedance) and low_pass_found==False):
        print ("First low pass when tau = ", tau)
        low_pass_found = True


Experiment.freeze_conductances(HH,index=0) #freeze FDR and SDR
passive_gdv = GDV(HH.body.impedance)

ax_gbwp1.plot(tau_fast,GBWP_a,'k')
ax_gbwp1.plot(tau_fast,GBWPf_a,'k--')
ax_gbwp1.vlines(x=tau_fast_original, ymin=1.001,ymax=max(GBWP_a),colors='r')
ax_gdv.plot(tau_fast,GDV_a,'k')
ax_gdv.plot(tau_fast,GDVf_a,'k--')
ax_gdv.plot(tau_fast,passive_gdv*ones_like(tau_fast),'k:')

ax_gbwp1.set_xlabel("FDR activation time constant (ms)")
ax_gbwp1.set_ylabel("Relative GBWP")

ax_gdv.set_xlabel("FDR activation time constant (ms)")
ax_gdv.set_ylabel("Std of group delay (ms)")

show()
