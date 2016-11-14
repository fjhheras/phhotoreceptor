#!/usr/bin/python3

from pylab import *
from numpy import *
import FlyFactory
from phhotoreceptor.DepolarisePhotoreceptor import DepolarisePhotoreceptor
import phhotoreceptor.Experiment as Experiment

option_debugging = False
depolarise_with_light = True #If depolarise with current, all cost calculations are not biological

HH  = FlyFactory.DrosophilaR16()
E_K = HH.body.voltage_channels[0].E
E_lL = HH.body.light_conductance.E


phasearray = vectorize (lambda z : angle(z))

def calculate_bandwidth_of_passive_photoreceptor(photoreceptor,low_limit_frequency = 0):
    C = photoreceptor.body.C
    R_a = photoreceptor.body.resistance()
    bw_when_low_freq_is_0 = 1000/( 2*pi*R_a*C) #Hz
    correction = 2*low_limit_frequency**2 #Hz2
    return sqrt(bw_when_low_freq_is_0**2 + correction)

def Calculate_Bandwidth(Z,f):
    Z_max = abs(Z[0])
    for i,ff in enumerate(f):
        if abs(Z[i]) > Z_max:
            Z_max = abs(Z[i])
        if Z_max/sqrt(2) > abs(Z[i]):
            return ff

def GBWP(Z,f):
    Z_max = abs(Z[0])
    for i,ff in enumerate(f):
        if abs(Z[i]) > Z_max:
            Z_max = abs(Z[i])
        if Z_max/sqrt(2) > abs(Z[i]):
            return ff*Z_max


####### BODY STARTS HERE

f_low = 0 #Hz
f_medium = 2 #Hz

fig1 = figure(1)
ax_GBWP=[]
for ii in range(3):
    ax_GBWP.append(fig1.add_subplot(3,1,ii+1))

### ONLY CERTAIN VOLTAGES BUT CONTINUOUS ACROSS FREQUENCIES


Vr=range(-68,-30,8)
Vr_continuous = arange(-68,-36,0.1)
delta_f = 0.01
#f = arange(.2,500,delta_f)
delta_V = 0.01
f_from_medium = arange(f_medium,300,delta_V)

colour_graph=['y','b','g','r','c']

GBWP_continuous_ = zeros_like(Vr_continuous)
GBWP_RC_continuous_ = zeros_like(Vr_continuous)
GBWP_selective_continuous_ = zeros((3,len(Vr_continuous)))

for i,V in enumerate(Vr_continuous):

    if depolarise_with_light:
        DepolarisePhotoreceptor.WithLight(HH,V)
    else:
        DepolarisePhotoreceptor.WithCurrent(HH,V)

    Z_ = HH.body.impedance(f_from_medium) #Only frequencies above f_medium, to calculate bandwidth
    Experiment.freeze_conductances(HH)
    Z_fixed_ = HH.body.impedance(f_from_medium)
    Experiment.unfreeze_conductances(HH)

    GBWP_continuous_[i] = GBWP(Z_,f_from_medium)/1e3
    GBWP_RC_continuous_[i] = GBWP(Z_fixed_,f_from_medium)/1e3

    for ii in range(3):
        for iii in range(3):
            if ii!=iii:
                Experiment.freeze_conductances(HH,index=iii)
        Z_ = HH.body.impedance(f_from_medium)
        GBWP_selective_continuous_[ii,i] = GBWP(Z_,f_from_medium)/1e3
        Experiment.unfreeze_conductances(HH)

for ii in range(3):
    ax_GBWP[ii].plot(Vr_continuous,GBWP_continuous_/1e3,'k',zorder=0,alpha=0.2)
    ax_GBWP[ii].plot(Vr_continuous,GBWP_selective_continuous_[ii,:]/1e3,'k',zorder=0)
    ax_GBWP[ii].plot(Vr_continuous,GBWP_RC_continuous_/1e3,'k--',zorder=0,alpha=0.2)





GBWP_ = zeros_like(Vr)
GBWP_RC_ = zeros_like(Vr)
GBWP_selective_ = zeros((3,len(Vr)))

for i,V in enumerate(Vr):

    if depolarise_with_light:
        DepolarisePhotoreceptor.WithLight(HH,V)
    else:
        DepolarisePhotoreceptor.WithCurrent(HH,V)

    Z_ = HH.body.impedance(f_from_medium) #Only frequencies above f_medium, to calculate bandwidth
    Experiment.freeze_conductances(HH)
    Z_fixed_ = HH.body.impedance(f_from_medium)
    Experiment.unfreeze_conductances(HH)

    GBWP_[i] = GBWP(Z_,f_from_medium)/1e3
    GBWP_RC_[i] = GBWP(Z_fixed_,f_from_medium)/1e3

    for ii in range(3):
        for iii in range(3):
            if ii!=iii:
                Experiment.freeze_conductances(HH,index=iii)
        Z_ = HH.body.impedance(f_from_medium)
        GBWP_selective_[ii,i] = GBWP(Z_,f_from_medium)/1e3
        Experiment.unfreeze_conductances(HH)

    for ii in range(3):
        ax_GBWP[ii].plot(V,GBWP_selective_[ii,i]/1e3,colour_graph[i] + '.',markersize=15)
        ax_GBWP[ii].plot(V,GBWP_[i]/1e3,colour_graph[i] + '.',markersize=15,alpha=0.5)
        ax_GBWP[ii].plot(V,GBWP_RC_[i]/1e3,colour_graph[i] + '.',markersize=15,alpha=0.5)

#for ii in range(3):
#    ax_GBWP[ii].plot(Vr,GBWP_,'k',zorder=0,alpha=0.2)
#    ax_GBWP[ii].plot(Vr,GBWP_selective_[ii,:],'k',zorder=0)
#    ax_GBWP[ii].plot(Vr,GBWP_RC_,'k--',zorder=0,alpha=0.2)




show()
