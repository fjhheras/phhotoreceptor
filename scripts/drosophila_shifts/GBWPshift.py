#!/usr/bin/python3
#from pylab import *
import numpy as np
import matplotlib.pyplot as plt
import copy
import FlyFactory
from phhotoreceptor.DepolarisePhotoreceptor import DepolarisePhotoreceptor
import phhotoreceptor.Experiment as Experiment
from GBWPutils import GBWP
import ShiftConductances

option_debugging = False
depolarise_with_light = True #If depolarise with current, all cost calculations are not biological
change_LIC_to_keep_depolarisation = True

HH  = FlyFactory.DrosophilaR16()
f_medium = 2 #Hz

fig1, ax_GBWP = plt.subplots(3,1)
plt.subplots_adjust(hspace=0.45)

for i in range(3):
    ax_GBWP[i].tick_params(direction='in')
    ax_GBWP[i].set_ylim([3,4.5])

Vr = np.arange(-68.0, -30.0, 8)
deltaV = 0.5
Vr_continuous = np.arange(-68, -36 + deltaV, deltaV)

colour_graph=['y','b','g','r','c']

##Continuous across depolarisations

GBWP_continuous_ = np.zeros_like(Vr_continuous)
GBWP_RC_continuous_ = np.zeros_like(Vr_continuous)
GBWP_shift_continuous_ = np.zeros((3,len(Vr_continuous)))

for i,V in enumerate(Vr_continuous):

    if depolarise_with_light:
        DepolarisePhotoreceptor.WithLight(HH,V)
    else:
        DepolarisePhotoreceptor.WithCurrent(HH,V)

    GBWP_continuous_[i] = GBWP(HH.body.impedance, f_min = f_medium) #  GBWP(Z_,f_from_medium)/1e3
    Experiment.freeze_conductances(HH)
    GBWP_RC_continuous_[i] = GBWP(HH.body.impedance, f_min = f_medium) #GBWP(Z_fixed_,f_from_medium)/1e3
    Experiment.unfreeze_conductances(HH)

    for ii in range(3):
        HH_shifted = copy.deepcopy(HH)
        if ii==0:
            ShiftConductances.WithLight(HH_shifted, change_LIC_to_keep_depolarisation = change_LIC_to_keep_depolarisation)
        elif ii==2:
            ShiftConductances.WithSerotonin(HH_shifted, change_LIC_to_keep_depolarisation = change_LIC_to_keep_depolarisation)
        elif ii==1:
            Experiment.modify_conductance(HH_shifted, "Shab", .5, change_LIC_to_keep_depolarisation = change_LIC_to_keep_depolarisation)
        else:
            raise Exception
        GBWP_shift_continuous_[ii,i] = GBWP(HH_shifted.body.impedance, f_min = f_medium) #Z_,f_from_medium)/1e3

for ii in range(3):
    ax_GBWP[ii].plot(Vr_continuous,GBWP_continuous_/1e3,'k',zorder=0,alpha=0.2)
    ax_GBWP[ii].plot(Vr_continuous,GBWP_shift_continuous_[ii,:]/1e3,'k',zorder=0)
    ax_GBWP[ii].plot(Vr_continuous,GBWP_RC_continuous_/1e3,'k--',zorder=0,alpha=0.2)



GBWP_ = np.zeros_like(Vr)
GBWP_RC_ = np.zeros_like(Vr)
GBWP_shift_ = np.zeros((3,len(Vr)))

for i,V in enumerate(Vr):

    if depolarise_with_light:
        DepolarisePhotoreceptor.WithLight(HH,V)
    else:
        DepolarisePhotoreceptor.WithCurrent(HH,V)

    GBWP_[i] = GBWP(HH.body.impedance, f_min = f_medium) #Z_,f_from_medium)/1e3
    Experiment.freeze_conductances(HH)
    GBWP_RC_[i] = GBWP(HH.body.impedance, f_min = f_medium) #Z_fixed_,f_from_medium)/1e3
    Experiment.unfreeze_conductances(HH)

    for ii in range(3):
        HH_shifted = copy.deepcopy(HH)
        if ii==0:
            ShiftConductances.WithLight(HH_shifted, change_LIC_to_keep_depolarisation = change_LIC_to_keep_depolarisation)
        elif ii==2:
            ShiftConductances.WithSerotonin(HH_shifted, change_LIC_to_keep_depolarisation = change_LIC_to_keep_depolarisation)
        elif ii==1:
            Experiment.modify_conductance(HH_shifted, "Shab", .5, change_LIC_to_keep_depolarisation = change_LIC_to_keep_depolarisation)
        else:
            raise Exception

        GBWP_shift_[ii,i] = GBWP(HH_shifted.body.impedance, f_min = f_medium) #Z_,f_from_medium)/1e3

    for ii in range(3):
        ax_GBWP[ii].plot(V,GBWP_shift_[ii,i]/1e3,colour_graph[i] + '.',markersize=15)
        ax_GBWP[ii].plot(V,GBWP_[i]/1e3,colour_graph[i] + '.',markersize=15,alpha=0.5)
        ax_GBWP[ii].plot(V,GBWP_RC_[i]/1e3,colour_graph[i] + '.',markersize=15,alpha=0.5)

for ii in range(3):
    ax_GBWP[ii].set_ylabel('GBWP (GOhm Hz)')
ax_GBWP[2].set_xlabel('Membrane voltage (mV)')

plt.show()
