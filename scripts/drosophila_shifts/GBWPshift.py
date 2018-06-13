#!/usr/bin/python3
#from pylab import *
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
import copy
import FlyFactory
from phhotoreceptor.DepolarisePhotoreceptor import DepolarisePhotoreceptor
import phhotoreceptor.Experiment as Experiment
from GBWPutils import GBWP
import ShiftConductances

rc('font',**{'family':'serif'}) #,'serif':['Liberation Serif']})

change_LIC_to_keep_depolarisation = False

HH  = FlyFactory.DrosophilaR16()
f_medium = 2 #Hz

fig1, ax_GBWP = plt.subplots(3, 1, figsize=(6,8))
label = '(a) (b) (c)'.split()
plt.subplots_adjust(hspace=0.45, left=0.1, right=0.95, bottom=0.05, top=0.95)

for i in range(3):
    ax_GBWP[i].tick_params(direction='in', top=True, right=True)
    ax_GBWP[i].set_ylim([3,4.5])
    ax_GBWP[i].set_xlim([-69,-32])
    ax_GBWP[i].text(-0.04, 1.2, label[i], transform=ax_GBWP[i].transAxes,
                    fontsize=14, va='top', ha='right')


Vr = np.arange(-68.0, -30.0, 8)
deltaV = 0.5
Vr_continuous = np.arange(-68, -36 + deltaV, deltaV)

colour_graph=['y','b','g','r','c']

##Continuous across depolarisations

GBWP_continuous_ = np.zeros_like(Vr_continuous)
GBWP_RC_continuous_ = np.zeros_like(Vr_continuous)
GBWP_shift_continuous_ = np.zeros((3,len(Vr_continuous)))
Vr_shift_continuous = np.zeros((3,len(Vr_continuous)))

for i,V in enumerate(Vr_continuous):
    DepolarisePhotoreceptor.WithLight(HH,V)
    GBWP_continuous_[i] = GBWP(HH.body.impedance, f_min = f_medium)
    Experiment.freeze_conductances(HH)
    GBWP_RC_continuous_[i] = GBWP(HH.body.impedance, f_min = f_medium)
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
        GBWP_shift_continuous_[ii,i] = GBWP(HH_shifted.body.impedance, f_min = f_medium)
        Vr_shift_continuous[ii,i] = HH_shifted.body.V_m

for ii in range(3):
    ax_GBWP[ii].plot(Vr_continuous,GBWP_continuous_/1e3,'k',zorder=0,alpha=0.2)
    ax_GBWP[ii].plot(Vr_shift_continuous[ii],GBWP_shift_continuous_[ii,:]/1e3,'k',zorder=0)
    ax_GBWP[ii].plot(Vr_continuous,GBWP_RC_continuous_/1e3,'k--',zorder=0,alpha=0.2)



GBWP_ = np.zeros_like(Vr)
GBWP_RC_ = np.zeros_like(Vr)
GBWP_shift_ = np.zeros((3,len(Vr)))

for i,V in enumerate(Vr):
    DepolarisePhotoreceptor.WithLight(HH,V)
    GBWP_[i] = GBWP(HH.body.impedance, f_min = f_medium)
    Experiment.freeze_conductances(HH)
    GBWP_RC_[i] = GBWP(HH.body.impedance, f_min = f_medium)
    Experiment.unfreeze_conductances(HH)
    HH_shifted = [copy.deepcopy(HH) for i in range(3)]

    for ii in range(3):
        if ii==0:
            ShiftConductances.WithLight(HH_shifted[ii], change_LIC_to_keep_depolarisation = change_LIC_to_keep_depolarisation)
        elif ii==2:
            ShiftConductances.WithSerotonin(HH_shifted[ii], change_LIC_to_keep_depolarisation = change_LIC_to_keep_depolarisation)
        elif ii==1:
            Experiment.modify_conductance(HH_shifted[ii], "Shab", .5, change_LIC_to_keep_depolarisation = change_LIC_to_keep_depolarisation)
        else:
            raise Exception

        GBWP_shift_[ii,i] = GBWP(HH_shifted[ii].body.impedance, f_min = f_medium)

    for ii in range(3):
        ax_GBWP[ii].plot(HH_shifted[ii].body.V_m ,GBWP_shift_[ii,i]/1e3,colour_graph[i] + '.',markersize=15)
        ax_GBWP[ii].plot(V,GBWP_[i]/1e3,colour_graph[i] + '.',markersize=15,alpha=0.5)
        ax_GBWP[ii].plot(V,GBWP_RC_[i]/1e3,colour_graph[i] + '.',markersize=15,alpha=0.5)

for ii in range(3):
    ax_GBWP[ii].set_ylabel(r'GBWP (G$\Omega$ Hz)')
ax_GBWP[2].set_xlabel('Membrane voltage (mV)')

plt.show()
