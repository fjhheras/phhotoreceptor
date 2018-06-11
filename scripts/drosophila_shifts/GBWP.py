#!/usr/bin/python3

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
import FlyFactory
from phhotoreceptor.DepolarisePhotoreceptor import DepolarisePhotoreceptor
import phhotoreceptor.Experiment as Experiment
from GBWPutils import GBWP

rc('font',**{'family':'serif'}) #,'serif':['Liberation Serif']})

option_debugging = False
depolarise_with_light = True #If depolarise with current, all cost calculations are not biological

HH  = FlyFactory.DrosophilaR16()
f_medium = 2 #Hz

fig1, ax_GBWP = plt.subplots(3, 1, figsize=(6,8))
label = '(a) (b) (c)'.split()
plt.subplots_adjust(hspace=0.45, left=0.1, right=0.95, bottom=0.05, top=0.95)
for i in range(3):
    ax_GBWP[i].tick_params(direction='in', top=True, right=True)
    ax_GBWP[i].set_ylim([2.75,4.25])
    ax_GBWP[i].text(-0.04, 1.2, label[i], transform=ax_GBWP[i].transAxes,
                    fontsize=14, va='top', ha='right')


#fig1 = plt.figure(1)
#ax_GBWP=[]
#for ii in range(3):
#    ax_GBWP.append(fig1.add_subplot(3,1,ii+1))

Vr=np.arange(-68.0,-30.0,8)
deltaV = 0.5
Vr_continuous = np.arange(-68,-36 + deltaV, deltaV)

colour_graph=['y','b','g','r','c']

GBWP_continuous_ = np.zeros_like(Vr_continuous)
GBWP_RC_continuous_ = np.zeros_like(Vr_continuous)
GBWP_selective_continuous_ = np.zeros((3,len(Vr_continuous)))

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
        for iii in range(3):
            if ii!=iii:
                Experiment.freeze_conductances(HH,index=iii)
        GBWP_selective_continuous_[ii,i] = GBWP(HH.body.impedance, f_min = f_medium) #Z_,f_from_medium)/1e3
        Experiment.unfreeze_conductances(HH)

for ii in range(3):
    ax_GBWP[ii].plot(Vr_continuous,GBWP_continuous_/1e3,'k',zorder=0,alpha=0.2)
    ax_GBWP[ii].plot(Vr_continuous,GBWP_selective_continuous_[ii,:]/1e3,'k',zorder=0)
    ax_GBWP[ii].plot(Vr_continuous,GBWP_RC_continuous_/1e3,'k--',zorder=0,alpha=0.2)





GBWP_ = np.zeros_like(Vr)
GBWP_RC_ = np.zeros_like(Vr)
GBWP_selective_ = np.zeros((3,len(Vr)))

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
        for iii in range(3):
            if ii!=iii:
                Experiment.freeze_conductances(HH,index=iii)
        GBWP_selective_[ii,i] = GBWP(HH.body.impedance, f_min = f_medium) #Z_,f_from_medium)/1e3
        Experiment.unfreeze_conductances(HH)

    for ii in range(3):
        ax_GBWP[ii].plot(V,GBWP_selective_[ii,i]/1e3,colour_graph[i] + '.',markersize=15)
        ax_GBWP[ii].plot(V,GBWP_[i]/1e3,colour_graph[i] + '.',markersize=15,alpha=0.5)
        ax_GBWP[ii].plot(V,GBWP_RC_[i]/1e3,colour_graph[i] + '.',markersize=15,alpha=0.5)


for ii in range(3):
    ax_GBWP[ii].set_ylabel('GBWP (GOhm Hz)')
ax_GBWP[2].set_xlabel('Membrane voltage (mV)')

plt.show()
