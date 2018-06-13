#!/usr/bin/python3

import matplotlib.pyplot as plt
from matplotlib import rc
from numpy import *
import FlyFactory
from phhotoreceptor.DepolarisePhotoreceptor import DepolarisePhotoreceptor
from GBWPutils import GBWP, Gain_Bandwidth

rc('font',**{'family':'serif'}) #,'serif':['Liberation Serif']})

depolarise_with_light = True #If depolarise with current, all cost calculations are not biological

####### BODY STARTS HERE

f_medium = 2 #2Hz
delta_cost = 1e-8 #If change, change also labels
calculate_bandwidth_without_inactivation = False

fig1, axes = plt.subplots(1, 3, figsize=(8,3))
label = '(a) (b) (c)'.split()
plt.tight_layout()
for i in range(3):
    axes[i].tick_params(direction='in', top=True, right=True)
    axes[i].text(-0.11, 1.04, label[i], transform=axes[i].transAxes,
                    fontsize=14, va='top', ha='right')


ax_no_inactivation_bw, ax_no_inactivation_cost, ax_no_inactivation_combined = axes


fig2 = plt.figure(6,figsize=[9,5])
ax_Z = fig2.add_subplot(1,1,1)


###### CONTINUOUS ACROSS VOLTAGES
Vr = arange(-68.0,-30.5,0.5)
Z_array = zeros(len(Vr))
Z_fixed_array = zeros(len(Vr))
total_resistance_array = zeros_like(Vr)
total_K_conductance = zeros_like(Vr)

HH = FlyFactory.DrosophilaR16(shift="none")
for i,V in enumerate(Vr): #Impedances at lowest frequency
    if depolarise_with_light:
        DepolarisePhotoreceptor.WithLight(HH,V)
    else:
        DepolarisePhotoreceptor.WithCurrent(HH,V)
    total_resistance_array[i] = HH.body.resistance()
    total_K_conductance[i] = HH.body.total_K_conductance()

### ONLY CERTAIN VOLTAGES BUT CONTINUOUS ACROSS FREQUENCIES

Vr=arange(-68.0,-30.5,8)
delta_f = 0.1
f = arange(1.5,900,delta_f)

colour_graph=['y','b','g','r','c']
Bandwidth = zeros_like(Vr)
Cost = zeros_like(Vr)
Cost_no_inactivation = zeros_like(Vr)
HH_no_inactivation = []

for i,V in enumerate(Vr):
    HH = FlyFactory.DrosophilaR16(shift="none")
    DepolarisePhotoreceptor.WithLight(HH,V,verbose=2)
    C = HH.body.C
    Z = HH.body.impedance(f) #All frequencies
    Cost[i] = HH.energy_consumption()*delta_cost
    if calculate_bandwidth_without_inactivation:
        HH.body.freeze_inactivation(index = None)
        _, Bandwidth[i] = Gain_Bandwidth(HH.body.impedance, f_min=f_medium)
        HH.body.unfreeze_all_conductances()
    else:
        _, Bandwidth[i] = Gain_Bandwidth(HH.body.impedance, f_min=f_medium)
    HH.body.freeze_inactivation(index = None) #Shaker 0, Shab 1
    HH_no_inactivation.append(HH)
    Cost_no_inactivation[i] = HH_no_inactivation[i].energy_consumption()*delta_cost
    Z_no_inactivation = HH_no_inactivation[i].body.impedance(f)#*HH.body.impedance(f_medium)/HH_no_inactivation[i].body.impedance(f_medium) #All frequencies


    label_str = str(V) + ' mV'
    ax_Z.loglog(f,abs(Z)/1000,colour_graph[i],linewidth=2,label = label_str)
    ax_Z.loglog(f,abs(Z_no_inactivation)/1000,colour_graph[i]+':',linewidth=1)

    ax_no_inactivation_bw.plot(V,Bandwidth[i],colour_graph[i] + '.',markersize=15)
    ax_no_inactivation_cost.plot(V,Cost[i],colour_graph[i] + '.',markersize=15)
    ax_no_inactivation_cost.plot(V,Cost_no_inactivation[i],colour_graph[i] + '.',markersize=15, markerfacecolor='None')
    ax_no_inactivation_combined.plot(Bandwidth[i],Cost[i],colour_graph[i] + '.',markersize=15)
    ax_no_inactivation_combined.plot(Bandwidth[i], Cost_no_inactivation[i], colour_graph[i] + '.', markersize=15, markerfacecolor='None')

#figure(6)
ax_Z.set_xlabel("Frequency (Hz)")
ax_Z.set_ylabel(r"Impedance (M$\Omega$)")
#ax_Z.set_ylim([10, 600])
ax_Z.legend(loc=3,prop={'size':12})


Cost_no_inactivation_a = zeros([len(Vr),len(Vr)])
Bandwidth_no_inactivation_a = zeros([len(Vr),len(Vr)])

for ii_no_inactivation in range(len(Vr)):
    for i,V in enumerate(Vr): #Highest and lowest
        DepolarisePhotoreceptor.WithLight(HH_no_inactivation[ii_no_inactivation],V)
        Cost_no_inactivation_a[i,ii_no_inactivation] = HH_no_inactivation[ii_no_inactivation].energy_consumption()*delta_cost
        _, Bandwidth_no_inactivation_a[i,ii_no_inactivation] = Gain_Bandwidth(HH_no_inactivation[ii_no_inactivation].body.impedance, f_min=f_medium)

for ii_no_inactivation in range(len(Vr)):
    ax_no_inactivation_cost.plot(Vr,Cost_no_inactivation_a[:,ii_no_inactivation],colour_graph[ii_no_inactivation], zorder=0)
    ax_no_inactivation_bw.plot(Vr,Bandwidth_no_inactivation_a[:,ii_no_inactivation],colour_graph[ii_no_inactivation], zorder=0)
    ax_no_inactivation_combined.plot(Bandwidth_no_inactivation_a[:,ii_no_inactivation],Cost_no_inactivation_a[:,ii_no_inactivation],colour_graph[ii_no_inactivation], zorder=0)

ax_no_inactivation_cost.plot(Vr,Cost,'k',zorder=0)
ax_no_inactivation_cost.set_xlabel("Membrane Voltage (mV)")
ax_no_inactivation_cost.set_ylabel(r"Energy Consumption ($10^8$ ATP/s)")
ax_no_inactivation_bw.plot(Vr,Bandwidth,'k',zorder=0)
ax_no_inactivation_bw.set_xlabel("Membrane Voltage (mV)")
ax_no_inactivation_bw.set_ylabel("Bandwidth (Hz)")
ax_no_inactivation_combined.plot(Bandwidth,Cost,'k',zorder=0)
ax_no_inactivation_combined.set_ylabel(r"Energy Consumption ($10^8$ ATP/s)")
ax_no_inactivation_combined.set_xlabel("Bandwidth (Hz)")


plt.show()

