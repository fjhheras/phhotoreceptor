#!/usr/bin/python3

from pylab import *
from numpy import *
import FlyFactory
from phhotoreceptor.DepolarisePhotoreceptor import DepolarisePhotoreceptor
import phhotoreceptor.Experiment as Experiment
from GBWPutils import Gain_Bandwidth

option_debugging = False
depolarise_with_light = True   #If depolarise with current, all cost calculations are not biological

HH = FlyFactory.DrosophilaR16()

####### BODY STARTS HERE

f_low = 0  #Hz
f_medium = 2  #Hz

fig2 = figure(5, figsize=[9, 5])
ax_Z = fig2.add_subplot(1, 3, 1)
ax_bw = fig2.add_subplot(2, 3, 3)
ax_cost = fig2.add_subplot(2, 3, 6)
ax_Z_V_low = fig2.add_subplot(2, 3, 2)
ax_Z_V_medium = fig2.add_subplot(2, 3, 5)

###### CONTINUOUS ACROSS VOLTAGES
Vr = arange(-70.0,-30.0,0.5)
Z_array = zeros(len(Vr))
Z_fixed_array = zeros(len(Vr))
total_resistance_array = zeros_like(Vr)
total_K_conductance = zeros_like(Vr)
for i,V in enumerate(Vr): #Impedances at lowest frequency
    if depolarise_with_light:
        DepolarisePhotoreceptor.WithLight(HH,V)
    else:
        DepolarisePhotoreceptor.WithCurrent(HH,V)
    Z_array[i] = abs(HH.body.impedance(f_low))
    Experiment.freeze_conductances(HH)
    Z_fixed_array[i] = abs(HH.body.impedance(f_low))
    Experiment.unfreeze_conductances(HH)
    total_resistance_array[i] = HH.body.resistance()
    total_K_conductance[i] = HH.body.total_K_conductance()

ax_Z_V_low.plot(Vr,Z_array/1000,'k',linewidth=2,label = "Fixed conductances")
ax_Z_V_low.plot(Vr,Z_fixed_array/1000,'k--',linewidth=2, label = "Non-fixed conductances")
if option_debugging:
    ax_Z_V_low.plot(Vr,total_resistance_array/1000,'k:',linewidth=2) #total resistance,should coincide with Z_fixed
    ax_Z_V_low.plot(Vr,1/total_K_conductance/1000,'r:',linewidth=2) #K channel mediated part of the resistance


for i,V in enumerate(Vr): #Impedances at a medium frequency (low but above the amplification)
    if depolarise_with_light:
        DepolarisePhotoreceptor.WithLight(HH,V)
    else:
        DepolarisePhotoreceptor.WithCurrent(HH,V)
    Z_array[i] = abs(HH.body.impedance(f_medium))
    Experiment.freeze_conductances(HH)
    Z_fixed_array[i] = abs(HH.body.impedance(f_medium))
    Experiment.unfreeze_conductances(HH)

ax_Z_V_medium.plot(Vr,Z_array/1000,'k',linewidth=2,label = "Fixed conductances")
ax_Z_V_medium.plot(Vr,Z_fixed_array/1000,'k--',linewidth=2, label = "Non-fixed conductances")

#ax_Z_V_low.set_xlabel("Voltage (mV)")
ax_Z_V_low.set_ylabel("Impedance (MOhms)")
ax_Z_V_low.set_yscale('log')
ax_Z_V_low.set_ylim([30, 500])
ax_Z_V_low.set_xticklabels([])
ax_Z_V_medium.set_xlabel("Membrane voltage (mV)")
ax_Z_V_medium.set_ylabel("Impedance (MOhms)")
ax_Z_V_medium.set_yscale('log')
ax_Z_V_medium.set_ylim([30, 500])



### ONLY CERTAIN VOLTAGES BUT CONTINUOUS ACROSS FREQUENCIES


Vr=range(-68,-30,8)
delta_f = 0.1
f = arange(.2,200,delta_f)
f_from_medium = arange(f_medium,200,.01)

colour_graph=['y','b','g','r','c']
Bandwidth = zeros_like(Vr)
Bandwidth_fixed = zeros_like(Vr)
Cost = zeros_like(Vr)
Cost_RC = zeros_like(Vr)
HH_RC = []

for i,V in enumerate(Vr):

    if depolarise_with_light:
        DepolarisePhotoreceptor.WithLight(HH,V,verbose=2)
    else:
        DepolarisePhotoreceptor.WithCurrent(HH,V)
    C = HH.body.C
    Z = HH.body.impedance(f) #All frequencies
    pippo,Bandwidth[i] = Gain_Bandwidth(HH.body.impedance,f_min = f_medium)
    Cost[i] = HH.energy_consumption()
    print ("Cost is ", Cost[i], "ATP/s")


    Experiment.freeze_conductances(HH)
    Z_fixed = HH.body.impedance(f)
    pippo, Bandwidth_fixed[i] = Gain_Bandwidth(HH.body.impedance, f_min = f_medium)
    Experiment.unfreeze_conductances(HH)


    HH_RC.append(FlyFactory.PassiveDrosophilaR16WithBandwidth(Bandwidth[i],V,low_limit_frequency=f_medium))
    if depolarise_with_light:
        DepolarisePhotoreceptor.WithLight(HH_RC[i],V)
    else:
        DepolarisePhotoreceptor.WithCurrent(HH_RC[i],V)
    Cost_RC[i] = HH_RC[i].energy_consumption()
    Z_RC = HH_RC[i].body.impedance(f)


    label_str = str(V) + ' mV'
    ax_Z.loglog(f,abs(Z)/1000,colour_graph[i],linewidth=2,label = label_str)
    ax_Z.loglog(f,abs(Z_fixed)/1000,colour_graph[i]+'--',linewidth=2)
    #ax_Z.loglog(f,abs(Z_RC)/1000,colour_graph[i]+':',linewidth=1)

    ax_bw.plot(V,Bandwidth[i],colour_graph[i] + '.',markersize=15)
    ax_bw.plot(V, Bandwidth_fixed[i], colour_graph[i] + '.', markersize=15)
    ax_cost.plot(V,Cost[i],colour_graph[i] + '.',markersize=15)
    ax_cost.plot(V,Cost_RC[i],colour_graph[i] + '.',markersize=15)


    ## Add bullet points in continuous graph across frequencies
    Z_low = abs(HH.body.impedance(f_low))
    Z_medium = abs(HH.body.impedance(f_medium))
    Experiment.freeze_conductances(HH)
    Z_low_fixed = abs(HH.body.impedance(f_low))
    Z_medium_fixed = abs(HH.body.impedance(f_medium))
    Experiment.unfreeze_conductances(HH)
    ax_Z_V_low.plot(V,Z_low/1000,colour_graph[i] + '.',markersize=15)
    ax_Z_V_low.plot(V,Z_low_fixed/1000,colour_graph[i] + '.',markersize=15)
    ax_Z_V_medium.plot(V,Z_medium/1000,colour_graph[i] + '.',markersize=15)
    ax_Z_V_medium.plot(V,Z_medium_fixed/1000,colour_graph[i] + '.',markersize=15)


#figure(5)
ax_Z.set_xlabel("Frequency (Hz)")
ax_Z.set_ylabel("Impedance (MOhms)")
ax_Z.set_ylim([10, 400])
ax_Z.set_xlim([0.1,300])
ax_Z.legend(loc=3,prop={'size':12})
ax_bw.set_ylabel("Bandwidth (Hz)")
ax_bw.plot(Vr,Bandwidth,'k',zorder=0)
ax_bw.plot(Vr, Bandwidth_fixed, 'k--', zorder=0)
ax_bw.set_xlim([-72, -28])
ax_bw.yaxis.set_label_position("right")
ax_bw.yaxis.tick_right()
ax_bw.yaxis.set_ticks_position('both')
ax_bw.set_xticklabels([])


ax_cost.set_xlabel("Membrane voltage (mV)")
ax_cost.set_ylabel("Cost (ATP/s)")
#ax_cost.set_ylim([5e7, 2e9])
ax_cost.set_xlim([-72, -28])
ax_cost.plot(Vr,Cost,'k',zorder=0)
#ax_cost.plot(Vr,Cost_RC,'k:',zorder=0)
#ax_cost.set_yscale('log')
ax_cost.yaxis.set_label_position("right")
ax_cost.yaxis.tick_right()
ax_cost.yaxis.set_ticks_position('both')


show()
