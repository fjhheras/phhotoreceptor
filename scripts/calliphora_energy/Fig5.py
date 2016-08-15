#!/usr/bin/python3
# Figure 1 -> Fig 5 in paper
# Figure 2 -> GBWP of active membranes (upper points), and the corresponding passive membranes (lower points)

from pylab import *
from numpy import *
import FlyFactory
from phhotoreceptor.DepolarisePhotoreceptor import DepolarisePhotoreceptor
import phhotoreceptor.Experiment as Experiment
from GBWPutils import GBWP, Gain_Bandwidth

option_debugging = False
HH  = FlyFactory.CalliphoraR16(channel_choice = "Anderson") #

####### BODY STARTS HERE

f_low = 0 #Hz
f_medium = 0 #1Hz

fig1 = figure(1,figsize=[9,5])
ax_Z = fig1.add_subplot(1,2,1)
ax_bw = fig1.add_subplot(2,2,4)
ax_Z_V_low = fig1.add_subplot(2,2,2)

fig2 = figure(2, figsize=[9,5])
ax_GBWP = fig2.add_subplot(1,1,1)
ax_GBWP.set_xlabel("Membrane potential (mV)")
ax_GBWP.set_ylabel("GBWP (MOhm Hz)")

###### CONTINUOUS ACROSS VOLTAGES
Vr = arange(-60.0,-35.0,0.5)
Z_array = zeros(len(Vr))
Z_fixed_array = zeros(len(Vr))
total_resistance_array = zeros_like(Vr)
total_K_conductance = zeros_like(Vr)
for i,V in enumerate(Vr): #Impedances at lowest frequency
    DepolarisePhotoreceptor.WithLight(HH,V)
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


#ax_Z_V_low.set_xlabel("Voltage (mV)")
ax_Z_V_low.set_ylabel("Impedance (MOhms)")
#ax_Z_V_low.set_yscale('log')
#ax_Z_V_low.set_ylim([30, 500])
ax_Z_V_low.set_xticklabels([])
ax_Z_V_low.yaxis.set_label_position("right")
ax_Z_V_low.yaxis.tick_right()
### ONLY CERTAIN VOLTAGES BUT CONTINUOUS ACROSS FREQUENCIES

Vr=array([-60,-52,-44,-37])
delta_f = 0.1
f = arange(1.5,900,delta_f)
f_from_medium = arange(f_medium,500,1)

colour_graph=['b','g','r','c']
Bandwidth = zeros_like(Vr)
Bandwidth_RC = zeros_like(Vr)

GBWP_ = zeros_like(Vr)
GBWP_RC_ = zeros_like(Vr)

for i,V in enumerate(Vr):

    DepolarisePhotoreceptor.WithLight(HH,V)
    Z = HH.body.impedance(f) #All frequencies
    GBWP_[i] = GBWP(HH.body.impedance, f_min =f_medium)
    pippo,Bandwidth[i] = Gain_Bandwidth(HH.body.impedance, f_min = f_medium)

    Experiment.freeze_conductances(HH)
    Z_fixed = HH.body.impedance(f)
    GBWP_RC_[i] = GBWP(HH.body.impedance, f_min =f_medium)
    pippo,Bandwidth_RC[i] = Gain_Bandwidth(HH.body.impedance, f_min = f_medium)
    Experiment.unfreeze_conductances(HH)

    label_str = str(V) + ' mV'
    ax_Z.loglog(f,abs(Z)/1000,colour_graph[i],linewidth=2,label = label_str)
    ax_Z.loglog(f,abs(Z_fixed)/1000,colour_graph[i]+'--',linewidth=2)

    ax_bw.plot(V,Bandwidth[i],colour_graph[i] + '.',markersize=15)
    ax_bw.plot(V,Bandwidth_RC[i],colour_graph[i] + '.',markersize=15)

    ax_GBWP.plot(V,GBWP_[i],colour_graph[i] + '.',markersize=15)
    ax_GBWP.plot(V,GBWP_RC_[i],colour_graph[i] + '.',markersize=15)

    ## Add bullet points in continuous graph across frequencies
    Z_low = abs(HH.body.impedance(f_low))
    Experiment.freeze_conductances(HH)
    Z_low_fixed = abs(HH.body.impedance(f_low))
    Experiment.unfreeze_conductances(HH)
    ax_Z_V_low.plot(V,Z_low/1000,colour_graph[i] + '.',markersize=15)
    ax_Z_V_low.plot(V,Z_low_fixed/1000,colour_graph[i] + '.',markersize=15)



ax_Z.set_xlabel("Frequency (Hz)")
ax_Z.set_ylabel("Impedance (MOhms)")
#ax_Z.set_ylim([10, 600])
ax_Z.legend(loc=1,prop={'size':12})
ax_bw.set_ylabel("Bandwidth (Hz)")
ax_bw.plot(Vr,Bandwidth,'k',zorder=0)
ax_bw.plot(Vr,Bandwidth_RC,'k--',zorder=0)
#ax_bw.set_xlim([-72, -28])
ax_bw.yaxis.set_label_position("right")
ax_bw.yaxis.tick_right()
ax_bw.yaxis.set_ticks_position('both')
#ax_bw.set_xticklabels([])


show()
