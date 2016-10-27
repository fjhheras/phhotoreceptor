#!/usr/bin/python3

from pylab import *
from numpy import *
import FlyFactory
from phhotoreceptor.DepolarisePhotoreceptor import DepolarisePhotoreceptor
import phhotoreceptor.Experiment as Experiment
from GBWPutils import GBWP, Gain_Bandwidth

option_debugging = False
depolarise_with_light = True #If depolarise with current, all cost calculations are not biological

HH  = FlyFactory.CalliphoraR16(channel_choice = "Anderson") #

def calculate_bandwidth_of_passive_photoreceptor(photoreceptor,low_limit_frequency = 0):
    ## Not used, but can be used to debug
    C = photoreceptor.body.C
    R_a = photoreceptor.body.resistance()
    bw_when_low_freq_is_0 = 1000/( 2*pi*R_a*C) #Hz
    correction = 2*low_limit_frequency**2 #Hz2
    return sqrt(bw_when_low_freq_is_0**2 + correction)

####### BODY STARTS HERE

f_low = 0 #Hz
f_medium = 0 #1Hz

fig1 = figure(6,figsize=[9,5]) # Plot comparing against RC
ax_RC_bw = fig1.add_subplot(131)
ax_RC_cost = fig1.add_subplot(132)
ax_RC_combined = fig1.add_subplot(133)

fig2 = figure(5,figsize=[9,5])
ax_Z = fig2.add_subplot(1,2,1)
ax_cost = fig2.add_subplot(2,2,4)
ax_R_0 = fig2.add_subplot(2, 2, 2)


###### CONTINUOUS ACROSS VOLTAGES
Vr = arange(-60.0,-35.0,0.5)
Z_array = zeros(len(Vr))
Z_fixed_array = zeros(len(Vr))
total_resistance_array = zeros_like(Vr)
total_K_conductance = zeros_like(Vr)
for i,V in enumerate(Vr): #Impedances at lowest frequency
    if depolarise_with_light:
        DepolarisePhotoreceptor.WithLight(HH,V)
    else:
        DepolarisePhotoreceptor.WithCurrent(HH,V)
    Experiment.freeze_conductances(HH)
    Experiment.unfreeze_conductances(HH)
    total_resistance_array[i] = HH.body.resistance()
    total_K_conductance[i] = HH.body.total_K_conductance()

ax_R_0.plot(Vr, total_resistance_array / 1000, 'k--', linewidth=2, label ="Non-fixed conductances")


#ax_Z_V_low.set_xlabel("Voltage (mV)")
ax_R_0.set_ylabel("Impedance (MOhms)")
#ax_Z_V_low.set_yscale('log')
#ax_Z_V_low.set_ylim([30, 500])
ax_R_0.set_xticklabels([])
ax_R_0.yaxis.set_label_position("right")
ax_R_0.yaxis.tick_right()
### ONLY CERTAIN VOLTAGES BUT CONTINUOUS ACROSS FREQUENCIES

Vr=array([-60,-52,-44,-37])
delta_f = 0.1
f = arange(1.5,900,delta_f)
f_from_medium = arange(f_medium,500,1)

colour_graph=['b','g','r','c']
Bandwidth = zeros_like(Vr)
Cost = zeros_like(Vr)
Cost_RC = zeros_like(Vr)
HH_RC = []

for i,V in enumerate(Vr):

    DepolarisePhotoreceptor.WithLight(HH,V,verbose=2)
    C = HH.body.C
    Z = HH.body.impedance(f) #All frequencies
    Cost[i] = HH.energy_consumption()
    pippo, Bandwidth[i] = Gain_Bandwidth(HH.body.impedance, f_min=f_medium)
    HH_RC.append(FlyFactory.PassiveCalliphoraR16WithBandwidth(Bandwidth[i],V,low_limit_frequency=f_medium))
    DepolarisePhotoreceptor.WithLight(HH_RC[i],V)
    Cost_RC[i] = HH_RC[i].energy_consumption()
    Z_RC = HH_RC[i].body.impedance(f)#*HH.body.impedance(f_medium)/HH_RC[i].body.impedance(f_medium) #All frequencies


    label_str = str(V) + ' mV'
    ax_Z.loglog(f,abs(Z)/1000,colour_graph[i],linewidth=2,label = label_str)
    ax_Z.loglog(f,abs(Z_RC)/1000,colour_graph[i]+':',linewidth=1)

    ax_cost.plot(V,Cost[i],colour_graph[i] + '.',markersize=15)
    ax_cost.plot(V,Cost_RC[i],colour_graph[i] + '.',markersize=15, markerfacecolor='None')

    ## Add bullet points in continuous graph across frequencies
    R_0 = HH.body.resistance()
    R_passive_0 = HH_RC[i].body.resistance()
    ax_R_0.plot(V, R_0 / 1000, colour_graph[i] + '.', markersize=15)
    ax_R_0.plot(V, R_passive_0 / 1000, colour_graph[i] + '.', markersize=15, markerfacecolor='None')

    ax_RC_bw.plot(V,Bandwidth[i],colour_graph[i] + '.',markersize=15) # Cost and Bandwidth of the photoreceptors, then to compare with RC
    ax_RC_cost.plot(V,Cost[i],colour_graph[i] + '.',markersize=15)
    ax_RC_cost.plot(V,Cost_RC[i],colour_graph[i] + '.',markersize=15, markerfacecolor='None')
    ax_RC_combined.plot(Bandwidth[i],Cost[i],colour_graph[i] + '.',markersize=15)
    ax_RC_combined.plot(Bandwidth[i], Cost_RC[i], colour_graph[i] + '.', markersize=15, markerfacecolor='None')

#figure(5)
ax_Z.set_xlabel("Frequency (Hz)")
ax_Z.set_ylabel("Impedance (MOhms)")
#ax_Z.set_ylim([10, 600])
ax_Z.legend(loc=2,prop={'size':12})
ax_cost.set_xlabel("Membrane voltage (mV)")
ax_cost.set_ylabel("Cost (ATP/s)")
#ax_cost.set_ylim([5e7, 2e9])
#ax_cost.set_xlim([-72, -28])
ax_cost.plot(Vr,Cost,'k',zorder=0)
#ax_cost.plot(Vr,Cost_RC,'k:',zorder=0)
#ax_cost.set_yscale('log')
ax_cost.yaxis.set_label_position("right")
ax_cost.yaxis.tick_right()
ax_cost.yaxis.set_ticks_position('both')

# (fixed) RC membranes driven to different voltages. Figure(6)

Cost_RC_a = zeros([len(Vr),len(Vr)])
Bandwidth_RC_a = zeros([len(Vr),len(Vr)])

for ii_RC in range(len(Vr)): #Four different RC membranes
    for i,V in enumerate(Vr): #Highest and lowest
        DepolarisePhotoreceptor.WithLight(HH_RC[ii_RC],V)
        Cost_RC_a[i,ii_RC] = HH_RC[ii_RC].energy_consumption()
        #Bandwidth_RC_a[i,ii_RC] = calculate_bandwidth_of_passive_photoreceptor(HH_RC[ii_RC],low_limit_frequency=f_medium)
        pippo, Bandwidth_RC_a[i,ii_RC] = Gain_Bandwidth(HH_RC[ii_RC].body.impedance, f_min=f_medium)

for ii_RC in range(len(Vr)): #Four different RC membranes
    ax_RC_cost.plot(Vr,Cost_RC_a[:,ii_RC],colour_graph[ii_RC])
    ax_RC_bw.plot(Vr,Bandwidth_RC_a[:,ii_RC],colour_graph[ii_RC])
    ax_RC_combined.plot(Bandwidth_RC_a[:,ii_RC],Cost_RC_a[:,ii_RC],colour_graph[ii_RC])

ax_RC_cost.plot(Vr,Cost,'k',zorder=0)
ax_RC_cost.set_xlabel("Membrane Potential (mV)")
ax_RC_cost.set_ylabel("Cost (ATP/s)")
ax_RC_bw.plot(Vr,Bandwidth,'k',zorder=0)
ax_RC_bw.set_xlabel("Membrane Potential (mV)")
ax_RC_bw.set_ylabel("Bandwidth (Hz)")
ax_RC_combined.plot(Bandwidth,Cost,'k',zorder=0)
ax_RC_combined.set_ylabel("Cost (ATP/s)")
ax_RC_combined.set_xlabel("Bandwidth (Hz)")


show()
