#!/usr/bin/python3
## Option 1:
## Calculates the impedance, bandwith and cost of the membrane with light shifted
## Shab channel (solid line) and membrane with non light shifted Shab (dashed lines)
## Relevant reference is Krause et al. 2008, which shows that Shab suffers a 10 mV shift
## to more negative values with light caused by PIP2 decrease
## Option 2:
## Calculates the impedance, bandwith and cost of the membrane and the membrane shifted by
## 5-HT (serotonin) as explained in Hevers and Hardie 1995 (continuous line)
## Option 3:
## Decreases Shab conductance by 50%

from pylab import *
from numpy import *
import copy
import FlyFactory
from phhotoreceptor.DepolarisePhotoreceptor import DepolarisePhotoreceptor
import phhotoreceptor.Experiment as Experiment
import ShiftConductances
from GBWPutils import Gain_Bandwidth

option_debugging = False

####### BODY STARTS HERE

f_medium = 2 #Hz
option = 2
change_LIC_to_keep_depolarisation = False
photoreceptor = FlyFactory.DrosophilaR16(shift="none")

### ONLY CERTAIN VOLTAGES BUT CONTINUOUS ACROSS FREQUENCIES

Vr=arange(-68.0,-30.0,8.0)
V_rest = photoreceptor.V_rest()
delta_f = 0.2
f = arange(.5,200,delta_f)
f_from_medium = arange(f_medium,200,1)

fig1 = figure(1)
ax_Z = fig1.add_subplot(1,2,1)
ax_bw = fig1.add_subplot(2,2,2)
ax_cost = fig1.add_subplot(2,2,4)

fig2 = figure(2)
ax_bw_cost = fig2.add_subplot(1,2,1)
ax_gain_cost = fig2.add_subplot(1,2,2)

fig3 = figure(3)
ax_gain = fig3.add_subplot(1,2,1)
ax_gain_vs_bw = fig3.add_subplot(1,2,2)

fig6 = figure(6,figsize=[9,5])
ax_bwprod = fig6.add_subplot(121)

colour_graph=['y','b','g','r','c']
Bandwidth = zeros_like(Vr)
Bandwidth_shift = zeros_like(Vr)
Cost = zeros_like(Vr)
Cost_shift = zeros_like(Vr)
gain_max = zeros_like(Vr)
gain_max_shift = zeros_like(Vr)
gain_f= zeros_like(Vr)
gain_f_shift= zeros_like(Vr)
Vr_new = zeros_like(Vr)

for i,V in enumerate(Vr):
    label_str = str(V) + ' mV'
    DepolarisePhotoreceptor.WithLight(photoreceptor,V=V)
    Z = photoreceptor.body.impedance(f) #All frequencies

    ax_Z.loglog(f,abs(Z)/1000,colour_graph[i]+'--',linewidth=2,label = label_str)
    pippo,Bandwidth[i] = Gain_Bandwidth(photoreceptor.body.impedance,f_min = f_medium)
    Cost[i] = photoreceptor.energy_consumption()

    ax_bw.plot(V,Bandwidth[i],colour_graph[i] + '.',markersize=15, alpha=0.5)
    ax_cost.plot(V,Cost[i],colour_graph[i] + '.',markersize=15, alpha=0.5)
    ax_bw_cost.plot(Bandwidth[i], Cost[i] ,colour_graph[i] + '.',markersize=15, alpha=0.5)

    if V > V_rest :
        gain = abs(photoreceptor.body.voltage_contrast_gain(f))
        gain_max[i] = max(abs(photoreceptor.body.voltage_contrast_gain(f_from_medium)))
        ax_gain.loglog(f,abs(gain),colour_graph[i]+'--',linewidth=2,label = label_str)
        ax_bwprod.plot(Vr[i],gain_max[i]*Bandwidth[i],colour_graph[i] + '.',markersize=15)
        ax_gain_cost.plot(gain_max[i],Cost[i],colour_graph[i] + '.',markersize=15, alpha=0.5)
        ax_gain_vs_bw.plot(gain_max[i],Bandwidth[i],colour_graph[i] + '.',markersize=15, alpha=0.5)

    photoreceptor_shifted = copy.deepcopy(photoreceptor)
    if option==1:
        ShiftConductances.WithLight(photoreceptor_shifted, change_LIC_to_keep_depolarisation = change_LIC_to_keep_depolarisation)
    elif option==2:
        ShiftConductances.WithSerotonin(photoreceptor_shifted, change_LIC_to_keep_depolarisation = change_LIC_to_keep_depolarisation)
    else:
        Experiment.modify_conductance(photoreceptor_shifted, "Shab", .5, change_LIC_to_keep_depolarisation = change_LIC_to_keep_depolarisation)

    Vr_new[i] = photoreceptor_shifted.body.V_m

    Z = photoreceptor_shifted.body.impedance(f) #All frequencies
    Z_cut = photoreceptor_shifted.body.impedance(f_from_medium) #Only frequencies above f_medium, to calculate bandwidth

    gain = abs(photoreceptor_shifted.body.voltage_contrast_gain(f))
    gain_max_shift[i] = max(abs(photoreceptor_shifted.body.voltage_contrast_gain(f_from_medium)))

    ax_Z.loglog(f,abs(Z)/1000,colour_graph[i],linewidth=2)
    pippo,Bandwidth_shift[i] = Gain_Bandwidth(photoreceptor_shifted.body.impedance,f_min = f_medium)
    Cost_shift[i] = photoreceptor_shifted.energy_consumption()
    ax_bw.plot(Vr_new[i],Bandwidth_shift[i],colour_graph[i] + '.',markersize=15)
    ax_cost.plot(Vr_new[i],Cost_shift[i],colour_graph[i] + '.',markersize=15)
    ax_bw_cost.plot(Bandwidth_shift[i], Cost_shift[i],colour_graph[i] + '.',markersize=15)

    if V > V_rest :
        ax_gain.loglog(f,abs(gain),colour_graph[i],linewidth=2)
        ax_bwprod.plot(Vr_new[i],gain_max_shift[i]*Bandwidth_shift[i],colour_graph[i] + '.',markersize=15)
        ax_gain_cost.plot(gain_max_shift[i], Cost_shift[i],colour_graph[i] + '.',markersize=15, alpha=0.5)
        ax_gain_vs_bw.plot(gain_max_shift[i], Bandwidth_shift[i],colour_graph[i] + '.',markersize=15, alpha=0.5)


#figure(1)
ax_Z.set_xlabel("Frequency (Hz)")
ax_Z.set_ylabel("Impedance (MOhms)")
ax_Z.legend(loc=3,prop={'size':12})
ax_Z.set_ylim([10,400])
ax_Z.set_xlim([0.1,300])
ax_bw.set_ylabel("Bandwidth (Hz)")
ax_bw.plot(Vr,Bandwidth,'k--',zorder=0)
ax_bw.plot(Vr_new,Bandwidth_shift,'k',zorder=0)
ax_bw.set_xlim([-72, -30])
ax_bw.yaxis.set_label_position("right")
ax_bw.yaxis.tick_right()
ax_bw.yaxis.set_ticks_position('both')
ax_cost.set_xlabel("Membrane voltage (mV)")
ax_cost.set_ylabel("Cost (ATP/s)")
#ax_cost.set_ylim([5e7, 1e9])
ax_cost.set_xlim([-72, -30])
ax_cost.plot(Vr,Cost,'k--',zorder=0)
ax_cost.plot(Vr_new,Cost_shift,'k',zorder=0)
#ax_cost.set_yscale('log')
ax_cost.yaxis.set_label_position("right")
ax_cost.yaxis.tick_right()
ax_cost.yaxis.set_ticks_position('both')

#figure(2)
ax_bw_cost.set_xlabel("Bandwidth (Hz)")
ax_bw_cost.set_ylabel("Cost (ATP/s)")
#ax_bw_cost.set_yscale('log')
ax_bw_cost.plot(Bandwidth,Cost,'k--',zorder=0)
ax_bw_cost.plot(Bandwidth_shift, Cost_shift,'k',zorder=0)
ax_gain_cost.set_xlabel("Gain (mV)")
ax_gain_cost.set_ylabel("Cost (ATP/s)")
#ax_gain_cost.set_yscale('log')
ax_gain_cost.plot(gain_max[1:],Cost[1:],'k--',zorder=0)
ax_gain_cost.plot(gain_max_shift[1:], Cost_shift[1:],'k',zorder=0)

#figure(3)
ax_gain.set_xlabel("Frequency (Hz)")
ax_gain.set_ylabel("Contrast gain (mV)")
ax_gain.set_xlim([.1,300])
ax_gain.legend(loc=3,prop={'size':12})

ax_gain_vs_bw.plot(gain_max[1:],Bandwidth[1:],'k--',zorder=0)
ax_gain_vs_bw.plot(gain_max_shift[1:], Bandwidth_shift[1:],'k',zorder=0)
ax_gain_vs_bw.set_ylabel("Bandwidth (Hz)")
ax_gain_vs_bw.set_xlabel("Gain (mV)")

ax_bwprod.set_xlabel("Membrane voltage (mV)")
ax_bwprod.set_ylabel("Gain-bandwidth product (mV*Hz)")
ax_bwprod.plot(Vr,gain_max*Bandwidth,'k--')
ax_bwprod.plot(Vr_new,gain_max_shift*Bandwidth_shift,'k')

#ax_bw_cost.yaxis.set_label_position("right")
#ax_bw_cost.yaxis.tick_right()
#ax_bw_cost.yaxis.set_ticks_position('both')


if option == 1:
    print("Continuous line is Shab shifted by light")
    fig1.suptitle("Effect of light dependent shift in Shab")
elif option == 2:
    print("Continuous line is Shab and Shaker shifted by serotonin")
    fig1.suptitle("Effect of shift in channel properties by serotonin")
else:
    print("Continuous line is Shab current decreased 50%")
    fig1.suptitle("Effect of 50% decrease in Shab conductance by changes in calmodulin")
show()
