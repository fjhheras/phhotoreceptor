#!/usr/bin/python3

## Plots a) Bandwidth vs (maximum) gain, b) GBWP vs cost
## For a Drosophila photoreceptor with unshifted conductances
## and for three different shifts

from pylab import *
from numpy import *
from matplotlib import rc
import FlyFactory
from phhotoreceptor.DepolarisePhotoreceptor import DepolarisePhotoreceptor
import phhotoreceptor.Experiment as Experiment
import ShiftConductances

#rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
rc('font',**{'family':'serif'}) #,'serif':['Liberation Serif']})
#rc('text', usetex=True)

option_debugging = False

phasearray = vectorize (lambda z : angle(z))
def Calculate_Bandwidth(Z,f):
    Z_max = abs(Z[0])
    for i,ff in enumerate(f):
        if abs(Z[i]) > Z_max:
            Z_max = abs(Z[i])
        if Z_max/sqrt(2) > abs(Z[i]):
            return ff

####### BODY STARTS HERE

f_medium = 2 #Hz
delta_cost = 1e-8 #If change, change also labels
change_LIC_to_keep_depolarisation = False

### ONLY CERTAIN VOLTAGES BUT CONTINUOUS ACROSS FREQUENCIES

Vcont=arange(-62.0,-35.0,0.2)
Vr=arange(-68.0,-30.0,8.0)
delta_f = 0.1
epsilon =1e-6
f_from_medium = arange(f_medium,500,delta_f)

colour_graph=['y','b','g','r','c']
#conditions=['Shab-50','Shab+50','Serotonin','PIP2']
#condition_line=['g','y','b','r','k--']
conditions=['Serotonin','PIP2']
condition_line=['b','r','k--']

N=len(conditions)

Bandwidth_r = zeros((len(conditions) + 2, len(Vr)))
Cost_r = zeros((len(conditions) + 2, len(Vr)))
gain_max_r = zeros((len(conditions) + 2, len(Vr)))
Vr_new_r = zeros((len(conditions) + 2, len(Vr)))

Bandwidth = zeros((len(conditions) + 2, len(Vcont)))
Cost = zeros((len(conditions) + 2, len(Vcont)))
gain_max = zeros((len(conditions) + 2, len(Vcont)))
Vr_new = zeros((len(conditions) + 2, len(Vcont)))

fig, axes = plt.subplots(1, 3, figsize=(9,4))
labels = '(a) (b) (c)'.split()
plt.subplots_adjust(wspace=0.4, left=0.1, right=0.95, bottom=0.1, top=0.95)
for ax, label in zip(axes, labels):
    ax.tick_params(direction='in', top=True, right=True)
    ax.text(-0.12, 1.02, label, transform=ax.transAxes,
                    fontsize=14, va='top', ha='right')



ax_bw_gain, ax_cost_bw, ax_cost_gain = axes

## Only the values in Vr to plot fat circles

for i,V in enumerate(Vr):
    label_str = str(V) + ' mV'
    photoreceptor = FlyFactory.DrosophilaR16(shift="none")
    DepolarisePhotoreceptor.WithLight(photoreceptor,V=V)

    Z_cut = photoreceptor.body.impedance(f_from_medium) #Only frequencies above f_medium, to calculate bandwidth
    gain_max_r[N, i] = max(abs(photoreceptor.body.voltage_contrast_gain(f_from_medium)))
    Bandwidth_r[N, i] = Calculate_Bandwidth(Z_cut, f_from_medium)
    Cost_r[N, i] = photoreceptor.energy_consumption()*delta_cost

    if gain_max_r[N, i]>epsilon:
        #ax_gain_vs_bw.plot(gain_max_r[N, i], Bandwidth_r[N, i], colour_graph[i] + '.', markersize=15)
        #ax_cost_vs_gbwp.plot(Cost_r[N, i], gain_max_r[N, i] * Bandwidth_r[N, i], colour_graph[i] + '.', markersize=15)
        ax_bw_gain.plot(gain_max_r[N, i], Bandwidth_r[N, i], colour_graph[i] + '.', markersize=15)
        ax_cost_bw.plot(Bandwidth_r[N, i], Cost_r[N, i], colour_graph[i] + '.', markersize=15)
        ax_cost_gain.plot(gain_max_r[N, i], Cost_r[N, i], colour_graph[i] + '.', markersize=15)




    for ii,condition in enumerate(conditions):
        photoreceptor = FlyFactory.DrosophilaR16(shift="none")
        DepolarisePhotoreceptor.WithLight(photoreceptor,V=V)
        if condition=='PIP2':
            ShiftConductances.WithLight(photoreceptor, change_LIC_to_keep_depolarisation = change_LIC_to_keep_depolarisation)
        elif condition=='Serotonin':
            ShiftConductances.WithSerotonin(photoreceptor, change_LIC_to_keep_depolarisation = change_LIC_to_keep_depolarisation)
        elif condition=='Shab-50':
            Experiment.modify_conductance(photoreceptor, "Shab", .5, change_LIC_to_keep_depolarisation = change_LIC_to_keep_depolarisation)
        elif condition=='Shab+50':
            Experiment.modify_conductance(photoreceptor, "Shab", 1.5, change_LIC_to_keep_depolarisation = change_LIC_to_keep_depolarisation)
        else:
            print("Warning, option not recognised")

        Vr_new_r[ii, i] = photoreceptor.body.V_m

        Z_cut = photoreceptor.body.impedance(f_from_medium) #Only frequencies above f_medium, to calculate bandwidth
        gain_max_r[ii, i] = max(abs(photoreceptor.body.voltage_contrast_gain(f_from_medium)))
        Bandwidth_r[ii, i] = Calculate_Bandwidth(Z_cut, f_from_medium)
        Cost_r[ii, i] = photoreceptor.energy_consumption()*delta_cost

        if gain_max_r[ii, i]>epsilon:
            #ax_gain_vs_bw.plot(gain_max_r[ii, i], Bandwidth_r[ii, i], colour_graph[i] + '.', markersize=15)
            #ax_cost_vs_gbwp.plot(Cost_r[ii, i], gain_max_r[ii, i] * Bandwidth_r[ii, i], colour_graph[i] + '.', markersize=15)
            ax_bw_gain.plot(gain_max_r[ii, i], Bandwidth_r[ii, i], colour_graph[i] + '.', markersize=15)
            ax_cost_bw.plot(Bandwidth_r[ii, i], Cost_r[ii, i], colour_graph[i] + '.', markersize=15)
            ax_cost_gain.plot(gain_max_r[ii, i], Cost_r[ii, i], colour_graph[i] + '.', markersize=15)



## Now calculating more fine array of frequencies Vcont

for i,V in enumerate(Vcont):
    photoreceptor = FlyFactory.DrosophilaR16(shift="none")
    DepolarisePhotoreceptor.WithLight(photoreceptor,V=V)

    Z_cut = photoreceptor.body.impedance(f_from_medium) #Only frequencies above f_medium, to calculate bandwidth
    gain_max[N, i] = max(abs(photoreceptor.body.voltage_contrast_gain(f_from_medium)))
    Bandwidth[N, i] = Calculate_Bandwidth(Z_cut, f_from_medium)
    Cost[N, i] = photoreceptor.energy_consumption()*delta_cost

    for ii,condition in enumerate(conditions):
        photoreceptor = FlyFactory.DrosophilaR16(shift="none")
        DepolarisePhotoreceptor.WithLight(photoreceptor,V=V)
        if condition=='PIP2':
            ShiftConductances.WithLight(photoreceptor, change_LIC_to_keep_depolarisation = change_LIC_to_keep_depolarisation)
        elif condition=='Serotonin':
            ShiftConductances.WithSerotonin(photoreceptor, change_LIC_to_keep_depolarisation = change_LIC_to_keep_depolarisation)
        elif condition=='Shab-50':
            Experiment.modify_conductance(photoreceptor, "Shab", .5, change_LIC_to_keep_depolarisation = change_LIC_to_keep_depolarisation)
        elif condition=='Shab+50':
            Experiment.modify_conductance(photoreceptor, "Shab", 1.5, change_LIC_to_keep_depolarisation = change_LIC_to_keep_depolarisation)
        else:
            print("Warning, option not recognised")

        Vr_new[ii, i] = photoreceptor.body.V_m

        Z_cut = photoreceptor.body.impedance(f_from_medium) #Only frequencies above f_medium, to calculate bandwidth
        gain_max[ii, i] = max(abs(photoreceptor.body.voltage_contrast_gain(f_from_medium)))
        Bandwidth[ii, i] = Calculate_Bandwidth(Z_cut, f_from_medium)
        Cost[ii, i] = photoreceptor.energy_consumption()*delta_cost


## Plotting condition N (no shift)

ax_bw_gain.plot(gain_max[N, 1:], Bandwidth[N, 1:], condition_line[N], zorder=0)
ax_cost_bw.plot(Bandwidth[N, 1:], Cost[N, 1:], condition_line[N], zorder=0)
ax_cost_gain.plot(gain_max[N, 1:], Cost[N, 1:], condition_line[N], zorder=0)


## Plotting other conditions

for ii,condition in enumerate(conditions):
    ax_bw_gain.plot(gain_max[ii, 1:], Bandwidth[ii, 1:], condition_line[ii], zorder=0)
    ax_cost_bw.plot(Bandwidth[ii, 1:], Cost[ii, 1:], condition_line[ii], zorder=0)
    ax_cost_gain.plot(gain_max[ii, 1:], Cost[ii, 1:], condition_line[ii], zorder=0)


ax_bw_gain.set_xlabel("Peak voltage contrast gain (mV)")
ax_cost_bw.set_xlabel("Bandwidth (Hz)")
ax_cost_gain.set_xlabel("Peak voltage contrast gain (mV)")
ax_bw_gain.set_ylabel("Bandwidth (Hz)")
ax_cost_bw.set_ylabel(r"Cost ($10^8$ ATP/s)")
ax_cost_gain.set_ylabel("Cost ($10^8$ ATP/s)")

show()
