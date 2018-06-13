#!/usr/bin/python3

## Plots two figures
## First:
## a) Bandwidth vs (maximum) gain, b) GBWP vs cost
## For a Drosophila photoreceptor with unshifted conductances
## and for three different shifts
## Second: Energy efficiency

import numpy as np
import matplotlib.pyplot as plt
import FlyFactory
from phhotoreceptor.DepolarisePhotoreceptor import DepolarisePhotoreceptor
import phhotoreceptor.Experiment as Experiment
import ShiftConductances

option_debugging = False

def Calculate_Bandwidth(Z, f):
    Z_max = abs(Z[0])
    for i, ff in enumerate(f):
        if abs(Z[i]) > Z_max:
            Z_max = abs(Z[i])
        if Z_max/np.sqrt(2) > abs(Z[i]):
            return ff

####### BODY STARTS HERE

f_medium = 2 # Hz
change_LIC_to_keep_depolarisation = False

### ONLY CERTAIN VOLTAGES BUT CONTINUOUS ACROSS FREQUENCIES

Vcont = np.arange(-62.0, -35.0, 0.2)
Vr = np.arange(-68.0, -30.0, 8.0)
delta_f = 0.1
epsilon = 1e-6
f_from_medium = np.arange(f_medium, 500, delta_f)

colour_graph = list('ybgrc')
conditions = 'Shab-50 Shab+50 Serotonin PIP2'.split()
condition_line = 'g y b r k--'.split()
N = len(conditions)

Bandwidth_r = np.zeros((len(conditions) + 2, len(Vr)))
Cost_r = np.zeros((len(conditions) + 2, len(Vr)))
gain_max_r = np.zeros((len(conditions) + 2, len(Vr)))
Vr_new_r = np.zeros((len(conditions) + 2, len(Vr)))

Bandwidth = np.zeros((len(conditions) + 2, len(Vcont)))
Cost = np.zeros((len(conditions) + 2, len(Vcont)))
gain_max = np.zeros((len(conditions) + 2, len(Vcont)))
Vr_new = np.zeros((len(conditions) + 2, len(Vcont)))

fig1 = plt.figure(1)
ax_gain_vs_bw = fig1.add_subplot(2,1,1)
ax_cost_vs_gbwp = fig1.add_subplot(2,1,2)

fig2 = plt.figure(2)
ax_eff = fig2.add_subplot(1,1,1)

## Only the values in Vr to plot fat circles

for i,V in enumerate(Vr):
    label_str = str(V) + ' mV'
    photoreceptor = FlyFactory.DrosophilaR16(shift="none")
    DepolarisePhotoreceptor.WithLight(photoreceptor,V=V)

    Z_cut = photoreceptor.body.impedance(f_from_medium) #Only frequencies above f_medium, to calculate bandwidth
    gain_max_r[N, i] = max(abs(photoreceptor.body.voltage_contrast_gain(f_from_medium)))
    Bandwidth_r[N, i] = Calculate_Bandwidth(Z_cut, f_from_medium)
    Cost_r[N, i] = photoreceptor.energy_consumption()

    if gain_max_r[N, i]>epsilon:
        ax_gain_vs_bw.plot(gain_max_r[N, i], Bandwidth_r[N, i], colour_graph[i] + '.', markersize=15)
        ax_cost_vs_gbwp.plot(Cost_r[N, i], gain_max_r[N, i] * Bandwidth_r[N, i], colour_graph[i] + '.', markersize=15)

    for ii,condition in enumerate(conditions):
        photoreceptor = FlyFactory.DrosophilaR16(shift="none")
        DepolarisePhotoreceptor.WithLight(photoreceptor,V=V)
        if condition == 'PIP2':
            ShiftConductances.WithLight(photoreceptor, change_LIC_to_keep_depolarisation = change_LIC_to_keep_depolarisation)
        elif condition == 'Serotonin':
            ShiftConductances.WithSerotonin(photoreceptor, change_LIC_to_keep_depolarisation = change_LIC_to_keep_depolarisation)
        elif condition == 'Shab-50':
            Experiment.modify_conductance(photoreceptor, "Shab", .5, change_LIC_to_keep_depolarisation = change_LIC_to_keep_depolarisation)
        elif condition == 'Shab+50':
            Experiment.modify_conductance(photoreceptor, "Shab", 1.5, change_LIC_to_keep_depolarisation = change_LIC_to_keep_depolarisation)
        else:
            print("Warning, option not recognised")

        Vr_new_r[ii, i] = photoreceptor.body.V_m

        Z_cut = photoreceptor.body.impedance(f_from_medium) #Only frequencies above f_medium, to calculate bandwidth
        gain_max_r[ii, i] = max(abs(photoreceptor.body.voltage_contrast_gain(f_from_medium)))
        Bandwidth_r[ii, i] = Calculate_Bandwidth(Z_cut, f_from_medium)
        Cost_r[ii, i] = photoreceptor.energy_consumption()

        if gain_max_r[ii, i]>epsilon:
            ax_gain_vs_bw.plot(gain_max_r[ii, i], Bandwidth_r[ii, i], colour_graph[i] + '.', markersize=15)
            ax_cost_vs_gbwp.plot(Cost_r[ii, i], gain_max_r[ii, i] * Bandwidth_r[ii, i], colour_graph[i] + '.', markersize=15)

## Now calculating more fine array of frequencies Vcont

for i,V in enumerate(Vcont):
    photoreceptor = FlyFactory.DrosophilaR16(shift="none")
    DepolarisePhotoreceptor.WithLight(photoreceptor,V=V)

    Z_cut = photoreceptor.body.impedance(f_from_medium) #Only frequencies above f_medium, to calculate bandwidth
    gain_max[N, i] = max(abs(photoreceptor.body.voltage_contrast_gain(f_from_medium)))
    Bandwidth[N, i] = Calculate_Bandwidth(Z_cut, f_from_medium)
    Cost[N, i] = photoreceptor.energy_consumption()

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
        Cost[ii, i] = photoreceptor.energy_consumption()


## Plotting condition N (no shift)

ax_gain_vs_bw.plot(gain_max[N, 1:], Bandwidth[N, 1:], condition_line[N], zorder=0)
ax_cost_vs_gbwp.plot(Cost[N, 1:], gain_max[N, 1:] * Bandwidth[N, 1:], condition_line[N], zorder=0)
ax_eff.plot(Cost[N, 1:], gain_max[N, 1:] * Bandwidth[N, 1:]/Cost[N, 1:], condition_line[N], zorder=0)


## Plotting other conditions

for ii,condition in enumerate(conditions):
    ax_gain_vs_bw.plot(gain_max[ii, 1:], Bandwidth[ii, 1:], condition_line[ii], zorder=0)
    ax_cost_vs_gbwp.plot(Cost[ii, 1:], gain_max[ii, 1:] * Bandwidth[ii, 1:], condition_line[ii], zorder=0)
    ax_eff.plot(Cost[ii, 1:], gain_max[ii, 1:] * Bandwidth[ii, 1:] / Cost[ii, 1:], condition_line[ii], zorder=0)

ax_gain_vs_bw.set_ylabel("Bandwidth (Hz)")
ax_gain_vs_bw.set_xlabel("Gain (mV)")
ax_cost_vs_gbwp.set_ylabel("cGBWP (mV Hz)")
ax_cost_vs_gbwp.set_xlabel("Cost (ATP/s)")

ax_eff.set_xlabel("Cost (ATP/s)")
ax_eff.set_ylabel("cGBWP efficiency (mV / ATP)")

plt.show()
