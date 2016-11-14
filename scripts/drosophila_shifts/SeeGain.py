#!/usr/bin/python3

from copy import deepcopy
from pylab import *
from numpy import *
import FlyFactory
from phhotoreceptor.DepolarisePhotoreceptor import DepolarisePhotoreceptor
import phhotoreceptor.Experiment as Experiment
from GBWPutils import Gain_Bandwidth


plot_freq_lims = [.1,300]

phasearray = vectorize (lambda z : angle(z))

HH  = FlyFactory.DrosophilaR16()

####### BODY STARTS HERE

f_low = 0 #Hz
f_medium = 2 #Hz
f_reference_gain = -1 #Hz If negative, max gain is taken

fig1 = figure(6,figsize=[9,5])
ax_RC_gain = fig1.add_subplot(121)
ax_bwprod = fig1.add_subplot(122)


### ONLY CERTAIN VOLTAGES BUT CONTINUOUS ACROSS FREQUENCIES

Vr=range(-68,-30,8)
V_rest = HH.V_rest()
delta_f = 0.1
f = arange(.2,200,delta_f)

colour_graph=['y','b','g','r','c']
Bandwidth = zeros_like(Vr)
Bandwidth_fixed = zeros_like(Vr)
gain_max = zeros_like(Vr)
gain_max_fixed = zeros_like(Vr)

gain_bandwidth_product = zeros_like(Vr)
gain_bandwidth_product_fixed = zeros_like(Vr)

HH_RC = []

for i,V in enumerate(Vr):

    label_str = str(V) + ' mV'

    DepolarisePhotoreceptor.WithLight(HH,V)

    gain_max[i],Bandwidth[i] = Gain_Bandwidth(HH.body.voltage_contrast_gain, f_min = f_medium)
    gain = abs(HH.body.voltage_contrast_gain(f))
    gain_bandwidth_product[i] = gain_max[i]*Bandwidth[i]

    Experiment.freeze_conductances(HH)

    gain_max_fixed[i],Bandwidth_fixed[i] = Gain_Bandwidth(HH.body.voltage_contrast_gain, f_min = f_medium)
    gain_fixed = abs(HH.body.voltage_contrast_gain(f))
    gain_bandwidth_product_fixed[i] = gain_max_fixed[i]*Bandwidth_fixed[i]

    Experiment.unfreeze_conductances(HH)

    ax_bwprod.plot(V,gain_bandwidth_product[i],colour_graph[i] + '.',markersize=15)
    ax_bwprod.plot(V,gain_bandwidth_product_fixed[i],colour_graph[i] + '.',markersize=15)

    if V>V_rest :
        ax_RC_gain.loglog(f,gain,colour_graph[i],linewidth=2,label = label_str)
        ax_RC_gain.loglog(f,gain_fixed,colour_graph[i]+':',linewidth=2,label = label_str)

ax_RC_gain.set_xlabel("Frequency (Hz)")
ax_RC_gain.set_ylabel("Gain (mV)")

ax_bwprod.set_xlabel("Frequency (Hz)")
ax_bwprod.set_ylabel("CGBWP (mV Hz)")
ax_bwprod.plot(Vr,gain_bandwidth_product,'k')
ax_bwprod.plot(Vr,gain_bandwidth_product_fixed,'k--')

show()
