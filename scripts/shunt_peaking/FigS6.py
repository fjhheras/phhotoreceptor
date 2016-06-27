#!/usr/bin/python3

from pylab import *
from numpy import *
import FlyFactory as FlyFactory
from phhotoreceptor.DepolarisePhotoreceptor import DepolarisePhotoreceptor
import phhotoreceptor.Experiment as Experiment
from GBWPutils import GBWP

__author__ = 'Francisco J. H. Heras'

HH  = FlyFactory.CalliphoraR16(channel_choice = "Weckstrom")
phasearray = vectorize (lambda z : angle(z))

### ONLY CERTAIN VOLTAGES BUT CONTINUOUS ACROSS FREQUENCIES

Vr=array([-60,-40])
delta_f = 0.1
f = arange(0.5,500,delta_f)

colour_graph=['b','r']

for i,V in enumerate(Vr):

    DepolarisePhotoreceptor.WithLight(HH,V)

    Z = HH.body.impedance(f) #All frequencies
    print("GBWP is ", GBWP(HH.body.impedance))

    Experiment.freeze_conductances(HH)
    Z_fixed = HH.body.impedance(f)
    print("Passive GBWP is ", GBWP(HH.body.impedance))
    Experiment.unfreeze_conductances(HH)

    figure(2) #phase lag
    plot(f,-phasearray(Z)/(2*pi*f),colour_graph[i])
    plot(f,-phasearray(Z_fixed)/(2*pi*f),colour_graph[i]+'--')
    xscale('log')

    figure(3) #group delay
    plot(f[1:],-convolve(phasearray(Z),[1,-1],'valid')/(2*pi*delta_f),colour_graph[i])
    plot(f[1:],-convolve(phasearray(Z_fixed),[1,-1],'valid')/(2*pi*delta_f),colour_graph[i]+'--') #resistance fitted to low frequencies + m.c.
    xscale('log')


figure(2)
xlabel("Frequency")
ylabel("Phase delay (s)")

figure(3)
xlabel("Frequency")
ylabel("Group delay (s)")


show()
