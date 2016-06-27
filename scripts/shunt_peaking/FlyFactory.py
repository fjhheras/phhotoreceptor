from __future__ import division
from math import pi,sqrt

import CalliphoraConductance as CalliphoraConductance
from phhotoreceptor import FlyPhotoreceptor

__author__ = 'Francisco J. H. Heras'

def CalliphoraR16(channel_choice="Weckstrom", cascade_choice="delta"):
    reversal_potentials = {'L': 5, 'K' : -85}
    axon = FlyPhotoreceptor.Axon(l=35e-4, r=1e-4) #cm,cm
    # photon flux I to voltage V: Vm*(RS*I)**a/((RS*I)**phi + 1), Vm = 24.3, RS=1/2e4,f=0.42
    if channel_choice == "Weckstrom":
        #body.add_leak_conductance(ion_name="K",g=1e-16) #mS,mV
        lic_current = [0.0, 1.0, 0.0] #All LIC current is Na, so there is no exchanger
        body = FlyPhotoreceptor.CellBody(l=250e-4, r=2.5e-4, area=1.3e-4, reversal_potentials=reversal_potentials,lic_current=lic_current) #cm,cm,cm2
        #body.add_leak_conductance(ion_name="K",g=0.5e-6) #mS,mV
        body.add_voltage_channel(CalliphoraConductance.FastWeckstrom91(30e-6)) #mS
        body.add_voltage_channel(CalliphoraConductance.SlowWeckstrom91(30e-6)) #mS
    if channel_choice == "WeckstromFDR":
        #body.add_leak_conductance(ion_name="K",g=1e-16) #mS,mV
        lic_current = [0.0, 1.0, 0.0] #All LIC current is Na, so there is no exchanger
        body = FlyPhotoreceptor.CellBody(l=250e-4, r=2.5e-4, area=1.3e-4, reversal_potentials=reversal_potentials,lic_current=lic_current) #cm,cm,cm2
        #body.add_leak_conductance(ion_name="K",g=0.5e-6) #mS,mV
        body.add_voltage_channel(CalliphoraConductance.FastWeckstrom91(30e-6)) #mS

    if channel_choice == "passive":
        body = FlyPhotoreceptor.CellBody(l=250e-4, r=2.5e-4, area=1.45e-4, reversal_potentials=reversal_potentials)  # cm,cm,cm2

    return FlyPhotoreceptor.FlyPhotoreceptor(V_rest=-60, body=body, axon=axon)

def PassiveCalliphoraR16WithBandwidth (bandwidth,V,low_limit_frequency = 0):
    new_photoreceptor = CalliphoraR16(channel_choice="passive")
    C = new_photoreceptor.body.C
    R_same_bandwidth = 1/( 2*pi*C*sqrt(bandwidth**2-2*(low_limit_frequency)**2)  ) #Hz -> MOhm
    new_photoreceptor.set_steady_state(new_photoreceptor.V_rest(),VZ_0=[V, R_same_bandwidth*1000]) #MOhm -> KOhm
    return new_photoreceptor