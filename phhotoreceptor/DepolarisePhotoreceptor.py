from __future__ import division
from numpy import *


class DepolarisePhotoreceptor:
    """Different modes of depolarising a photoreceptor"""
    def WithLight (photoreceptor,V=None,g_light=None,verbose=0):
        photoreceptor.depolarise_with_light(V=V,g_light=g_light) #this is a member of the phhotoreceptor class
        if verbose >= 1:
            print ("Depolarising with light to V=",photoreceptor.body.V_m, " mV")
        if verbose >= 2:
            print("Pump current is ", photoreceptor.body.I_pump*1e3, " nA")
            print("Energy consumption is ", photoreceptor.energy_consumption(), " ATP/s")
            print("Light conductance is ", photoreceptor.body.light_conductance.g_max , " mS")
            print("Light leak conductance is ", photoreceptor.body.leak_conductances['L'].g_max, " mS")
    def WithCurrent (photoreceptor,V):
        photoreceptor.reset_voltage(V)

        
