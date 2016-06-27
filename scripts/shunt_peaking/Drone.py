from __future__ import division

from numpy import *
import phhotoreceptor.Conductance as Conductance
from phhotoreceptor import FlyPhotoreceptor
from math import exp

__author__ = 'Francisco J.H. Heras'

## Admittedly not a fly, but...
## Drone conductances

class HH_Na(Conductance.VoltageDependent): ## K Shab channel, aka the "slow delayed rectifier"
    m_order = 3
    def alpha_m(self,V):
        return 0.1*(-V + self.V_r + 25)/(exp((-V + self.V_r + 25)/10) - 1) *self.T_correction /self.m_time_multiplier
    def beta_m(self, V):
        return 4*exp(-(V-self.V_r)/18) *self.T_correction /self.m_time_multiplier

    inactivation_modes = 1
    def alpha_h (self,V):
        return array([ 0.07*exp(-(V-self.V_r)/20) ]) * self.k_h *self.T_correction /self.h_time_multiplier
    def beta_h (self,V):
        return array([ 1/(exp((-V + self.V_r + 30)/10) + 1) ]) * self.k_h *self.T_correction /self.h_time_multiplier

    def __init__(self,g_max,parent=None,k_h=0.1):
        Conductance.VoltageDependent.__init__(self,g_max,ion_name='Na', name = "Na",parent=parent)
        self.inactivation_fraction_mode = [1]
        self.h_order = array([1])
        self.h_time_multiplier = 1
        self.m_time_multiplier = 1
        self.k_h = k_h #factor to scale the time constant of inactivation (as in Vallet et al 1992)
        self.V_r = -55.5 #rest to include into the HH equations
        self.Temperature = 20 #Celsius
        self.T_correction = 1 #No correction used in Vallet et al
        #self.T_correction = power(3,0.1*(self.Temperature - 6.3)) #-> This is used in G.C. Taylor, J.A. Coles, J.C. Eilbeck (1995)

## Drone factory

def Vallet92(k_h = .1):
    reversal_potentials = {'L': 0, 'K' : -66, 'Na' : 57}
    area = 1.3e-4 #cm2, completely arbitrary
    axon = FlyPhotoreceptor.Axon(l=35e-4, r=1e-4) #cm,cm
    lic_current = [0.0, 1.0, 0.0] #All LIC current is Na, so there is no exchanger
    body = FlyPhotoreceptor.CellBody(l=250e-4, r=2.5e-4, area=area, reversal_potentials=reversal_potentials,lic_current=lic_current,is_there_pump=False) #cm,cm,cm2
    body.add_leak_conductance(ion_name="K",g=0.2*area) #0.2 mS cm-2
    body.add_voltage_channel(HH_Na(g_max=4*area, k_h=k_h)) #4.0 mS cm-2

    return FlyPhotoreceptor.FlyPhotoreceptor(V_rest=-55.5, body=body, axon=axon)