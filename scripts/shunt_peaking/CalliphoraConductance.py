from __future__ import division

from numpy import *
import phhotoreceptor.Conductance as Conductance
from math import exp

__author__ = 'Francisco J.H. Heras'


## Calliphora R1-6 conductance

class FastWeckstrom91(Conductance.DelayedRectifier):
    m_order = 1
    Vps = 15 #mV patch shift
    time_multiplier = 1
    def m_inf(self,V):
        return 1/(1+exp(-(V-self.Vps+65)/8.5))
    def m_time(self,V):
        return self.time_multiplier/( 3*exp((V-self.Vps+15)/24.4) + 9.4e-8*exp(-(V-self.Vps+15)/7.8) )

class SlowWeckstrom91(Conductance.DelayedRectifier):
    m_order = 1
    Vps = 15 #mV patch shift
    time_multiplier = 1
    def alpha_m(self,V):
        return 0.9*exp((V-self.Vps)/13)/self.time_multiplier
    def beta_m(self, V):
        return 0.0037*exp(-(V-self.Vps)/33.8)/self.time_multiplier
