from __future__ import division

from numpy import *
import phhotoreceptor.Conductance as Conductance
from math import exp,sqrt

__author__ = 'Francisco J.H. Heras'

##Global functions

bell_shaped = vectorize( lambda x,A,a,b,c,d : A/(1+b*(x-a)**2 + d*(x-c)**4))

## Drosophila types below

class ShabNiven03(Conductance.VoltageDependent): ## K Shab channel, aka the "slow delayed rectifier"
    g_single_channel = 32*1e-9 #30-35 pS -> mS ; from Hardie(91)
    m_order = 2
    def alpha_m(self,V):
        return (  1/(1+exp(-0.12102084*V-0.63052366))  )**(1/2)/bell_shaped(V,5.61542061e+00,-2.17412676e+01,1.02208600e-03,-2.84855405e+01,3.37269152e-15)
    def beta_m(self, V):
        return (1-  (  1/(1+exp(-0.12102084*V-0.63052366))  )**(1/2))/bell_shaped(V,5.61542061e+00,-2.17412676e+01,1.02208600e-03,-2.84855405e+01,3.37269152e-15)

    inactivation_modes = 1
    def alpha_h (self,V):
        return array([ (1/(1+exp(0.15751498*V+4.06180895)))/890 ])
    def beta_h (self,V):
        return array([ (1-  1/(1+exp(0.15751498*V+4.06180895))  )/890 ])
    def __init__(self,g_max,parent=None): ## It has to be like this to avoid passing original array by reference
        Conductance.VoltageDependent.__init__(self,g_max,name = "Shab",parent=parent)
        self.inactivation_fraction_mode = [1]
        self.h_order = array([1])

#Weckstrom

#class ParameterSubunitWeckstromTemplate(modes,):
#
#
#class ParameterWecstromTemplate()
#
#class WeckstromTemplate(Conductance.VoltageDependent):
#    def __init__(self,parameter,g_single_channel):
#        self.weckstrom_template_parameters = parameter


class ShabWeckstrom(Conductance.VoltageDependent): ## K Shab channel, aka the "slow delayed rectifier", from Weckstrom model
    g_single_channel = 32*1e-9 #30-35 pS -> mS ; from Hardie(91)
    m_order = 2
    def m_time(self,V):
        return 1/( 0.116258*exp((-V-25.6551)/32.1933) + 0.00659219*(-V-23.8032)/(exp((-V-23.8032)/1.34548)-1) )
    def m_inf(self, V):
        return 1/( 1+exp((-1-V)/9.1) )**(1/2)

    inactivation_modes = 2
#    def set_inactivation_fractions(self,V): # but only affect above -13mV, so not used elsewhere
#        w = min(max(.72-.0134*V,.1),.9);
#        self.inactivation_fraction_mode = [w, 1-w]
    def h_inf(self,V):
        return array([ 1/( 1+exp((-25.7-V)/-6.4) ),    # Fast inactivation
                       1/( 1+exp((-25.7-V)/-6.4) )  ]) # Slow inactivation (but same h_inf)
    def h_time(self,V):
        return array([ 335*exp(V/71.3)+73.2, #Fast inactivation
                      3000])                 #Slow inactivation
    def __init__(self,g_max,inactivation_fraction = 1,parent=None): ## It has to be like this to avoid passing original array by reference
        Conductance.VoltageDependent.__init__(self,g_max, name = "Shab",parent=parent)
        self.h_order = array([1,1])
        self.inactivation_fraction_mode = [0.7, 0.3]*inactivation_fraction # it will depend on the voltage --> In FUTURE WILL UPDATE USING set_inactivation_fractions(V)


class ShabWeckstrom_LightShifted(ShabWeckstrom): #From Krause et al 2008: (I shift both 10mV, and they say it's only the slow component)
    def m_inf(self, V):
        return ShabWeckstrom.m_inf(self,V+10)
    #def m_time(self, V):
    #    return ShabWeckstrom.m_time(self,V+10)


class ShabWeckstrom_serotonin(ShabWeckstrom): #From Hevers and Hardie 1995
    def m_inf(self, V):
        return 1/( 1+exp((12.7-V)/9.5) )**(1/2)
    def h_inf(self,V):
        return array([ 1/( 1+exp((-15.7-V)/-7.2) ),    # Fast inactivation
                       1/( 1+exp((-15.7-V)/-7.2) )  ]) # Slow inactivation (but same h_inf)
#    def m_inf(self,V):
#        return ShabWeckstrom.m_inf(self,V-13.7)
    def m_time(self,V):
        return ShabWeckstrom.m_time(self,V-13.7)
#    def h_inf(self,V):
#        return ShabWeckstrom.h_inf(self,V-10)
    def h_time(self,V):
        return ShabWeckstrom.h_time(self,V-10)


class ShakerWeckstrom(Conductance.VoltageDependent): ## K Shaker channel, A-type, from Weckstrom model
    g_single_channel = 11*1e-9 #11 pS -> mS for outward currents, 19pS for Inward ; from Hardie(91)
    m_order = 3
    def m_time(self,V):
        return 1/( 0.008174*exp((-V+1.61882)/24.6538)  +  0.058139*(-V-59.639)/(exp((-V-59.639)/4.50122)-1) )
    def m_inf(self, V):
        return 1/(  1 + exp((-23.7-V)/12.8)  )**(1/3)
    inactivation_modes = 1
    def h_time(self,V):
        return array([ 1/( 0.230299*exp((-V-192.973)/31.31961) + 0.0437316*(-V+13.4859)/(exp((-V+13.4859)/11.11)-1) )  ])
    def h_inf(self,V):
        return array([ 0.8/(1+exp((-55.3-V)/-3.9))+0.2/(1+exp((-74.8-V)/-10.7)) ])
    def __init__(self,g_max, inactivation_fraction = .83,parent=None): ## It has to be like this to avoid passing original array by reference
        Conductance.VoltageDependent.__init__(self,g_max,name = "Shaker",parent=parent)
        self.inactivation_fraction_mode = array([inactivation_fraction]) #0.83
        self.h_order = array([1])

class ShakerWeckstrom_serotonin(ShakerWeckstrom):
    def m_inf(self, V):
         return 1/(  1 + exp((11.7-V)/12.9)  )**(1/3)
    def h_inf(self,V):
        return array([ 0.8/(1+exp((-27.1-V)/-4.3))+0.2/(1+exp((-52.9-V)/-10.5)) ])
    def h_time(self,V):
        return ShakerWeckstrom.h_time(self,V-25)
#    def m_inf(self,V):
#        return ShabWeckstrom.m_inf(self,V-35.4)
    def m_time(self,V):
        return ShakerWeckstrom.m_time(self,V-35.4)


class NovelWeckstrom(Conductance.DelayedRectifier): ## K novel channel, from Weckstrom model
    g_single_channel = 7*1e-9 #Data from only one experiment 5-10 pS; from Hardie(91)
    m_order = 1
    def m_time(self,V):
        return 13 + (6232/(30*sqrt(pi/2)))*exp(-2*((V+19.4)/30)**2)
    def m_inf(self, V):
        return 1/( 1+exp((-14-V)/10.6) )
    def __init__(self,g_max,parent=None):
        Conductance.DelayedRectifier.__init__(self,g_max,name="Novel",parent=parent)
