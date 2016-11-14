from __future__ import division
from numpy import *
import scipy.optimize
import phhotoreceptor.Conductance as Conductance
from phhotoreceptor.ElectricalCompartment import ElectricalCompartment
from phhotoreceptor.Constants import *
#from scipy import *

__author__ = 'Francisco J. H. Heras'

class CellBody(ElectricalCompartment):
    """Variables describing photoreceptor cell body geometry and passive electrical properties, such as length, areas, capacitance..."""
    def __init__(self, l, r, area=None, Cm=1, RA=.08, reversal_potentials = None, lic_current=DEFAULT_LIC_CURRENT, is_there_pump = True):
        self.r = r # body radius (cm)
        self.l = l # body length (cm)
        self.non_microvillar_area = self.l*2*pi*self.r
        self.lic_current=array(lic_current)
        self.is_there_pump = is_there_pump #If True, calculation with pump and exchanger. If false, calculation without.

        if reversal_potentials == None:
            self.E = {'L': 5, 'K' : -85}
        else :
            self.E = reversal_potentials

        if area == None :
            ElectricalCompartment.__init__(self,area = self.non_microvillar_area * 25/4/pi, Cm = Cm, RA =RA)
        else :
            ElectricalCompartment.__init__(self,area = area, Cm = Cm, RA = RA)

        self.set_light_conductance(g=0)
        self.pump = array([-2, 3, 0]) # elements of the array are numbers of ions per cycle.
        self.exchanger = array([0, 3, -1])

    def set_light_conductance(self,g): #it tries to change light conductance. If it fails, creates a new one
        try:
            self.light_conductance.set_g_max_to(g)
        except:
            self.light_conductance = Conductance.LightInduced(g)
            self.light_conductance.add_parent(self)

    def set_steady_state(self,V_rest,VZ_0 = None):
        print("Resetting to dark steady state...")
        self.set_light_conductance(0.0)
        if self.is_there_pump:
            if VZ_0 == None:
                self._set_steady_state_given_leak(V_rest) #No Z_0 given -> only change light leak
            else:
                self._set_steady_state_given_Z_0(V_rest,VZ_0) #Z_0 given -> change light and K leak
        else:
            if VZ_0 == None:
                self._set_steady_state_given_leak_no_pump(V_rest) #TODO
            else:
                print("Not yet implemented... Sorry!")

    def _set_steady_state_given_leak(self,V_r):
        """Sets cell body in steady state, updating pump and the leak conductance with reversal potential of light. For the moment, it accepts no axon input"""
        self.reset_voltage(V_r) #set membrane potential to rest and also tell the channels
        #Short hand notation
        E_L = self.E['L']
        # Calculate sum of all currents (except light leak and light conductance)
        current = zeros_like(self.lic_current)
        for channel in self.voltage_channels:
            current += channel.current_ion(V_r)
        for conductance in self.leak_conductances.values():
            if conductance.ion_name != 'L':
                current += conductance.current_ion(V_r)

        pump_ion_current_fractions = ION_CHARGE*self.pump
        exchanger_ion_current_fractions = ION_CHARGE*self.exchanger
        exchanger_ion_current_fractions = exchanger_ion_current_fractions/sum(exchanger_ion_current_fractions) #in case it is not normalised
        B = matrix([pump_ion_current_fractions, exchanger_ion_current_fractions,self.lic_current]).getT()
        # Now we want to solve the system B*[I_pump,I_exchanger,I_light]' + current = 0, where the unknowns are the I_x
        I_x = B.getI().A.dot(-current)
        self.I_pump = I_x[0]
        self.I_exchanger = I_x[1]
        self.add_leak_conductance(ion_name='L',g=I_x[2]/(E_L-V_r))
        #print("g_L_leak must be ", I_x[2]/(E_L-V_r))

    def _set_steady_state_given_leak_no_pump(self,V_r):
        """Sets cell body with no pump nor exchanger in steady state, updating the leak conductance with reversal potential of light. For the moment, it accepts no axon input"""
        self.reset_voltage(V_r) #set membrane potential to rest and also tell the channels
        #Short hand notation
        E_L = self.E['L']
        # Calculate sum of all currents (except light leak and light conductance)
        current = 0
        for channel in self.voltage_channels:
            current += channel.current(V_r)
        for conductance in self.leak_conductances.values():
            if conductance.ion_name != 'L':
                current += conductance.current(V_r)

        self.I_pump = 0
        self.I_exchanger = 0
        self.add_leak_conductance(ion_name='L',g=-current/(E_L-V_r))
        #print("g_L_leak must be ", -current/(E_L-V_r))


    def _set_steady_state_given_Z_0(self,V_rest, VZ_0):
        """Sets L leak and K leak so cell in the dark has membrane potential self.V_m and low limit of impedance, Z_0 at voltage V. Warning: It clears all leak conductances!"""

        V = VZ_0[0]
        Z_0 =VZ_0[1]
        self.reset_voltage(V) #Change voltage to V in order to match impedances
        self.reset_leak_conductances()

        pump_ion_current_fractions = ION_CHARGE*self.pump
        exchanger_ion_current_fractions = ION_CHARGE*self.exchanger
        exchanger_ion_current_fractions = exchanger_ion_current_fractions/sum(exchanger_ion_current_fractions) #in case it is not normalised

        current_and_inverse_of_Z_0 = zeros(len(self.lic_current)+1)
        for channel in self.voltage_channels:
            current_and_inverse_of_Z_0 += concatenate((channel.current_ion(self.V_m),channel.admitance(f=array([0]),V=self.V_m)))
        current_and_inverse_of_Z_0[3] += -1/Z_0

        B = matrix([concatenate( (pump_ion_current_fractions,[0]) ),
                    concatenate( (exchanger_ion_current_fractions,[0]) ),
                    concatenate( (self.lic_current,[1/(self.E['L'] - self.V_m)]) ),
                    concatenate( (DEFAULT_K_CURRENT,[1/(self.E['K'] - self.V_m)]) )
            ]).getT()
        # B*[I_pump,I_exchanger,I_light,I_K_leak]' = [K_current,Na_current,Ca_current, 1/Z_0_leaks]
        # Now we want to solve the system B*[I_pump,I_exchanger,I_light,I_K_leak]' + [current_channels,1/Z_0_channels] = [0,1/Z_0], where the unknowns are the I_x
        I_x = B.getI().A.dot(-current_and_inverse_of_Z_0)

        self.add_leak_conductance(ion_name='K',g=I_x[3]/(self.E['K'] - self.V_m))
        #print("g_K_leak must be ", I_x[3]/(self.E['K'] - self.V_m))
        self._set_steady_state_given_leak(V_rest) #Once I know K leak, I use this to calculate other leaks

    def depolarise_with_light(self,V_m):
        self.reset_voltage(V_m) #set membrane potential to rest and also tell the channels
        #Short hand notation
        E_L = self.E['L']
        # Calculate sum of all currents (except light conductance)
        if self.is_there_pump:
            current = zeros_like(self.lic_current)
            for channel in self.voltage_channels:
                current += channel.current_ion(V_m)
            for conductance in self.leak_conductances.values():
                current += conductance.current_ion(V_m)

            pump_ion_current_fractions = ION_CHARGE*self.pump
            exchanger_ion_current_fractions = ION_CHARGE*self.exchanger
            exchanger_ion_current_fractions = exchanger_ion_current_fractions/sum(exchanger_ion_current_fractions) #in case it is not normalised
            B = matrix([pump_ion_current_fractions, exchanger_ion_current_fractions,self.lic_current]).getT()
            # Now we want to solve the system B*[I_pump,I_exchanger,I_light]' + current = 0, where the unknowns are the I_x
            I_x = B.getI().A.dot(-current)

            self.I_pump = I_x[0]
            self.I_exchanger = I_x[1]
            self.set_light_conductance(g=I_x[2]/(E_L-V_m))
            #print("g_L must be ", I_x[2]/(E_L-V_m), " that is, ", I_x[2]/(E_L-V_m)/self.leak_conductances['L'].g())
            #print ("Depolarising with light to V=",V_m, " mV")
        else:
            current = 0
            for channel in self.voltage_channels:
                current += channel.current(V_m)
            for conductance in self.leak_conductances.values():
                current += conductance.current(V_m)
            self.set_light_conductance(g=-current/(E_L-V_m))
            #print("g_light ", -current/(E_L-V_m))

    ## Calculate various resistances and currents

    def internally_generated_inward_current(self): #Used in simulations. Keep simple.
        return ElectricalCompartment.internally_generated_inward_current(self) + self.light_conductance.current(self.V_m)

    def inward_current_not_coming_from_conductances(self):
        return self.I_pump + self.I_exchanger

    def inward_current_from_conductances_ion_inf(self,V, dark = False):
        """Steady-state currents from all the conductances. If dark is True, LIC is considered 0"""
        current = self.light_conductance.current_ion_inf(V)
        if dark:
            current = zeros_like(current)
        for channel in self.voltage_channels:
            current += channel.current_ion_inf(V)
        for conductance in self.leak_conductances.values():
            current += conductance.current_ion_inf(V)
        return current

    def total_voltage_independent_conductance(self):
        return self.total_leak_conductance() + self.light_conductance.g()

    def voltage_contrast_gain(self,f):
        """Returns gain from a change in light contrast to voltage"""
        return self.impedance(f) * (self.light_conductance.E - self.V_m) * self.light_conductance.signal_conductance_gain(f)

    def noise_power(self,f,photon_flux):
        return abs( self.impedance(f) * (self.light_conductance.E - self.V_m) * self.light_conductance.g() * self.light_conductance.bump_filter(f) )**2 / photon_flux


class Axon(ElectricalCompartment) :
    """Passive variables of the photoreceptor axon - Not used yet"""
    def __init__(self, l, r, Cm=1, RA=.08, Rm=10, Rt=1000e3):
        self.r = r
        self.l = l
        ElectricalCompartment.__init__(self,area = l*2*pi*r, Cm=Cm, RA=RA)
        self.Rm = Rm               # Axon membrane resistance - 8-50 kOhm cm2
        self.Rt = Rt               #Terminal resistance, kOhm - Hateren estimates 1000kOhm for blowfly

class FlyPhotoreceptor :
    """Fly photoreceptor. Contains constants and an array of conductance"""
    def __init__(self, V_rest, body, axon = None):
        #self.V_rest = V_rest
        self.body = body
        self.axon = axon # Not implemented yet
        self.set_steady_state(V_rest) #Sets pump and leak currents so rest voltage is V_rest

    def set_steady_state(self, V_rest,VZ_0=None) :
        """Set photoreceptor in steady state, updating pump and the leak conductance with reversal potential of light. For the moment, it only accepts K and L leak conductances"""
        self.body.set_steady_state(V_rest, VZ_0)
        #to do: self.axon.set_steady_state
        #self.energy_consumption = - self.body.I_pump* 1e-6 * 6.241e18 # uA -> A -> ATP/s
        #self.surface_taken_by_pump = self.I_pump*1e-6 / (32e-18) * (self.geometry.area_light_insensitive_membrane*1e-8) # cm2 - From Sengupta 2013 (effect cell size...)

    def change_voltage(self,V_m):
        self.body.change_voltage(V_m)

    def reset_voltage(self,V_m):
        self.body.reset_voltage(V_m)

    def what_is_steady_state(self, dark=False):
        def current_ion_inf(p):
            """p[0]=V,p[1]=I_pump,p[2]=I_exchanger"""
            current = self.body.inward_current_from_conductances_ion_inf(p[0], dark = dark)
            current += p[1]*array(ION_CHARGE)*self.body.pump + p[2]*array(ION_CHARGE)*self.body.exchanger
            return current
        p_0 = array([self.body.V_m,self.body.I_pump,self.body.I_exchanger])
        return scipy.optimize.fsolve(current_ion_inf,p_0) #Voltage,I_pump, I_exchanger

    def V_rest(self):
        V_rest = self.what_is_steady_state(dark=True)[0]
        print("V_rest is ", V_rest)
        return V_rest

    def relax_to_steady_state(self):
        p_new = self.what_is_steady_state()
        self.reset_voltage(p_new[0])
        self.body.I_pump = p_new[1]
        self.body.I_exchanger = p_new[2]

    def depolarise_with_light(self,V=None,g_light=None):
        if (g_light == None) and (V != None):
            self.body.depolarise_with_light(V)
        elif (g_light != None) and (V == None):
            self.body.light_conductance.g_max = g_light
            self.relax_to_steady_state()
        else:
            print ("Error! No voltage and no conductance given in depolarization by light command")

    def energy_consumption(self):
        return - self.body.I_pump* 1e-6 * 6.241e18 # uA -> A -> ATP/s


