from __future__ import division
from scipy import pi
import phhotoreceptor.Conductance as Conductance
__author__ = 'Francisco J. H. Heras'

class ElectricalCompartment:
    """Electrical compartment. It is the basic building block, and then it is inherited in axons and bodies. Can have conductances, conductance, etc"""
    def __init__(self, area, Cm, RA):
        self.area = area #cm2
        self.Cm = Cm #uF/cm2
        self.C = Cm * self.area #uF
        self.RA = RA # specific intracellular resistivity (kOhm*cm) .1-> standard; .08 -> Hateren&Laughlin 90 (for LMC)
        self.reset_voltage_channels()
        self.reset_leak_conductances() #It should be a diccionary of Conductance

    def reset_voltage_channels(self): #Clears the array of conductances
        self.voltage_channels = []
    def add_voltage_channel(self,voltage_channel): #Adds a new voltage channel to the photoreceptor
        self.voltage_channels.append(voltage_channel)
        voltage_channel.add_parent(self) #Parent can be also appended when the channel is created

    def reset_leak_conductances(self): #Clears the dictionary of leak conductances
        self.leak_conductances = {}
    def add_leak_conductance(self,ion_name,g): #Adds a new leak conductance to the dictionary
        leak_conductance = Conductance.VoltageIndependent(ion_name,g)
        self.leak_conductances[ion_name] = leak_conductance
        leak_conductance.add_parent(self)

    def update_membrane_voltage(self, I_total, dt):
        dV = I_total / self.C * dt
        self.V_m += dV
        return self.V_m

    def change_voltage(self,V_m):
        self.V_m = V_m #change membrane potential

    def reset_voltage(self,V_m):
        self.V_m = V_m #change membrane potential
        for channel in self.voltage_channels:
            channel.initialise_mh() # and reset all channels to their rest value at that voltage

    def freeze_conductance(self,index):
        if index == None:
            self.freeze_all_conductances()
        else :
            self.voltage_channels[index].freeze_conductance()
    def freeze_inactivation(self,index):
        if index == None:
            self.freeze_all_inactivations()
        else :
            self.voltage_channels[index].freeze_inactivation()
    def freeze_all_conductances(self):
        for channel in self.voltage_channels:
            channel.freeze_conductance()
    def unfreeze_all_conductances(self):
        for channel in self.voltage_channels:
            channel.unfreeze_conductance()
    def freeze_all_inactivations(self):
        for channel in self.voltage_channels:
            channel.freeze_inactivation()
    def unfreeze_all_inactivations(self):
        for channel in self.voltage_channels:
            channel.unfreeze_inactivation()



    ## Calculate various resistances and currents

    def voltage_dependent_inward_current(self):
        current = 0
        for channel in self.voltage_channels:
            current += channel.current(self.V_m)
        return current

    def inward_current_not_coming_from_conductances(self):
        return 0

    def internally_generated_inward_current(self):
        current = self.inward_current_not_coming_from_conductances()
        for channel in self.voltage_channels:
            current += channel.current(self.V_m)
        for conductance in self.leak_conductances.values():
            current += conductance.current(self.V_m)
        return current

    def total_leak_conductance(self):
        total = 0
        for leak in self.leak_conductances.values():
            total += leak.g()
        return total #mS

    def total_voltage_independent_conductance(self):
        return self.total_leak_conductance()

    def total_voltage_dependent_conductance(self):
        total = 0
        for channel in self.voltage_channels:
            total += channel.g() #mS
        return total

    def total_K_conductance(self): #temporary: in use in some scripts
        total = 0
        for channel in self.voltage_channels:
            total += channel.g()*channel.ion_current_fraction[0] #mS
        for conductance in self.leak_conductances.values():
            total += conductance.g()*conductance.ion_current_fraction[0] #mS
        return total

    def total_voltage_dependent_conductance_inf(self,V): # at rest
        total = 0
        for channel in self.voltage_channels:
            total += channel.g_inf(V) #mS
        return total

    def resistance(self):
        return 1/(self.total_voltage_dependent_conductance() + self.total_voltage_independent_conductance())

    def impedance(self,f): #f-> Frequency in Hz, V-> Membrane voltage
        """Calculates impedances (using the RrL approximation)"""
        omega_ms = f*2*pi/1000
        g_array_K = self.C*omega_ms*1j + self.total_voltage_independent_conductance() #mS from light+leak and capacitance(uF)
        for channel in self.voltage_channels:
            g_array_K += channel.admitance(f, self.V_m)
        return 1/g_array_K #kOhm
