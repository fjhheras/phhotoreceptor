from __future__ import division
import types
from math import pi,sqrt

import DrosophilaConductance 
from phhotoreceptor import FlyPhotoreceptor

def DrosophilaR16(channel_choice="Vahasoyrinki06MediumLeak", shift="none"):
    reversal_potentials = {'L': 5, 'K' : -85}
    body = FlyPhotoreceptor.CellBody(l=100e-4, r=2e-4, area = 50e-6, reversal_potentials=reversal_potentials) #cm,cm
    axon = None
    if channel_choice == "Vahasoyrinki06LowLeak":
        body.add_leak_conductance(ion_name="K",g=1e-6) #mS,mV
        g_max_Shaker = 12e-6
        g_max_Shab = 46e-6
        g_max_Novel = 1.7e-6

        Shaker={"none": DrosophilaConductance.ShakerWeckstrom(g_max_Shaker),
                "light": DrosophilaConductance.ShakerWeckstrom(g_max_Shaker),
                "serotonin": DrosophilaConductance.ShakerWeckstrom_serotonin(g_max_Shaker)}
        Shab  ={"none": DrosophilaConductance.ShabWeckstrom(g_max_Shab),
                "light": DrosophilaConductance.ShabWeckstrom_LightShifted(g_max_Shab),
                "serotonin": DrosophilaConductance.ShabWeckstrom_serotonin(g_max_Shab)}
        body.add_voltage_channel(Shaker[shift]) #mS
        body.add_voltage_channel(Shab[shift]) #mS
        body.add_voltage_channel(DrosophilaConductance.NovelWeckstrom(g_max_Novel)) #mS

    if channel_choice == "Vahasoyrinki06MediumLeak":
        body.add_leak_conductance(ion_name="K",g=2.1e-6) #mS,mV
        g_max_Shaker = 12e-6
        g_max_Shab = 46e-6
        g_max_Novel = 1.7e-6

        Shaker={"none": DrosophilaConductance.ShakerWeckstrom(g_max_Shaker),
                "light": DrosophilaConductance.ShakerWeckstrom(g_max_Shaker),
                "serotonin": DrosophilaConductance.ShakerWeckstrom_serotonin(g_max_Shaker)}
        Shab  ={"none": DrosophilaConductance.ShabWeckstrom(g_max_Shab),
                "light": DrosophilaConductance.ShabWeckstrom_LightShifted(g_max_Shab),
                "serotonin": DrosophilaConductance.ShabWeckstrom_serotonin(g_max_Shab)}
        body.add_voltage_channel(Shaker[shift]) #mS
        body.add_voltage_channel(Shab[shift]) #mS
        body.add_voltage_channel(DrosophilaConductance.NovelWeckstrom(g_max_Novel)) #mS

    if channel_choice == "Vahasoyrinki06HighLeak":
        body.add_leak_conductance(ion_name="K",g=600e-6) #mS,mV
        g_max_Shaker = 12e-6
        g_max_Shab = 46e-6
        g_max_Novel = 1.7e-6

        Shaker={"none": DrosophilaConductance.ShakerWeckstrom(g_max_Shaker),
                "light": DrosophilaConductance.ShakerWeckstrom(g_max_Shaker),
                "serotonin": DrosophilaConductance.ShakerWeckstrom_serotonin(g_max_Shaker)}
        Shab  ={"none": DrosophilaConductance.ShabWeckstrom(g_max_Shab),
                "light": DrosophilaConductance.ShabWeckstrom_LightShifted(g_max_Shab),
                "serotonin": DrosophilaConductance.ShabWeckstrom_serotonin(g_max_Shab)}
        body.add_voltage_channel(Shaker[shift]) #mS
        body.add_voltage_channel(Shab[shift]) #mS
        body.add_voltage_channel(DrosophilaConductance.NovelWeckstrom(g_max_Novel)) #mS

    if channel_choice == "WeckstromLowLeak":
        body.add_leak_conductance(ion_name="K",g=1e-6) #mS,mV
        g_max_Shaker = 8e-6
        g_max_Shab = 40e-6
        g_max_Novel = 2e-6

        Shaker={"none": DrosophilaConductance.ShakerWeckstrom(g_max_Shaker),
                "light": DrosophilaConductance.ShakerWeckstrom(g_max_Shaker),
                "serotonin": DrosophilaConductance.ShakerWeckstrom_serotonin(g_max_Shaker)}
        Shab  ={"none": DrosophilaConductance.ShabWeckstrom(g_max_Shab),
                "light": DrosophilaConductance.ShabWeckstrom_LightShifted(g_max_Shab),
                "serotonin": DrosophilaConductance.ShabWeckstrom_serotonin(g_max_Shab)}
        body.add_voltage_channel(Shaker[shift]) #mS
        body.add_voltage_channel(Shab[shift]) #mS
        body.add_voltage_channel(DrosophilaConductance.NovelWeckstrom(g_max_Novel)) #mS

    if channel_choice == "Weckstrom":
        body.add_leak_conductance(ion_name="K",g=6.3e-6) #mS,mV
        g_max_Shaker = 8e-6
        g_max_Shab = 40e-6
        g_max_Novel = 2e-6

        Shaker={"none": DrosophilaConductance.ShakerWeckstrom(g_max_Shaker),
                "light": DrosophilaConductance.ShakerWeckstrom(g_max_Shaker),
                "serotonin": DrosophilaConductance.ShakerWeckstrom_serotonin(g_max_Shaker)}
        Shab  ={"none": DrosophilaConductance.ShabWeckstrom(g_max_Shab),
                "light": DrosophilaConductance.ShabWeckstrom_LightShifted(g_max_Shab),
                "serotonin": DrosophilaConductance.ShabWeckstrom_serotonin(g_max_Shab)}
        body.add_voltage_channel(Shaker[shift]) #mS
        body.add_voltage_channel(Shab[shift]) #mS
        body.add_voltage_channel(DrosophilaConductance.NovelWeckstrom(g_max_Novel)) #mS

    if channel_choice == "Niven03":
        a_nm = body.non_microvillar_area #cm2
        body.add_leak_conductance(ion_name="K",g=0.082*a_nm) #mS,mV
        body.add_voltage_channel(DrosophilaConductance.ShakerWeckstrom(1.6*a_nm)) #mS
        body.add_voltage_channel(DrosophilaConductance.ShabWeckstrom(3.5*a_nm)) #mS

    if channel_choice == "passive":
        pass

    return FlyPhotoreceptor.FlyPhotoreceptor(V_rest=-68, body=body, axon=axon)

def PassiveDrosophilaR16WithSameZ_0(photoreceptor):
    '''Create passive photoreceptor with same impedance and V_m rest potential'''
    new_photoreceptor = DrosophilaR16(channel_choice="passive")
    new_photoreceptor.set_steady_state(VZ_0=[photoreceptor.body.V_m, abs(photoreceptor.body.impedance(0))])
    return new_photoreceptor

def PassiveDrosophilaR16WithSameR(photoreceptor):
    '''Create passive photoreceptor with same resistance and V_m rest potential'''
    new_photoreceptor = DrosophilaR16(channel_choice="passive")
    new_photoreceptor.body.add_leak_conductance('K',photoreceptor.body.total_K_conductance(),new_photoreceptor.body.E['K'])
    new_photoreceptor.set_steady_state()
    return new_photoreceptor

def PassiveDrosophilaR16WithBandwidth (bandwidth,V,low_limit_frequency = 0):
    new_photoreceptor = DrosophilaR16(channel_choice="passive")
    C = new_photoreceptor.body.C
    R_same_bandwidth = 1/( 2*pi*C*sqrt(bandwidth**2-2*(low_limit_frequency)**2)  ) #Hz -> MOhm
    new_photoreceptor.set_steady_state(new_photoreceptor.V_rest(), VZ_0=[V, R_same_bandwidth*1000]) #MOhm -> KOhm
    return new_photoreceptor

