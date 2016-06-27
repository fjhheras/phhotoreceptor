from pylab import *
from numpy import *

__author__ = 'Francisco J. H. Heras'

def freeze_conductances(photoreceptor,index=None): #If index is None -> freeze all conductances
    photoreceptor.body.freeze_conductance(index=index)
def unfreeze_conductances(photoreceptor):
    photoreceptor.body.unfreeze_all_conductances()
def freeze_inactivations(photoreceptor,index=None):
    photoreceptor.body.freeze_inactivation(index=index)
def unfreeze_inactivations(photoreceptor):
    photoreceptor.body.unfreeze_all_inactivations()

def modify_conductance(photoreceptor,channel_name,multiplier,change_LIC_to_keep_depolarisation=False):
    """Multiplies the conductance of name channel_name by a fraction (multiplier). Then, it either changes LIC or relaxes voltage to steady state."""
    for conductance in photoreceptor.body.voltage_channels:
        if conductance.channel_name == channel_name:
            conductance.set_g_max_to(conductance.g_max * multiplier)
    if change_LIC_to_keep_depolarisation:
        photoreceptor.depolarise_with_light(photoreceptor.body.V_m)
    else:
        photoreceptor.relax_to_steady_state()


def inject_current(photoreceptor,injected_current,dt):

    V = zeros_like(injected_current)
    V[0] = photoreceptor.body.V_m
    g_Ch = []
    for i,channel in enumerate(photoreceptor.body.voltage_channels) :
        g_Ch.append(zeros_like(injected_current))
        g_Ch[i][0] = channel.g()

    print("Initial voltage = ", V[0])

    for i in range(1,len(injected_current)):
        V[i] = photoreceptor.body.update_membrane_voltage(photoreceptor.body.internally_generated_inward_current() + injected_current[i-1],dt)
        for ii,channel in enumerate(photoreceptor.body.voltage_channels):
            channel.update(V[i],dt)
            g_Ch[ii][i] = channel.g()


    return V, g_Ch


def voltage_clamp(photoreceptor,voltage_command,dt):

    I = zeros_like(voltage_command)
    print("Initial voltage = ", voltage_command[0])
    photoreceptor.body.reset_voltage(voltage_command[0])
    I[0] = - photoreceptor.body.internally_generated_inward_current() #clamp must oppose the sum of all internal currents
    g_Ch = []
    for i,channel in enumerate(photoreceptor.body.voltage_channels) :
        g_Ch.append(zeros_like(voltage_command))
        g_Ch[i][0] = channel.g()

    for i in range(1,len(voltage_command)):
        photoreceptor.body.change_voltage(voltage_command[i])
        I[i] = - photoreceptor.body.internally_generated_inward_current() #clamp must oppose the sum of all internal currents
        #V[i] = photoreceptor.body.update_membrane_voltage(photoreceptor.body.internally_generated_inward_current() + injected_current[i-1],dt)
        for ii,channel in enumerate(photoreceptor.body.voltage_channels):
            channel.update(voltage_command[i],dt)
            g_Ch[ii][i] = channel.g()

    return I, g_Ch

def voltage_clamp_leak_substracted(photoreceptor,voltage_command,dt):

    I = zeros_like(voltage_command)
    print("Initial voltage = ", voltage_command[0])
    photoreceptor.body.reset_voltage(voltage_command[0])
    I[0] = - photoreceptor.body.voltage_dependent_inward_current() #clamp must oppose the sum of all internal currents
    g_Ch = []
    for i,channel in enumerate(photoreceptor.body.voltage_channels) :
        g_Ch.append(zeros_like(voltage_command))
        g_Ch[i][0] = channel.g()

    for i in range(1,len(voltage_command)):
        photoreceptor.body.change_voltage(voltage_command[i])
        I[i] = - photoreceptor.body.voltage_dependent_inward_current() #clamp must oppose the sum of all internal currents
        #V[i] = photoreceptor.body.update_membrane_voltage(photoreceptor.body.internally_generated_inward_current() + injected_current[i-1],dt)
        for ii,channel in enumerate(photoreceptor.body.voltage_channels):
            channel.update(voltage_command[i],dt)
            g_Ch[ii][i] = channel.g()

    return I, g_Ch

