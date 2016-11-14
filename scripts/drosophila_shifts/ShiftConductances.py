from __future__ import division
from numpy import *
import DrosophilaConductance


    #Different modes of shifting conductances in Drosophila photoreceptor

def WithSerotonin (photoreceptor, change_LIC_to_keep_depolarisation = False):
    for i,conductance in enumerate(photoreceptor.body.voltage_channels):
        if conductance.channel_name == 'Shaker':
            photoreceptor.body.voltage_channels[i] = DrosophilaConductance.ShakerWeckstrom_serotonin(conductance.g_max,parent=photoreceptor.body)
        elif conductance.channel_name == 'Shab':
            photoreceptor.body.voltage_channels[i] = DrosophilaConductance.ShabWeckstrom_serotonin(conductance.g_max,parent=photoreceptor.body)
    if change_LIC_to_keep_depolarisation:
        photoreceptor.depolarise_with_light(photoreceptor.body.V_m)
    else:
        photoreceptor.relax_to_steady_state()

def WithLight (photoreceptor, change_LIC_to_keep_depolarisation = False):
    for i,conductance in enumerate(photoreceptor.body.voltage_channels):
        if conductance.channel_name == 'Shab':
            photoreceptor.body.voltage_channels[i] = DrosophilaConductance.ShabWeckstrom_LightShifted(conductance.g_max,parent=photoreceptor.body)
    if change_LIC_to_keep_depolarisation:
        photoreceptor.depolarise_with_light(photoreceptor.body.V_m)
    else:
        photoreceptor.relax_to_steady_state()

        
