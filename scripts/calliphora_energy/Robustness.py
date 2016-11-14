#!/usr/bin/python3

from pylab import *
from numpy import *
import FlyFactory
from phhotoreceptor.DepolarisePhotoreceptor import DepolarisePhotoreceptor
from GBWPutils import Gain_Bandwidth

HH  = FlyFactory.CalliphoraR16(channel_choice = "Anderson") #

bw_factor = True # True -> plots bandwidth improvement, False -> plots resistance ratio instead

####### BODY STARTS HERE

fig1 = figure(1,figsize=[9,5])
ax_dark = fig1.add_subplot(1,2,1)
ax_light = fig1.add_subplot(1,2,2)

V_rest  = -60
V_light = -37

factor_dark = []
factor_light = []
factor_bw_dark = []
factor_bw_light = []
input_r_dark = []
input_r_light = []

factor_dark_u = []
factor_light_u = []
factor_bw_dark_u = []
factor_bw_light_u = []
input_r_dark_u = []
input_r_light_u = []

perturbed_ch_ = ["a","b","m_order","g_max","tau"]
perturbation_ = array([0.8,0.85,0.9,0.95,1.05,1.1,1.15,1.2])

for channel in HH.body.voltage_channels:
    for perturbed in perturbed_ch_:
        original = getattr(channel, perturbed)
        for perturbation in perturbation_:
            setattr(channel,perturbed,original*perturbation)

            HH.set_steady_state(V_rest)
            input_r = HH.body.impedance(0)
            R = HH.body.resistance()
            pippo, Bandwidth = Gain_Bandwidth(HH.body.impedance, f_min=0)
            factor_bw_dark.append(2*pi*R*HH.body.C* 1e-3*Bandwidth)
            factor_dark.append(R/input_r)
            input_r_dark.append(input_r/1000) #kOhm -> MOhm

            DepolarisePhotoreceptor.WithLight(HH,V_light)
            input_r = HH.body.impedance(0)
            R = HH.body.resistance()
            pippo, Bandwidth = Gain_Bandwidth(HH.body.impedance, f_min=0)
            factor_bw_light.append(2*pi*R*HH.body.C* 1e-3*Bandwidth)
            factor_light.append(R / input_r)
            input_r_light.append(input_r/1000) #kOhm -> MOhm

        setattr(channel, perturbed, original)

original = HH.body.leak_conductances['K'].g_max
for perturbation in perturbation_:
    HH.body.leak_conductances['K'].g_max = original*perturbation

    HH.set_steady_state(V_rest)
    input_r = HH.body.impedance(0)
    R = HH.body.resistance()
    pippo, Bandwidth = Gain_Bandwidth(HH.body.impedance, f_min=0)
    factor_bw_dark.append(2 * pi * R * HH.body.C * 1e-3* Bandwidth)
    factor_dark.append(R/input_r)
    input_r_dark.append(input_r/1000) #kOhm -> MOhm

    DepolarisePhotoreceptor.WithLight(HH,V_light)
    input_r = HH.body.impedance(0)
    R = HH.body.resistance()
    pippo, Bandwidth = Gain_Bandwidth(HH.body.impedance, f_min=0)
    factor_bw_light.append(2 * pi * R * HH.body.C * 1e-3* Bandwidth)
    factor_light.append(R / input_r)
    input_r_light.append(input_r/1000) #kOhm -> MOhm

HH.body.leak_conductances['K'].g_max = original

HH.set_steady_state(V_rest)
input_r = HH.body.impedance(0)
R = HH.body.resistance()
pippo, Bandwidth = Gain_Bandwidth(HH.body.impedance, f_min=0)
factor_bw_dark_u.append(2 * pi * R * HH.body.C * 1e-3 * Bandwidth)
factor_dark_u.append(R / input_r)
input_r_dark_u.append(input_r / 1000)  # kOhm -> MOhm

DepolarisePhotoreceptor.WithLight(HH, V_light)
input_r = HH.body.impedance(0)
R = HH.body.resistance()
pippo, Bandwidth = Gain_Bandwidth(HH.body.impedance, f_min=0)
factor_bw_light_u.append(2 * pi * R * HH.body.C * 1e-3* Bandwidth)
factor_light_u.append(R / input_r)
input_r_light_u.append(input_r / 1000)  # kOhm -> MOhm


if bw_factor:
    ax_dark.plot(input_r_dark,factor_bw_dark,'+')
    ax_light.plot(input_r_light,factor_bw_light,'+')
    ax_dark.plot(input_r_dark_u,factor_bw_dark_u,'r+')
    ax_light.plot(input_r_light_u,factor_bw_light_u,'r+')
    ax_dark.set_ylabel('Bandwidth improvement factor')
    ax_light.set_ylabel('Bandwidth improvement factor')

else:
    ax_dark.plot(input_r_dark,factor_dark,'+')
    ax_light.plot(input_r_light,factor_light,'+')
    ax_dark.plot(input_r_dark_u,factor_dark_u,'r+')
    ax_light.plot(input_r_light_u,factor_light_u,'r+')
    ax_dark.set_ylabel('Ratio between membrane resistance and input resistance')
    ax_light.set_ylabel('Ratio between membrane resistance and input resistance')

ax_dark.set_title('Dark adapted')
ax_light.set_title('Light adapted')
ax_dark.set_xlabel('Input resistance (MOhm)')
ax_light.set_xlabel('Input resistance (MOhm)')

show()

