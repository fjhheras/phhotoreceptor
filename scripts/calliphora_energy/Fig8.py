#!/usr/bin/python3

from pylab import *
from numpy import *
import FlyFactory
from phhotoreceptor.DepolarisePhotoreceptor import DepolarisePhotoreceptor
from GBWPutils import GBWP, Gain_Bandwidth

def classic_cost_estimation(R,V,E_K,E_L):
    '''Estimates cost as in e.g. Niven, Anderson and Laughlin (2007), using an estimate
    of the membrane resistance R, the membrane voltage V and the reversal potentials
    of potassium E_K and the light-induced current E_L'''
    g_total_K = -(V- E_L)/(V/2-E_K*3/2 + E_L) / R; #mS
    I_pump = g_total_K*(V - E_K)/2 #uA
    return I_pump* 1e-6 * 6.241e18 # uA -> A -> ATP/s

####### BODY STARTS HERE

option_debugging = False

HH  = FlyFactory.CalliphoraR16(channel_choice = "Anderson") #
E_K = HH.body.E['K']
E_L = HH.body.E['L']

figure(1,figsize=[3.307,2])
xlim([-62, -35])
ylim([8e8, 2e10])


###### CONTINUOUS ACROSS VOLTAGES
Vr = arange(-60.0,-36.9,0.1)
cost = zeros_like(Vr)
good_estimation = zeros_like(Vr)
estimated_cost_from_peak = zeros_like(Vr)
estimated_cost_from_ir = zeros_like(Vr)
total_K_conductance = zeros_like(Vr)
for i,V in enumerate(Vr): #Impedances at lowest frequency
    DepolarisePhotoreceptor.WithLight(HH,V)
    cost[i] = HH.energy_consumption()
    R = HH.body.resistance() # Good estimation! (Cheating)
    good_estimation[i] = classic_cost_estimation(R,V,E_K,E_L)

    R,pippo = Gain_Bandwidth(HH.body.impedance) # R estimated as peak impedance. Wrong!
    estimated_cost_from_peak[i] = classic_cost_estimation(R,V,E_K,E_L)

    R = abs(HH.body.impedance(0)) # R estimated as input resistance. Even worse!
    estimated_cost_from_ir[i] = classic_cost_estimation(R,V,E_K,E_L)

semilogy(Vr,cost)
semilogy(Vr,estimated_cost_from_ir,'--')
semilogy(Vr,estimated_cost_from_peak,'--')
if option_debugging:
    semilogy(Vr,good_estimation,'--')
xlabel("Voltage (mV)")
ylabel("Cost (ATP/s)")
title("Photoreceptor cost at different light levels")

show()
