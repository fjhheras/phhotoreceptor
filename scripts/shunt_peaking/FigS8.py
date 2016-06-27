#!/usr/bin/python3
#INCOMPLETE

from pylab import *
from numpy import *
import FlyFactory as FlyFactory
from phhotoreceptor.DepolarisePhotoreceptor import DepolarisePhotoreceptor
import phhotoreceptor.Linearise as Linearise
from GBWPutils import GBWP

option_debugging = False

####### BODY STARTS HERE

fig1 = figure(1) # Plot comparing against RC
ax_max_gbwp = fig1.add_subplot(121)
ax_taumax_gbwp = fig1.add_subplot(122)

HH  = FlyFactory.CalliphoraR16(channel_choice = "Weckstrom")
passive_gbwp = 1/2/pi/HH.body.C

### DARK ADAPTED

RrLC = Linearise.give_RrLC_equivalent(HH.body)

tau_fast_a = linspace(.10,5,400) #ms
R_a = linspace(RrLC.R/20,RrLC.R,20) #mS

gbwp_max = zeros_like(R_a)
gbwp_maxi = zeros_like(R_a)
gbwp = zeros_like(tau_fast_a)

for ii,R in enumerate(R_a):
    RrLC.R = R
    print("R=",R)
    for i,tau in enumerate(tau_fast_a):
        RrLC.L[0] = RrLC.r[0]*tau
        Z = lambda f: RrLC.impedance(f)
        gbwp[i] = GBWP(Z)/passive_gbwp

    gbwp_max[ii] = max(gbwp)
    gbwp_maxi[ii] = tau_fast_a[argmax(gbwp)]

ax_max_gbwp.plot(R_a/1e3, gbwp_max)
ax_taumax_gbwp.plot(R_a/1e3, gbwp_maxi)

### LIGHT ADAPTED

V = -40
#HH  = FlyFactory.CalliphoraR16(channel_choice = "Weckstrom")
DepolarisePhotoreceptor.WithLight(HH,V)
RrLC = Linearise.give_RrLC_equivalent(HH.body)

tau_fast_a = linspace(.10,2,200) #ms
R_a = linspace(RrLC.R/20,RrLC.R,20) #mS


gbwp_max = zeros_like(R_a)
gbwp_maxi = zeros_like(R_a)
gbwp = zeros_like(tau_fast_a)

for ii,R in enumerate(R_a):
    RrLC.R = R
    print("R=",R)
    for i,tau in enumerate(tau_fast_a):
        RrLC.L[0] = RrLC.r[0]*tau
        Z = lambda f: RrLC.impedance(f)
        gbwp[i] = GBWP(Z)/passive_gbwp

    gbwp_max[ii] = max(gbwp)
    gbwp_maxi[ii] = tau_fast_a[argmax(gbwp)]

ax_max_gbwp.plot(R_a/1e3, gbwp_max,'r')
ax_taumax_gbwp.plot(R_a/1e3, gbwp_maxi,'r')

ax_max_gbwp.set_xlabel("Membrane resistance (MOhms)")
ax_taumax_gbwp.set_xlabel("Membrane resistance (MOhms)")
ax_max_gbwp.set_ylabel("Maximum relative GBWP")
ax_taumax_gbwp.set_ylabel("Optimum FDR activation time constant (ms)")


show()
