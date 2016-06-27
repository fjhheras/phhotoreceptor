#!/usr/bin/python3

from pylab import *
from numpy import *
import Drone as Drone
from phhotoreceptor.DepolarisePhotoreceptor import DepolarisePhotoreceptor

def Q_value(Z,f):
    Z_max = abs(Z[0])
    for i,ff in enumerate(f):
        if abs(Z[i]) > Z_max:
            Z_max = abs(Z[i])
    return Z_max/abs(Z[0])


####### BODY STARTS HERE

f_low = 0 #Hz
f_medium = 5 #Hz

fig1 = figure(1)
ax_Z = fig1.add_subplot(1,1,1)


V=-38.0
delta_f = 0.1
f = arange(.2,500,delta_f)
f_from_medium = arange(f_medium,500,1)

colour_graph=['r','b']
k_h_values = [1,0.1]

for i_k_h,k_h in enumerate(k_h_values):

    HH  = Drone.Vallet92(k_h=k_h)
    DepolarisePhotoreceptor.WithLight(HH,V)
    Z = HH.body.impedance(f) #All frequencies
    print ("When k_h=",k_h, ", Q value = ",Q_value(Z,f))

    label_str = 'k_h = ' + str(k_h)
    ax_Z.loglog(f,abs(Z)/1000,colour_graph[i_k_h],linewidth=2,label = label_str)

ax_Z.set_xlabel("Frequency (Hz)")
ax_Z.set_ylabel("Impedance (MOhms)")
ax_Z.legend(loc=3,prop={'size':12})


show()
