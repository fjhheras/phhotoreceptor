#!/usr/bin/python3
## Prints a cartoon of the impedance of a cell dominated by either a slow or a fast non-inactivating channel
import matplotlib.pyplot as plt
from numpy import *


option = 2 # 2 fast, 1 slow, 0 RC

R = 1
r = 2

if option == 1:
    tau_channel = 20
elif option == 2:
    tau_channel = 1
else:
    r = 0
    tau_channel = 1

tau_membrane = 1
C = tau_membrane/R
L = tau_channel*r

f_channel = 1/2/pi/tau_channel
f_membrane = 1/2/pi/tau_membrane
f = logspace(-3,0.5)
print(f)

if option == 0:
    Z = vectorize ( lambda f : abs(1/(1/R + 1j*2*pi*f*C)) )
else:
    Z = vectorize ( lambda f : abs(1/(1/R + 1/(r+1j*2*pi*f*L) + 1j*2*pi*f*C)) )
    Z_no_C = vectorize ( lambda f : abs(1/(1/R + 1/(r+1j*2*pi*f*L))) )
    y_no_C = Z_no_C(f)

y = Z(f)
fig, ax = plt.subplots(figsize=[3.307,2])
plt.plot(log10(f), log10(y), 'r', linewidth=2)

if option == 1:
    plt.plot(log10(f), log10(y_no_C), 'r--', linewidth=2)
    ax.set_xticks((log10(f_channel), log10(f_membrane)))
    ax.set_xticklabels(('$f_{channel}$', '$f_{membrane}$'))
    ax.set_yticks(( log10(R) , log10( 1/(1/R+1/r) )  ))
    ax.set_yticklabels(('$R$', '$R||r$'))
elif option == 2:
    plt.plot(log10(f), log10(y_no_C), 'r--', linewidth=2)
    ax.set_xticks([log10(f_channel)])
    ax.set_xticklabels(['$f_{channel}=f_{membrane}$'])
    ax.set_yticks(( log10(R) , log10( 1/(1/R+1/r) )  ))
    ax.set_yticklabels(('$R$', '$R||r$'))
else:
    ax.set_xticks([log10(f_membrane)])
    ax.set_xticklabels(['$f_{membrane}$'])
    ax.set_yticks([ log10(R) ])
    ax.set_yticklabels(['$R$'])

ax.xaxis.set_ticks_position('bottom')
ax.set_xlim([-3.2,0.5])
ax.set_xlabel("Frequency")

ax.yaxis.set_ticks_position('left')
ax.set_ylim([log10(R)-1,log10(R)+0.1])
ax.set_ylabel("Impedance")

ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
plt.show()
