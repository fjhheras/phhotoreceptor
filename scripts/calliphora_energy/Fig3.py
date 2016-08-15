#!/usr/bin/python3
# Simulation to show that the linearisation of the channel impedance works

from scipy import signal #or was it pylab
from pylab import rfft,figure,loglog,xlabel,ylabel,plot,show
from numpy import *
import phhotoreceptor.Linearise as Linearise
from phhotoreceptor.DepolarisePhotoreceptor import DepolarisePhotoreceptor
import phhotoreceptor.Experiment as Experiment
import FlyFactory as FlyFactory

__author__ = 'Francisco J. H. Heras'

T         = 1000  #ms
dt        = 0.01 #ms
time      = arange(0,T,dt) #ms
f         = arange(0,1000/2/dt+1000/T,1000/T) #Hz
repetitions = 30
sample = 2

cutoff_freq = 1000 #Hz

Vr = [-60,-52,-44,-37]
colour_graph=['b','g','r','c']
SD_noise  = 0.7e-3 #nA -> uA

phasearray = vectorize (lambda z : angle(z))

def butter_lowpass(lowcut, fs, order=6):
    nyq = 0.5 * fs
    low = lowcut / nyq
    print(low)
    b, a = signal.butter(order, low, btype='lowpass')
    return b, a, floor(lowcut)

def give_me_white_noise(standard_deviation, cutoff_freq, T, dt): #Hz,ms,ms
    time      = arange(0,T,dt) #ms
    I_nofilter= array(random.randn(len(time)))
    b, a, cut_off_i = butter_lowpass(cutoff_freq, 1000/dt) #40-th order Butterworth filter
    I=signal.lfilter(b, a, I_nofilter)
    I=I-mean(I)
    return I*standard_deviation/sqrt(I.var()) #uA

def run_WN(photoreceptor, I):

    Z_dictionary = {}
    V_dictionary = {}
    Vl_dictionary = {}
    I_dictionary = {}

    for i,V in enumerate(Vr):
        DepolarisePhotoreceptor.WithLight(photoreceptor,V)
        ## Filtering signal
        print("Processing V=",V, " mV")


        ## Run HH

        print("Running HH...")
        Vm, g = Experiment.inject_current(photoreceptor,I,dt)
        Vl = Linearise.inject_current(photoreceptor,I,dt)
        print("Variance of voltage is ",var(Vm), "mV2")
        ## Calculate impedance

        Vmm = (Vm-mean(Vm)) #without DC
        CPSD = zeros(round(len(time)/2+1),dtype=complex_)
        PSD = zeros(round(len(time)/2+1))
        window = hamming(len(time))
        print("len(window)",len(window))
        for ii in range(repetitions):
            I_cut = I[ii*len(time):(ii+1)*len(time)]
            V_cut = Vmm[ii*len(time):(ii+1)*len(time)]
            PSD += power(absolute(rfft(I_cut*window)),2)
            CPSD += rfft(V_cut*window)*rfft(I_cut*window).conjugate()

        Z = CPSD/PSD #MOhm

        Z_dictionary[V] = Z
        I_dictionary[V] = I[sample*len(time):(sample+1)*len(time)]
        V_dictionary[V] = Vm[sample*len(time):(sample+1)*len(time)]
        Vl_dictionary[V] = Vl[sample*len(time):(sample+1)*len(time)]
    return (Z_dictionary, V_dictionary,Vl_dictionary, I_dictionary)

if __name__ == '__main__':

    photoreceptor = FlyFactory.CalliphoraR16(channel_choice='Anderson')
    I = give_me_white_noise(SD_noise, cutoff_freq, T * repetitions, dt)
    Z_dictionary, V_dictionary, Vl_dictionary, I_dictionary = run_WN(photoreceptor,I)

    fig2 = figure(2)
    axV = fig2.add_subplot(2,1,1)
    axV.set_xlabel("Time (ms)")
    axV.set_ylabel("Membrane voltage (mV)")
    axI = fig2.add_subplot(2,1,2)

    for i,V in enumerate(Vr):
        DepolarisePhotoreceptor.WithLight(photoreceptor,V)
        figure(1)
        f_plot = f[1:len(f)/10]
        Z = photoreceptor.body.impedance(f_plot)
        Z_abs = abs(Z_dictionary[V])[1:len(f)/10]
        print(len(f_plot),len(Z_abs/1000))
        loglog(f_plot,Z_abs/1000,colour_graph[i])
        loglog(f_plot,abs(Z)/1000,'k'+'--')
        xlabel("Frequency (Hz)")
        ylabel("Impedance (MOhm)")



        axV.plot(time,V_dictionary[V],colour_graph[i])
        axV.plot(time,Vl_dictionary[V]+V,'k'+'--')

        axI.plot(time,I_dictionary[V]*1e3,'k') #uA->nA

show()

