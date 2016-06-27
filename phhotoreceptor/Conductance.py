from __future__ import division
__author__ = 'Francisco J. H. Heras'
from numpy import *
import types
from phhotoreceptor.Constants import *


class VoltageIndependent:
    """Conductance that is not voltage dependent"""
    def __init__(self, ion_name, g, parent=None) :
        self.g_max = g #conductance, mS
        self.ion_name = ion_name
        if parent!=None:
            self.add_parent(parent)
        self.ion_current_fraction = zeros_like(ION_CHARGE)
        if ion_name in ION2NUMBER :
            self.ion_current_fraction[ION2NUMBER[ion_name]]=1.0

    parent = None
    def add_parent(self,parent):
        """Links conductance to the parent (ElectricalCompartment). Used to see the membrane voltage and the reversal potentials."""
        self.parent = parent
        if self.ion_name == 'L':
            self.ion_current_fraction = parent.lic_current
        if self.ion_name in parent.E:
            self.E = parent.E[self.ion_name]
            #print("Adding ", self.ion_name, " reversal potential of ", parent.E[self.ion_name])
        else:
            print ("Error: Ion name not in reversal potential dictionary!")

    def g(self):
        return self.g_max
    def set_g_max_to(self,g):
        self.g_max = g
    def g_inf(self,V_m = None):
        return self.g_max
    def current(self,V_m = None):
        """Gives current in uA. Positive means inward"""
        if V_m == None:
            return (self.E - self.parent.V_m)*self.g()
        else:
            return (self.E - V_m)*self.g()
    def current_ion(self,V_m = None):
        """Gives current in uA separated among the different ions"""
        return self.current(V_m) * self.ion_current_fraction
    def current_inf(self,V_m):
        """Gives current in uA. Positive means inward"""
        return (self.E - V_m)*self.g_inf(V_m)
    def current_ion_inf(self,V_m):
        """Gives current in uA separated among the different ions"""
        return self.current_inf(V_m) * self.ion_current_fraction

class LightInduced(VoltageIndependent):
    """Light induced conductance"""
    def __init__(self, g) :
        VoltageIndependent.__init__(self,ion_name='L',g=g)
    def bump_filter(self,f):
        """Filtering due to the bump shape. It should be normalised so it is 1 at low frequencies"""
        return 1
    def delay_filter(self,f):
        """Filtering due to stochastiticy in bump delay. It should be normalised so it is 1 at low frequencies"""
        return 1
    def signal_conductance_gain(self,f):
        return self.g() * self.delay_filter(f) * self.bump_filter(f)


#Dictionary of functions to substitute VoltageDependent.update to freeze and unfreeze conductances
def _update_00(self,V,dt):
    self.update_m(V,dt)
    self.update_h(V,dt)
def _update_10(self,V,dt):
    self.update_h(V,dt)
def _update_01(self,V,dt):
    self.update_m(V,dt)
def _update_11(self,V,dt):
    return 0
_update_functions = {
    (False,False): _update_00,
    (True,False): _update_10,
    (False,True): _update_01,
    (True,True): _update_11
}

class VoltageDependent(VoltageIndependent):
    def __init__(self,g_max,ion_name='K',name=None,parent=None): #by default, it is a K channel
        VoltageIndependent.__init__(self,ion_name=ion_name,g=g_max,parent=parent) #mV -- for the moment we only use K conductance
        self.m = 0
        self.h = array([])    #Contains the h of different modes of inactivation
        self.frozen = (False,False)                   #(ActivationFrozen,InactivationFrozen)
        self.update = types.MethodType(_update_functions[self.frozen], self)
        self.channel_name = name

    ## There are two options to define conductance properties
    ## Both are defined circularly in function of the other in VoltageDependent
    ## In child classes, either inf and time is defined or beta and alpha.
    # Option 1: Define alpha, beta and calculate times and steady state conductances
    def m_inf(self,V):
        return self.alpha_m(V)/(self.alpha_m(V) + self.beta_m(V))
    def h_inf(self,V):
        return self.alpha_h(V)/(self.alpha_h(V) + self.beta_h(V))
    def m_time(self,V): #ms
        return 1/(self.alpha_m(V) + self.beta_m(V))
    def h_time(self,V): #ms
        return 1/(self.alpha_h(V) + self.beta_h(V))
    # Option 2: Define times and steady state conductances and calculate alpha beta
    def alpha_m(self,V):
        return self.m_inf(V) / self.m_time(V)
    def beta_m(self,V):
        return (1 - self.m_inf(V)) / self.m_time(V)
    def alpha_h(self,V):
        return self.h_inf(V) / self.h_time(V)
    def beta_h(self,V):
        return (1 - self.h_inf(V) ) / self.h_time(V)


    ## Derivatives for calculating impedances with the linear transform
    def alpha_m_prime(self,V): #alpha derivative, for impedance
        dV=1e-6
        return (self.alpha_m(V+dV)-self.alpha_m(V))/dV
    def beta_m_prime(self,V): #beta derivative, for impedance
        dV=1e-6
        return (self.beta_m(V+dV)-self.beta_m(V))/dV
    def alpha_h_prime(self,V): #alpha derivative, for impedance
        dV=1e-6
        return (self.alpha_h(V+dV)-self.alpha_h(V))/dV
    def beta_h_prime(self,V): #beta derivative, for impedance
        dV=1e-6
        return (self.beta_h(V+dV)-self.beta_h(V))/dV
    ##Channel conductance, steady state or instantanous
    def g_inf(self,V):
        """Steady state conductance at voltage V"""
        if self.frozen[0]:
            m_inf = self.m #If activation is frozen,take the m value
        else:
            m_inf = self.m_inf(V) #Otherwise m_inf describes it well

        if self.h.size == 0: # means the channel does not inactivate
            return self.g_max * (m_inf ** self.m_order)
        else:
            if self.frozen[1]:
                h_inf_array = self.h #If inactivation is frozen,take the current h value
            else:
                h_inf_array = self.h_inf(V) #Otherwise h_inf describes it well

            Gh = (1-sum(self.inactivation_fraction_mode)) #fraction that fails to inactivate
            for i,h_inf in enumerate(h_inf_array):
                #if self.h_order[i] > 0: # if there is inactivation
                Gh = Gh + self.inactivation_fraction_mode[i] * (h_inf ** self.h_order[i])
            return self.g_max * (m_inf ** self.m_order) * Gh
    def g(self):
        """Conductance with the current, assumed up-to-date, self.m and self.h"""
        if self.h.size == 0: # means the channel does not inactivate
            return self.g_max * (self.m ** self.m_order)
        else:
            h_array = self.h
            Gh = (1-sum(self.inactivation_fraction_mode)) #fraction that fails to inactivate
            for i,h in enumerate(h_array):
                #if self.h_order[i] > 0: # if there is inactivation
                Gh = Gh + self.inactivation_fraction_mode[i] * (h ** self.h_order[i])
            return self.g_max * (self.m ** self.m_order) * Gh

    def initialise_mh(self,V=None):
        if V==None:
            V = self.parent.V_m
        if self.frozen[0]==False:
            self.m = self.m_inf(V)
        if self.frozen[1]==False:
            self.h = self.h_inf(V)
        #print("h=",self.h)
    def update_m(self,V,dt):
        self.m = min(max(0,self.m + dt*(self.alpha_m(V)*(1 - self.m) - self.beta_m(V)*self.m)),1) #done here because activation faster than inact.
    def update_h(self,V,dt):
        self.h = self.h + dt*(self.alpha_h(V)*(1 - self.h) - self.beta_h(V)*self.h)

    def freeze_activation(self):
        if self.frozen[0]:
            print("Warning: freezing activation that was already frozen!")
        self.frozen = (True,self.frozen[1])
        self.update = types.MethodType(_update_functions[self.frozen], self)
    def unfreeze_activation(self):
        if self.frozen[0]==False:
            print("Warning: unfreezing activation that was not frozen!")
        self.frozen = (False,self.frozen[1])
        self.update = types.MethodType(_update_functions[self.frozen], self)
    def freeze_inactivation(self):
        if self.frozen[1]:
            print("Warning: freezing inactivation that was already frozen!")
        self.frozen= (self.frozen[0],True)
        self.update = types.MethodType(_update_functions[self.frozen], self)
    def unfreeze_inactivation(self):
        if self.frozen[1]==False:
            print("Warning: unfreezing inactivation that was not frozen!")
        self.frozen= (self.frozen[0],False)
        self.update = types.MethodType(_update_functions[self.frozen], self)
    def freeze_conductance(self):
        self.freeze_activation()
        self.freeze_inactivation()
    def unfreeze_conductance(self):
        self.unfreeze_activation()
        self.unfreeze_inactivation()

    def linearise(self,V = None, SI = False):
        """Calculates the RrL approximation of the channel"""
        if V == None:
            V = self.parent.V_m
        R = 1/self.g_inf(V)
        if self.frozen[0]:
            rp = 0
            Lp = 0
        else:
            rp = self.m_inf(V)/self.m_time(V)/self.m_order/self.g_inf(V)/(V-self.E)/( self.alpha_m_prime(V) - self.m_inf(V)*(self.alpha_m_prime(V) + self.beta_m_prime(V)) )
            Lp = rp*self.m_time(V)

        if self.h.size > 0 and self.frozen[1]==False:
            rq = 1/self.h_time(V)/self.inactivation_fraction_mode/self.h_order/self.g_max/self.m_inf(V)**self.m_order/self.h_inf(V)**(self.h_order-1)/(V-self.E)/( self.alpha_h_prime(V) - self.h_inf(V)*(self.alpha_h_prime(V) + self.beta_h_prime(V)) )
            Lq = rq*self.h_time(V)
        else:
            rq = []
            Lq = []

        if SI:
            return (R*1e3,rp*1e3,Lp,rq*1e3,Lq)
        return (R,rp,Lp,rq,Lq)

    def admitance(self,f,V=None):
        '''Inverse of the conductance impedance'''
        omega_ms = f*2*pi/1000 #Hz -> ms-1
        R,r_m,L_m,r_q,L_q = self.linearise(V)
        g_array_K = 1/R
        if abs(r_m) + abs(L_m) > 0:
            g_array_K += 1/(r_m+1j*omega_ms*L_m)
        for ii,x in enumerate(r_q):
            g_array_K += 1/(r_q[ii]+1j*omega_ms*L_q[ii])
        return g_array_K


    def current_variance(self,f,V=None):
        """Calculates the current variance due to voltage-gated channels in uA2/Hz"""
        if V == None:
            V = self.parent.V_m
        total_current_variance = zeros_like(f)
        N_channels = self.g_max/self.g_single_channel
        lorentzian = lambda coeff : 2*coeff/ (coeff**2 + (2*pi*f/1e3)**2) # f is Hz, so it converts it to ms-1
        ## first the non-inactivating part
        a = self.m_order     #Order of activation
        m = self.m_inf(V)    #probability activation gating at rest
        t = self.m_time(V)   #ms
        m_ = (1-m)/m
        suma = zeros(len(f))
        for ii in range(1,a+1):
            suma = suma + misc.comb(a,ii,exact=True)*m_**ii * lorentzian(ii/t)
        if self.h.size == 0:
            total_current_variance += 1e-3*suma*((m**a)*(V-self.E)*self.g_max)**2 / N_channels #All noise comes from activation, uA2*ms -> uA2/Hz
        else: ## then the rest
            total_current_variance += (1-sum(self.inactivation_fraction_mode))*1e-3*suma*((m**a)*(V-self.E)*self.g_max)**2 / N_channels #Non inactivation fraction
            print (self.h_order)
            for i,h_order in enumerate(self.h_order):
                print (h_order)
                if h_order == 1:
                    h = self.h_inf(V)[i]    #probability inactivation gating at rest
                    t_h = self.h_time(V)[i] #ms
                    h_ = (1-h)/h
                    suma = h_ * lorentzian(1/t_h)
                    for ii in range(1,a+1):
                        suma = suma + misc.comb(a,ii,exact=True)*m_**ii *( lorentzian(ii/t) + h_ * lorentzian(1/t_h+ii/t) )
                    total_current_variance += self.inactivation_fraction_mode[i]*1e-3*suma*((m**a)*h*(V-self.E)*self.g_max)**2 / N_channels #uA2*ms -> uA2/Hz
                else:
                    print ("Inactivation with order different than one non implemented. Sorry!")
        return total_current_variance


class DelayedRectifier(VoltageDependent):
    '''VoltageDependent without any kind of inactivation'''
    def update(self,V,dt):
        self.update_m(V,dt)
    def update_h(self,V,dt):
        return 0 #no inactivation
    def h_inf(self,V):
        return array([])
    def h_time(self,V):
        return array([])
    inactivation_modes = 0
