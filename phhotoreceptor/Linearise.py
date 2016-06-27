from numpy import array,zeros,ones_like,convolve,zeros_like,real, imag, pi,roots
from pylab import Arrow
__author__ = 'Francisco J. H. Heras'

class RrLC:
    def __init__(self,R,r,L,C):
        if len(r) == len(L):
            self.C = C #uF
            self.R = R #kOhm
            self.r = r #kOhm
            self.L = L #H
            self.n = len(r)+1 #number of resistor branches (including R)
        else:
            print("Error! r and L must have the same size")

    def impedance(self,f):
        omega_ms = f*2*pi/1000
        admittance = self.C*omega_ms*1j + 1/self.R #mS from light+leak and capacitance(uF)
        for branch in range(self.n-1):
            admittance += 1/(self.r[branch] + 1j*omega_ms*self.L[branch])
        return 1/admittance #kOhm

    def plot_admittance(self,f,ax,phenom_colour=['b','c']):
        '''Plot admittances at a frequency f (which has to be real and not array). Ax is the axes where it has to be plotted'''
        omega_ms = f*2*pi/1000
        scale = 1e6 #1e3 from mS -> pS
        admittance_C = scale*self.C*omega_ms*1j
        ax.plot([0,0],[0,imag(admittance_C)],'r',linewidth=4)
        G = scale*1/self.R #mS from light+leak and capacitance(uF)
        ax.plot([0,G],[0,0],'k',linewidth=4)
        admittance_ch = zeros(self.n-1,dtype=complex)
        admittance = G + admittance_C

        for branch in range(self.n-1):
            admittance_ch[branch] = scale*1/(self.r[branch] + 1j*omega_ms*self.L[branch])
            admittance += admittance_ch[branch]
            ax.plot([0,real(admittance_ch[branch])],[0,imag(admittance_ch[branch])],phenom_colour[branch],linewidth=4)

        #ax.axis([-abs(admittance),abs(admittance),-abs(admittance),abs(admittance)])
        ax.grid(True)


    def backup_plot_admittance(self,f,ax):
        '''Plot admittances at a frequency f (which has to be real and not array). Ax is the axes where it has to be plotted'''
        head_length=1e-8
        head_width=5e-9
        omega_ms = f*2*pi/1000
        admittance_C = self.C*omega_ms*1j
        ax.quiver(0, 0, 0, imag(admittance_C), head_width=head_width, head_length=head_length, fc='r', ec='r')
        G = 1/self.R #mS from light+leak and capacitance(uF)
        ax.arrow(0, 0, G, 0, head_width=head_width, head_length=head_length, fc='k', ec='k')
        admittance_ch = zeros(self.n-1,dtype=complex)
        admittance = G + admittance_C
        for branch in range(self.n-1):
            admittance_ch[branch] = 1/(self.r[branch] + 1j*omega_ms*self.L[branch])
            admittance += admittance_ch[branch]
            ax.arrow(0, 0, real(admittance_ch[branch]), imag(admittance_ch[branch]), head_width=head_width, head_length=head_length, fc='b', ec='b')

        ax.axis([-abs(admittance),abs(admittance),-abs(admittance),abs(admittance)])

    def print_poles(self):
        T = array(self.R)*array(self.C)/1e3 # ms -> Hz
        t = array(self.L)/array(self.r)/1e3 # ms -> Hz
        g = self.R*ones_like(self.r)/array(self.r)

        poly1=[1]
        for tt in t:
            poly1 = convolve(poly1,[tt,1])
        sumpoly2 = zeros(self.n+1)
        for i,tt in enumerate(t):
            poly2 = [0,0,g[i]]
            for ii,ttt in enumerate(t):
                if ii!=i:
                    poly2 = convolve(poly2,[ttt,1])
            sumpoly2 += poly2

        numerator = poly1
        denominator = sumpoly2 + convolve(poly1,[T,1])
        #print("zeros:", roots(numerator))
        #print("t:",-1/t)
        #print("max pole:", max(real(roots(denominator))))

        #if max(real(roots(denominator))) > 0:
        #    print ("Instability reached when taus = ", t*1e3, " ms")
        return max(real(roots(denominator)))

#def impedance(compartment):
#    '''Gives the impedance of a ElectricalCompartment'''
#    return (lambda f: compartment.impedance(f))

def give_RrLC_equivalent(compartment):
        C = compartment.C
        G = compartment.total_voltage_independent_conductance() #mS
        r=[]
        L=[]
        for channel in compartment.voltage_channels:
            Rc,rp,Lp,rq,Lq = channel.linearise()
            G+=1/Rc
            r.append(rp)
            L.append(Lp)
            r.extend(rq)
            L.extend(Lq)
        return RrLC(1/G,r,L,C)

def inject_current(photoreceptor,injected_current,dt):
    '''Generates a linearised photoreceptor and then solves the equations using forward Euler method'''


    linearised_ph = give_RrLC_equivalent(photoreceptor.body)

    C = linearised_ph.C
    R = linearised_ph.R
    r = linearised_ph.r
    L = linearised_ph.L
    n = linearised_ph.n

    A = zeros((n,n))
    b = zeros(n)


    A[0,0] = -1/R/C
    for i in range(1,n):
        A[0,i] = -1/C
        A[i,0] = 1/L[i-1]
        A[i,i] = -r[i-1]/L[i-1]
    b[0] = injected_current[0]/C

    V = zeros_like(injected_current)

    x = zeros(n)
    dx = zeros(n)
    for i in range(1,len(injected_current)):
        for ii in range(n):
            dx[ii] = (A[ii].dot(x) + b[ii])*dt
        x+=dx
        b[0] = injected_current[i]/C
        V[i] = x[0]


    return V

