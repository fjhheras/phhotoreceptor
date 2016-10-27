from numpy import sqrt, vectorize, linspace, angle, pi
import statistics
from scipy.optimize import minimize_scalar,bisect

__author__ = 'Francisco J. H. Heras'

#### TWO OPTIONS FOR GAIN AND BANDWIDTH
### Option 1; Gain defined as maximum value, Bandwidth defined as the frequency where Z falls to 1/Sqrt(2)
def Gain_Bandwidth_1(Z, f_cut_off_max = 5000, f_min = 0):
    Z_n = lambda f: -abs(Z(f))
    res = minimize_scalar(Z_n,bounds=(f_min,f_cut_off_max), method='bounded')
    f_max = res.x #gain taken at the maximum
    Z_s = lambda f: abs(Z(f)) + res.fun/sqrt(2)
    f_corner = bisect(Z_s,f_max,f_cut_off_max)
    return -res.fun,f_corner
def GBWP1(Z, f_min = 0):
    gain,bw = Gain_Bandwidth_1(Z, f_min=f_min)
    return gain*bw/1e3 #kOhm*Hz -> MOhm*Hz

### Option 2; Gain defined as DC value, Bandwidth defined as the frequency where Z falls to 1/Sqrt(2)
def Gain_Bandwidth_2(Z, f_cut_off_max = 5000, f_min = 0):
    Z_n = lambda f: -abs(Z(f))
    f_max = f_min #gain taken at the f-> 0 limit
    Z_s = lambda f: abs(Z(f)) + Z_n(f_max)/sqrt(2)
    f_cut_off = bisect(Z_s,f_max,f_cut_off_max)
    return abs(Z(f_max)),f_cut_off
def GBWP2(Z, f_min = 0):
    gain,bw = Gain_Bandwidth_2(Z, f_min = f_min)
    return gain*bw

### Choosing between them
which_gbwp = 1
if which_gbwp == 1:
    GBWP = GBWP1
    Gain_Bandwidth = Gain_Bandwidth_1
else:
    GBWP = GBWP2
    Gain_Bandwidth = Gain_Bandwidth_2


### UTILS

def Is_Band_Pass(Z,f_cut_off_max = 1000,df=0.1):
    '''Returns true if the function Z() is band_pass before f_cut_off'''
    Z_n = lambda f: -abs(Z(f))
    res = minimize_scalar(Z_n,bounds=(0,f_cut_off_max), method='bounded')
    return res.x > df

def GDD(Z,f=0,df=1e-4):
    '''Approximates Group Delay Dispersion, defined as the second derivative of Z() around the frequency f'''
    return (angle(Z(f+2*df)) - 2*angle(Z(f+df)) + angle(Z(f)))/(2*pi*df)**2  *1e6 #s2 -> ms2

def GDV(Z,f_down = 1,f_up = 100, n=250, df=1e-4):
    '''Aproximates Group Delay Dispersion, defined as the standard deviation of Z from frequency f_down up to frequency f_up'''
    group_delay = vectorize(lambda f : (angle(Z(f+df)) - angle(Z(f)))/(2*pi*df))
    f = linspace(f_down,f_up,n)
    return statistics.stdev(group_delay(f)) * 1e3 #s ->ms