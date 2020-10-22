
# Constants
X = 0.74
mH = 1.66e-24
G = 6.67428e-8
c = 2.99792458e10
gamma = 1.0 # isothermal gas
kB = 1.3806485279e-16 # in cgs  
yrins = 3.1556926e7
pc = 3.086e+18
Msun = 1.98847e33
def alpha_B_HII(T):
    # HII recombination rate
    # input  : T in K
    # output : HII recombination rate (in cm3 / s)
    l = 315614./T
    a = 2.753e-14 * l**1.5 / (1. + (l/2.74)**0.407)**2.242
    return a            

if __name__=="__main__":
    ci = 1e6
    Ti = ci**2 * (1.0/gamma)*mH/X/kB
    print("ALPHAB 10 km/s = ", alpha_B_HII(Ti))
    ci = 2e6
    Ti = ci**2 * (1.0/gamma)*mH/X/kB
    print("ALPHAB 20 km/s = ", alpha_B_HII(Ti))
