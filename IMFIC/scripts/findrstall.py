'''
Find rstall in the HIISN model
Sam Geen, July 2015
'''

from startup import *

def anarstall(Sstar,n0,r0,w,ci):
    '''
    Sstar - photon emission rate
    '''
    r0w = r0**w
    rho0 = n0 * mH / X
    Te = ci**2 * mH/X / (gamma * kB)
    alphaB = alpha_B_HII(Te)
    print "TE, ALPHAB", Te, alphaB
    vfact = 6.0/5.0
    rstall = ((3.0*Sstar / (4*np.pi * alphaB))**0.5 * \
              ((3.0-w) * ci**2 * mH/X) / (4*vfact * np.pi * G * rho0**2 * r0w**2))**(2.0/(7.0-4.0*w))
    return rstall/pcincm
    
if __name__=="__main__":
    # Read by eye on figure around first sink to make a star
    w = 1.13062013244 
    r0 = 2.49733076359 
    #w = 1.34532281129 
    #r0 = 1.63376570434
    #w = 1.38117192514 
    #r0 = 1.52347112405
    ci = 10 * 1e5 # 10 km/s in cm/s
    r0 *= pcincm
    n0 = 1e2
    for Sstar in 10.0**np.arange(47,51.00001,0.3):
        print Sstar, anarstall(Sstar,n0,r0,w,ci)
