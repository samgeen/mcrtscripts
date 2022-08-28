''' 
Test the Pittard model for resolving wind injection
Sam Geen, April 2022
'''

from startup import *

def getmomentum(mstar,ageins):
    dtins = 1.0
    ewind, mwind = singlestar.star_winds(mstar,ageins,dtins)
    pwind = np.sqrt(2.0 * ewind * mwind)
    print("pwind for m=",mstar,"at",ageins,"s = ",pwind)
    def printforPamb(Pamb):
        rinjmax = np.sqrt(pwind / 4.0 / np.pi / Pamb)
        rinjmax /= pcincm
        print("rinjmax for Pamb=",Pamb," in pc = ",rinjmax)
    printforPamb(1e-9) # erg/cm^3, from gradrays plots in Fiducial run
    printforPamb(1e-8) # erg/cm^3, from gradrays plots in Duchamel run
    
if __name__=="__main__":
    getmomentum(35.0,3.16e12)
