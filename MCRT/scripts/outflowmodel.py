'''
Functions taken from Spitzer-like solution with an external medium
WARNING - THE POWER LAW MODEL IN THIS DOES NOT USE OUTFLOW
          I HACKED THIS CODE HEAVILY, USE WITH CAUTION!!!
Sam Geen, February 2015
'''

from pymses.utils import constants as C
import profilesphere, rayprof

import numpy as np
import edgeradius
import testprofilefit

# NOTE!!!! THERE IS SOME HARD-CODING IN THIS, BE SUPER CAREFUL

# Some global properties
# BE VERY CAREFUL TO MAKE SURE THAT THESE ARE GOOOOOOOOD!
# Embedded -> star at centre of cloud (embedded = True)
# Blister -> star at edge of cloud (embedded = False)
tstart = 1.25
Te = 8400.0 # K
kB = 1.3806e-16 # erg / K
gamma = 1.4
X = 0.74
#mu = X*2.0 + 0.25*(1-X)*2.0 # Ionised hydrogen plus once-ionised He
mu = 0.61 # From Matzner 2002
mH = 1.67e-24 # g
cs = np.sqrt(gamma * kB * Te / (mH*mu))
print "USING cs = ", cs/1e5, "km/s"
def alpha_B_HII(T):
    # input  : T in K
    # output : HII recombination rate (in cm3 / s)
    l = 315614./T
    a = 2.753e-14 * l**1.5 / (1. + (l/2.74)**0.407)**2.242
    return a
beta2 = alpha_B_HII(Te)
#beta2 = 2e-10 * Te**(-0.75) # cm^3/s
G = 6.674e-8
pcincm = 3.08567758e18
Myrins = 3.15569e13
#profilemodule = rayprof
profilemodule = profilesphere

accretion_on = False
powerlaw = False

momback = 3.56571525e+42 # Initial momentum of cloud

def Findrstromgren(sim,dens=None):
    flux = FindFlux(sim)
    if dens == None:
        dens = FindnH(sim)
    Rs = (3.0/4.0/np.pi * flux /  dens**2 / beta2)**(1.0/3.0)
    print "Stromgren radius found of",Rs,"cm, (dens=",dens,"cm^-3)"
    return Rs

def Findtff(sim):
    dens = FindnH(sim)*mH/X
    print "Density", dens
    tff = np.sqrt(3.0/(2*np.pi*G*dens))
    print "Freefall time in seconds: ", tff
    # HACK !!
    tff = tstart*Myrins
    return tff
    
def unittimeMyr(sim):
    return sim.Snapshots()[0].RawData().info["unit_time"].express(C.Myr)

def FindNH22(sim):
    # Calculate the column density at t=1.25Myr the start of the sim
    unit_time_Myr = unittimeMyr(sim)
    start = sim.FindAtTime(tstart/unit_time_Myr)
    # 0.25 is the radius of the dense sphere
    NH22 = profilemodule.column_density(start,rcut=0.25)/1e22
    print "Using a value of N_H22 = ",NH22
    return NH22

def FindnH(sim,rcut=None):
    # Calculate the central density at t=1.25Myr the start of the sim
    unit_time_Myr = unittimeMyr(sim)
    start = sim.FindAtTime(tstart/unit_time_Myr)
    # 0.25 is the radius of the dense sphere
    #nH = profilemodule.mean_density(start)
    if rcut is not None:
        nH = profilemodule.central_density(start,rcut=rcut)
    else:
        nH = profilemodule.central_density(start)
    print "Using a value of n_H = ",nH
    return nH

def FindScaleRadius(sim,dens=None):
    #Find the scale radius r0 for rho = rho0 * (r/r0)**-2
    unit_time_Myr = unittimeMyr(sim)
    start = sim.FindAtTime(tstart/unit_time_Myr)
    # 0.25 is the radius of the dense sphere
    power, factor = profilemodule.powerlaw(start,"rho",rinner=0.1)
    factor = np.exp(factor)
    if dens is None:
        rho0 = profilemodule.central_density(start,rcut=0.1)
    else:
        rho0 = dens
    print "RHO0, POWER, FACTOR", rho0, power, factor
    r0 = (rho0 / factor)**(-1.0/power)
    print "Using a scale radius of", r0, "pc"
    return r0*pcincm

def FindProfilePower(sim):
    '''
    Find the power exponent of the density profile
    '''
    unit_time_Myr = unittimeMyr(sim)
    start = sim.FindAtTime(tstart/unit_time_Myr)
    # 0.25 is the radius of the dense sphere
    power, factor = profilemodule.powerlaw(start,"rho",rinner=0.1)
    return power

def FindFlux(sim):
    '''
    Find the flux. This uses the naming convention NXX_MY_BZZ
    '''
    if "noturb" in sim.Name():
        fluxstr = sim.Name()[6:8]
    else:
        fluxstr = sim.Name()[1:3]
    return 10.0**float(fluxstr)

def FindMass(sim):
    '''
    Find the cloud mass. This uses the naming convention NXX_MY_BZZ
    '''
    if noturb in sim.Name():
        return 1e4
    massstr = sim.Name()[5:6]
    return 10.0**float(massstr)

def FindriiOLD(sim):
    '''
    Find extent of HII region in pc with simulation snapshot time
    '''
    print "Finding Hennebelle radius for sim ", sim.Name()
    Myrins = 3.15569e13
    unit_time_Myr = unittimeMyr(sim)
    times = sim.Times()*unit_time_Myr - tstart
    times = times[times >= 0.0]
    times *= Myrins
    r0 = Findrstromgren(sim)
    cii = cs
    flux = FindFlux(sim)
    tff = Findtff(sim)
    print "r0, cii, times, tff", r0, cii, times, tff
    newt = 1.0 - 0.5*times / tff
    newt[newt < 0.0] = 0.0
    rii = r0 * (1.0 + 7.0/4.0 * cii / r0 * \
                times * newt)**(4.0/7.0)

    return times/Myrins, rii/pcincm

class Radius(object):
    '''
    Find the radius of the HII region
    '''

    Myrins = 3.15569e13
    rmin = 3.0*pcincm
    rmax = 13.0*pcincm

    def __init__(self,sim):
        self.sim = sim
        self.accretion = False
        self.powerlaw = False
        self.r0 = None
        self.cii = cs
        self.flux = FindFlux(sim)
        self.tff = Findtff(sim)
        self.n0 = None #FindnH(sim)
        self.nw = 1.0 # Force external medium to have a density of 1 at/cc
        self.powr0 = None
        self.pown = None
        unit_time_Myr = unittimeMyr(sim)
        times = sim.Times()*unit_time_Myr - tstart
        times = times[times >= 0.0]
        times *= self.Myrins
        self.times = times

    def _TCloudAccretion(self,r_cloud):
        # Find tcloud for the accretion solution
        # Find the components of the quadratic formula for the accretion soln
        a = 0.5/self.tff
        b = -1.0
        c = ( (r_cloud/self.r0)**(7.0/4.0) - 1.0 ) * \
            (4.0/7.0 * self.r0 / self.cii)
        det = b*b - 4.0*a*c
        if det < 0.0:
            # Front never escapes cloud if the sqrt part is -ve
            return None
        else:
            # Solve the quadratic formula
            # Only -sqrt(det) - we want the first time it reaches r_cloud
            return (-b - np.sqrt(det) ) / (2.0*a)

    def _TCloudSpitzer(self,r_cloud):
        # Time Spitzer solution reaches r_cloud
        t_cloud = 4.0/7.0 * self.r0 / self.cii * \
            ((r_cloud/self.r0)**(7.0/4.0) - 1.0)
        return t_cloud

    def _TCloudPowerLaw(self,r_cloud):
        # Time power law solution reaches cloud edge
        index = (7.0 - 2.0 * self.pown)/4.0
        fact = index * self.powr0 ** (-0.5*self.pown) * \
            self.r0 ** 0.75 * self.cii
        t_cloud = (r_cloud / fact) ** index
        return t_cloud

    def _FindriiRay(self,r_cloud):
        # Find the time at which the Spitzer solution reaches r_cloud
        #t_cloud = 4.0/7.0 * self.r0 / self.cii * \
        #    ((r_cloud/self.r0)**(7.0/4.0) - 1.0)
        t_cloud = None
        # Do we have an r_cloud?
        if not r_cloud is None:
            if self.accretion:
                t_cloud = self._TCloudAccretion(r_cloud)
            elif self.powerlaw:
                t_cloud = self._TCloudPowerLaw(r_cloud)
            else:
                t_cloud = self._TCloudSpitzer(r_cloud)
        if not t_cloud is None:
            # Split times into times before and after t_cloud
            oldtimes = self.times[self.times <= t_cloud]
            newtimes = self.times[self.times > t_cloud] - t_cloud
        else:
            # None means that the solution never escapes the cloud
            oldtimes = self.times
            newtimes = np.array([])
        #rold = self._FindriiSpitzer(oldtimes)
        if self.accretion and self.powerlaw:
            rold = self._FindriiPowerAccretion(oldtimes)
        elif self.accretion:
            rold = self._FindriiAccretion(oldtimes)
        elif self.powerlaw:
            rold = self._FindriiPowerLaw(oldtimes)
        else:
            rold = self._FindriiSpitzer(oldtimes)
        if not t_cloud is None:
            rnew = self._FindriiCloud(newtimes,r_cloud)
            radii = np.concatenate((rold,rnew))
        else:
            radii = rold
        return radii

    def _FindriiSpitzer(self,times):
        rii = self.r0 * (1.0 + 7.0/4.0 * self.cii / self.r0 * \
                        times)**(4.0/7.0)
        return rii

    def _FindriiAccretion(self,times):
        newt = 1.0 - 0.5*times / self.tff
        newt[newt < 0.0] = 0.0
        rii = self.r0 * (1.0 + 7.0/4.0 * self.cii / self.r0 * \
                        times * newt)**(4.0/7.0)
        return rii

    def _FindriiPowerLaw(self,times):
        index = (7.0 - 2.0 * self.pown)/4.0
        fact = index * self.powr0 ** (-0.5*self.pown) * \
            self.r0 ** 0.75 * self.cii
        rii = (fact * times)**(1.0/index)
        return rii

    def _FindriiPowerAccretion(self,times):
        newt = 1.0 - 0.5*times / self.tff
        newt[newt < 0.0] = 0.0
        index = (7.0 - 2.0 * self.pown)/4.0
        fact = index * self.powr0 ** (-0.5*self.pown) * \
            self.r0 ** 0.75 * self.cii
        rii = (fact * times * newt)**(1.0/index)
        return rii

    def _FindriiCloud(self,times,r_cloud):
        fact = (self.n0 / self.nw) ** 0.5 * (self.r0 / r_cloud)**(3.0/4.0)
        rii = r_cloud * (1.0 + 7.0/4.0 * fact * self.cii / r_cloud * \
                        times)**(4.0/7.0)
        return rii

    def FindScaleRadius(sim):
        ''' 
        Find the scale radius r0 for rho = rho0 * (r/r0)**-2
        '''
        unit_time_Myr = unittimeMyr(sim)
        start = sim.FindAtTime(tstart/unit_time_Myr)
        # 0.25 is the radius of the dense sphere
        power, factor = profilemodule.powerlaw(start,"rho")#,rinner=0.1)
        factor = np.exp(factor)
        rho0 = profilemodule.central_density(start)#,rcut=0.1)
        print "RHO0, POWER, FACTOR", rho0, power, factor
        r0 = (rho0 / factor)**(-1.0/power)
        print "Using a scale radius of", r0, "pc"
        return r0*pcincm

    def FindRadii(self):
        # Make the range of cloud radii
        #rng = self.rmax - self.rmin
        #rclouds = np.arange(self.rmin,self.rmax+1e-5,rng/100.0)
        self.r0 = Findrstromgren(self.sim)
        self.n0 = FindnH(self.sim)
        if self.powerlaw:
            # Power law ignores cloud cut-offs
            rhofit,n = testprofilefit.Run(self.sim.Name(),tstart)
            # Scale to r0, choose n0
            self.powr0 = np.exp((np.log(rhofit) - np.log(self.n0))/n)
            self.powr0 *= pcincm
            self.pown = n
            self.powrho0 = self.n0
            #rclouds *= 0.0 # HACK!
            #rclouds += 1e30 # HACK!
            radii = self._FindriiRay(None) # No r_cloud
        else:
            rclouds, dclouds = edgeradius.Run(self.sim.Name(),tstart,noplot=True)
            print "RCLOUDS", rclouds, "DCLOUDS", dclouds
            #rclouds = rclouds[rclouds > 0.0]
            # No detection? Just make really big
            rclouds[rclouds <= 0.0] = 1e9*rclouds.max()
            rclouds *= pcincm
            nrad = len(rclouds)
            # Average over all radii
            radii = np.zeros((len(self.times)))
            for i in range(0,nrad):
            #self.r0 = Findrstromgren(self.sim,dclouds[i])
            #self.n0 = dclouds[i]
                rcurr = self._FindriiRay(rclouds[i])
                radii += rcurr
            radii /= float(nrad)
        return self.times / Myrins, radii/pcincm

    def FindMomentum(self):
        '''
        NOTE THIS IS ONLY FOR POWER LAW
        '''
        if not self.powerlaw:
            print "MOMENTUM NEEDS POWER LAW FOR NOW!!!"
            raise ValueError
        self.r0 = Findrstromgren(self.sim)
        self.n0 = FindnH(self.sim)
        rhofit,n = testprofilefit.Run(self.sim.Name(),tstart)
        # Scale to r0, choose n0
        self.powr0 = np.exp((np.log(rhofit) - np.log(self.n0))/n)
        self.powr0 *= pcincm
        self.pown = n
        self.powrho0 = self.n0
        # Let alpha = 4 / (7-2*n)
        alpha = 4.0 / (7.0 - 2.0*self.pown)
        dummy, radii = self.FindRadii()
        radii *= pcincm
        rho0 = self.n0 * mH / X
        # Power law, so rdot = 1/alpha * r/t
        rdot = alpha * radii / self.times
        mass = 4 * np.pi / (3 - self.pown) * radii ** (3 - self.pown) * \
            self.r0 ** (self.pown) * rho0
        momentum = mass * rdot
        momentum += momback
        return self.times / Myrins, momentum

    def FindVelocity(self):
        '''
        NOTE THIS IS ONLY FOR POWER LAW
        '''
        if not self.powerlaw:
            print "VELOCITY NEEDS POWER LAW FOR NOW!!!"
            raise ValueError
        self.r0 = Findrstromgren(self.sim)
        self.n0 = FindnH(self.sim)
        rhofit,n = testprofilefit.Run(self.sim.Name(),tstart)
        # Scale to r0, choose n0
        self.powr0 = np.exp((np.log(rhofit) - np.log(self.n0))/n)
        self.powr0 *= pcincm
        self.pown = n
        self.powrho0 = self.n0
        # Let alpha = 4 / (7-2*n)
        alpha = 4.0 / (7.0 - 2.0*self.pown)
        dummy, radii = self.FindRadii()
        radii *= pcincm
        rho0 = self.n0 * mH / X
        # Power law, so rdot = 1/alpha * r/t
        rdot = alpha * radii / self.times
        return self.times / Myrins, rdot / 1e5 # km/s

    def FindMdest(self):
        '''
        NOTE THIS IS ONLY FOR POWER LAW
        '''
        if not self.powerlaw:
            print "MOMENTUM NEEDS POWER LAW FOR NOW!!!"
            raise ValueError
        self.r0 = Findrstromgren(self.sim)
        self.n0 = FindnH(self.sim)
        rhofit,n = testprofilefit.Run(self.sim.Name(),tstart)
        # Scale to r0, choose n0
        self.powr0 = np.exp((np.log(rhofit) - np.log(self.n0))/n)
        self.powr0 *= pcincm
        self.pown = n
        self.powrho0 = self.n0
        # Let alpha = 4 / (7-2*n)
        alpha = 4.0 / (7.0 - 2.0*self.pown)
        dummy, radii = self.FindRadii()
        radii *= pcincm
        rho0 = self.n0 * mH / X
        mdest = 4 * np.pi / 3.0 * (radii * self.r0) ** (1.5) * rho0
        Msolar = 1.9891e33
        return self.times / Myrins, mdest / Msolar

def Findrii(sim):
    print "Finding radius for sim ", sim.Name()
    rad = Radius(sim)
    rad.accretion = accretion_on
    rad.powerlaw = powerlaw
    return rad.FindRadii()

def FindriiAccretion(sim):
    print "Finding power-law density radius for sim ", sim.Name()
    rad = Radius(sim)
    rad.accretion = True
    rad.powerlaw = True # HACK powerlaw
    return rad.FindRadii()

def FindriiPowerLaw(sim):
    print "Finding power-law density radius for sim ", sim.Name()
    rad = Radius(sim)
    rad.accretion = accretion_on
    rad.powerlaw = True
    return rad.FindRadii()

def FindriiisothermalDEPRECATED(sim):
    '''
    Find extent of HII region in pc with simulation snapshot time
    Uses an isothermal density gradient
    '''
    print "Finding radius for sim ", sim.Name()
    Myrins = 3.15569e13
    unit_time_Myr = unittimeMyr(sim)
    times = sim.Times()*unit_time_Myr - tstart
    times = times[times >= 0.0]
    times *= Myrins
    rs = Findrstromgren(sim)
    r0 = Findrstromgren(sim)# FindScaleRadius(sim)
    cii = cs
    flux = FindFlux(sim)
    tff = Findtff(sim)
    print "r0, cii, times, tff", r0, cii, times, tff
    newt = 1.0 - 0.5*times / tff
    newt[newt < 0.0] = 0.0
    #rii = r0 * (3.0/4.0 * cii / r0 * \
    #            times * newt)**(4.0/3.0)
    rii = rs * (3.0/4.0 * cii / r0 * \
                    times * newt)**(4.0/3.0)
    # HACK
    print rii
    return times/Myrins, rii/pcincm

def FindriipowerprofileDEPRECATED(sim):
    '''
    Find extent of HII region in pc with simulation snapshot time
    Uses a power law density gradient
    '''
    print "Finding radius for sim ", sim.Name()
    Myrins = 3.15569e13
    unit_time_Myr = unittimeMyr(sim)
    times = sim.Times()*unit_time_Myr - tstart
    times = times[times >= 0.0]
    times *= Myrins
    rs = Findrstromgren(sim)
    r0 = FindScaleRadius(sim)
    cii = cs
    flux = FindFlux(sim)
    tff = Findtff(sim)
    power = FindProfilePower(sim)
    print "r0, cii, times, tff", r0, cii, times, tff
    newt = 1.0 - 0.5*times / tff
    newt[newt < 0.0] = 0.0
    # HACK
    newt*= 0.0
    newt += 1.0
    #rii = r0 * (3.0/4.0 * cii / r0 * \
    #            times * newt)**(4.0/3.0)
    newpower = (7.0 - 2.0*power)/4.0
    rii = (newpower * cii * r0**(-power/2.0) * rs ** (3.0/4.0) * \
                    times * newt)**(1.0/newpower)
    # HACK
    print rii
    return times/Myrins, rii/pcincm

def FindMom(sim):
    '''
    Find the momentum of the shell in g cm/s with simulation snapshot time
    '''
    print "Finding Matzner 2002 momentum for sim ", sim.Name()
    msolarkmps_cgs = 1.9891e38 # 1 Msolar km/s in g cm/s
    unit_time_Myr = unittimeMyr(sim)
    times = sim.Times()*unit_time_Myr - tstart
    times = times[times >= 0.0]
    NH22 = FindNH22(sim)
    mass = FindMass(sim)
    flux = FindFlux(sim)
    mom = momcoeff[embedded] * ( (times/3.7)**9 * (NH22/1.5)**(-1.5) * \
        (mass/1e6)**0.5 * (flux/1e49)**4 )**(1.0/7.0) * msolarkmps_cgs
    mom += 3.56571525e+42 # Initial momentum of cloud
    return times, mom

def FindMomPowerLaw(sim):
    print "Finding power-law density momentum for sim ", sim.Name()
    rad = Radius(sim)
    rad.accretion = accretion_on
    rad.powerlaw = True
    return rad.FindMomentum()

def FindRdotPowerLaw(sim):
    print "Finding power-law density shell velocity for sim ", sim.Name()
    rad = Radius(sim)
    rad.accretion = accretion_on
    rad.powerlaw = True
    return rad.FindVelocity()
    

def FindMdest(sim):
    '''
    Find the evaporated mass in msolar with simulation snapshot time
    Equation 19 plus the limit given in Equation 20
    NOTE/WARNING: Matzner 2002 only quotes for a blister region
    '''
    print "Finding Matzner 2002 ionised mass for sim ", sim.Name()
    unit_time_Myr = unittimeMyr(sim)
    times = sim.Times()*unit_time_Myr - tstart
    times = times[times >= 0.0]
    NH22 = FindNH22(sim)
    mass = FindMass(sim)
    flux = FindFlux(sim)
    mdest = mdestcoeff * ( (times/3.7)**9 * (NH22/1.5)**(-1.5) * \
        (mass/1e6)**0.5 * (flux/1e49)**4 )**(1.0/7.0)
    mmax = 4.6e4 * (NH22/1.5)**(-3.0/8.0) * (mass/1e6)**(7.0/8.0) * \
        (flux/1e49)**(1.0/4.0)
    mdest[mdest > mmax] = mmax
    return times, mdest

def FindMdestPowerLaw(sim):
    print "Finding power-law density ionised mass for sim ", sim.Name()
    rad = Radius(sim)
    rad.accretion = accretion_on
    rad.powerlaw = True
    return rad.FindMdest()

def FindMdestMax(sim):
    '''
    Equation 20, limit on ionised mass
    Returns a flat function in time at mdest,max
    '''
    unit_time_Myr = unittimeMyr(sim)
    times = sim.Times()*unit_time_Myr - tstart
    times = times[times >= 0.0]
    NH22 = FindNH22(sim)
    mass = FindMass(sim)
    flux = FindFlux(sim)
    mmax = 4.6e4 * (NH22/1.5)**(-3.0/8.0) * (mass/1e6)**(7.0/8.0) * \
        (flux/1e49)**(1.0/4.0)
    return times, times*0.0+mmax
    
if __name__=="__main__":

    import Hamu
    import matplotlib.pyplot as plt

    Hamu.Workspace("MCRT")
    
    sims = {}
    tstarts = [1.25,0.15625,0.15625]
    ts = {}
    '''
    for sim,tst in zip(["N48_M4_B02","N48_M4_B02_C2","N48_M4_B02_C"],tstarts):
        sims[sim] = Hamu.Simulation(sim)
        t = sims[sim].Times()
        ts[sim] = tst
    for s in sims.iterkeys():
        tstart = ts[s]
        Findtff(sims[s])
    '''
    s = Hamu.Simulation("N47_M4_B02")
    t,r = Findrii(s)
    #plt.plot(t,r)
    #plt.savefig("../plots/debugmodel.pdf")
    
