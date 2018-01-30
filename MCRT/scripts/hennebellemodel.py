'''
Functions taken from Spitzer-like solution derived by Patrick Hennebelle 2015
Sam Geen, January 2015
'''

from pymses.utils import constants as C
import profilesphere, rayprof

import numpy as np

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
mu = X*2.0 + 0.25*(1-X)*2.0 # Ionised hydrogen plus once-ionised He
mH = 1.67e-24 # g
cs = np.sqrt(gamma * kB * Te / (mH/mu))
beta2 = 2e-10 * Te**(-0.75) # cm^3/s
G = 6.674e-8
pcincm = 3.08567758e18

profilemodule = rayprof

def Findrstromgren(sim):
    flux = FindFlux(sim)
    dens = FindnH(sim)
    Rs = (3.0/4.0/np.pi * flux /  dens**2 / beta2)**(1.0/3.0)
    print "Stromgren radius found of",Rs,"cm"
    return Rs

def Findtff(sim):
    dens = FindnH(sim)*mH/X
    print "Density", dens
    tff = np.sqrt(3.0/(2*np.pi*G*dens))
    print "Freefall time in seconds: ", tff
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

def FindnH(sim):
    # Calculate the central density at t=1.25Myr the start of the sim
    unit_time_Myr = unittimeMyr(sim)
    start = sim.FindAtTime(tstart/unit_time_Myr)
    # 0.25 is the radius of the dense sphere
    nH = profilemodule.central_density(start)
    print "Using a value of n_H = ",nH
    return nH

def FindScaleRadius(sim):
    '''
    Find the scale radius r0 for rho = rho0 * (r/r0)**-2
    '''
    unit_time_Myr = unittimeMyr(sim)
    start = sim.FindAtTime(tstart/unit_time_Myr)
    # 0.25 is the radius of the dense sphere
    power, factor = profilemodule.powerlaw(start,"rho",rinner=0.1)
    factor = np.exp(factor)
    rho0 = profilemodule.central_density(start,rcut=0.1)
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
    fluxstr = sim.Name()[1:3]
    return 10.0**float(fluxstr)

def FindMass(sim):
    '''
    Find the cloud mass. This uses the naming convention NXX_MY_BZZ
    '''
    massstr = sim.Name()[5:6]
    return 10.0**float(massstr)

def Findrii(sim):
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

def Findriiisothermal(sim):
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

def Findriipowerprofile(sim):
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
    s = Hamu.Simulation("N47_M4_B02")
    t,r = Findriiisothermal(s)
    plt.plot(t,r)
    plt.savefig("../plots/debugmodel.pdf")
    
