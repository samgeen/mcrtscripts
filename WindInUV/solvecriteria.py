'''
Calculate point at which stars should overrun their photoionisation bubble with winds/radiation pressure
Sam Geen, December 2017 (lol, try August 2018)
'''

import sys

import customplot

import matplotlib.pyplot as plt

import numpy as np
import scipy.interpolate
import ionisedtemperatures

from consts import *

import stars
from stars import _Star # needed to make caching work

from collections import OrderedDict

# Stellar track stuff
windtracks = None
trackloc = "../StellarSources/"

Msolar = "M$_{\odot}$"

clabelfmt = '%1.1f'

# Colour maps for the various variables
cmaps = {}
cmaps["Rstall"] = "BrBG"
cmaps["Cb"] = "PiYG_r"
cmaps["Cw"] = "RdBu_r"
cmaps["Crp"] = "PuOr_r"

def SoundSpeed(cloud,star,t):
    #HACK - ignore change in Ti due to ionised gas density
    Ti = ionisedtemperatures.FindTemperature(star.Teff(t),cloud.metal)
    # Find ci
    # 2.0 factor is because of free electrons
    ci = np.sqrt(2.0 * gamma * kB * Ti * X / mH)
    return ci

def SoundSpeedFromT(temp):
    #HACK - ignore change in Ti due to ionised gas density
    Ti = temp
    # Find ci
    # 2.0 factor is because of free electrons
    ci = np.sqrt(2.0 * gamma * kB * Ti * X / mH)
    return ci

def FindalphaB(cloud,star,t):
    #HACK - ignore change in Ti due to ionised gas density
    Ti = ionisedtemperatures.FindTemperature(star.Teff(t),cloud.metal)
    return alpha_B_HII(Ti)

class IntegratorError(Exception):
    pass

def CloudFromr0n0(r0,n0,metal,accreting=True,gravityOn=True,ci=None):
    M0 = 4 * np.pi * n0 * r0**3.0 * mH/X
    Sigma0 = np.pi * n0 * r0 * mH / X
    return Cloud(M0,Sigma0,metal,accreting,gravityOn,ci)

class Cloud(object):
    def __init__(self,M0,Sigma0,metal,accreting=True,gravityOn=True,ci=None):
        self._M0 = M0
        self._Sigma0 = Sigma0
        self._metal = metal
        self._rho0 = 2.0/np.pi * (Sigma0**3 / M0)**0.5
        self._r0 = 0.5 * np.sqrt(M0/Sigma0) 
        self._tff0 = np.sqrt(3.0 * np.pi / (32.0 * G * self._rho0))
        vconst = 1.5
        if not accreting:
            vconst = 0.5
        self._v0 = np.sqrt(vconst * np.pi * G * self._rho0 * self._r0**2.0)
        self._accreting = accreting
        self._gravityOn = gravityOn
        self._ci = ci

    @property
    def M0(self):
        return self._M0

    @property
    def Sigma0(self):
        return self._Sigma0

    @property
    def metal(self):
        return self._metal

    @property
    def n0(self):
        return self._rho0*X/mH

    @property
    def r0(self):
        return self._r0

    def ci(self,star,t):
        if self._ci is None:
            return SoundSpeed(self,star,t)
        else:
            return self._ci

    @property
    def tff0(self):
        return self._tff0

    def alphaB(self,star,t):
        return FindalphaB(self,star,t)

    def v0(self,t):
        if self._gravityOn:
            return self._v0 * np.sqrt(self.AccretionFunc(t))
        else:
            return 0.0

    def AccretionFunc(self,t):
        return 1.0
        # NOTE: If Mdot is const, the density profile should be constant
        # The mass just flows from r0 to the star at r=0
        # if self._accreting:
        #     return t/self.tff0+1.0
        # else:
        #     return 1.0

    def MassInsideR(self,r):
        return self._M0 * r / self._r0

def TCool(t,star,cloud):
    Lw = star.WindLuminosity(t)
    L38 = Lw/1e38
    M0 = cloud.M0
    r0 = cloud.r0
    Sigma0 = cloud.Sigma0
    ci = cloud.ci(star,t)
    Ft = cloud.AccretionFunc(t)
    QH = star.NPhotons(t)
    alphaB = cloud.alphaB(star,t)
    ni = np.sqrt((3.0*QH/(alphaB*4*np.pi))/r0**3.0)
    tc = 16.0 * ni**(-8.0/11.0) * L38**(3.0/11.0) 
    return tc

def WindCriterion(t,star,cloud,cool=False):
    kw = np.sqrt(7.0/25.0) * (375.0*(gamma-1.0)/(28.0*np.pi*(9.0*gamma-4.0)))
    Lw = star.WindLuminosity(t)
    pw = star.WindMomentumRate(t)
    M0 = cloud.M0
    r0 = cloud.r0
    n0 = cloud.n0
    Sigma0 = cloud.Sigma0
    ci = cloud.ci(star,t)
    Ft = cloud.AccretionFunc(t)
    QH = star.NPhotons(t)
    alphaB = cloud.alphaB(star,t)
    #print kw, Lw, M0, Sigma0, ci, Ft, alphaB, QH
    f = 2.0 # CHANGE THIS!
    # Based on stalling
    #Cw = (kw*Lw*X*M0*Sigma0 * gamma/(mH*2*ci**2))**0.25 * (Ft/(16*np.pi**2*f*G))**0.5
    Cci = (ci**2)
    Cqh = (3.0*QH/(alphaB*4*np.pi)) 
    Clw = (kw*Lw*X/mH)
    Cpw = (pw*X/(4*np.pi*mH*n0*r0**2))
    Cn0 = (np.sqrt(M0*Sigma0)/(2.0*np.pi)*X/mH*Ft)
    # Based on ratio of effects
    Awadiabatic = Cci**(-0.75) * (Clw/Cn0)**(1.5)
    Awcool = (Cpw)**1.5
    if not cool:
        Aw = Awadiabatic
    else:
        Aw = Awcool
    Ai = Cqh * (Cci/Cn0)**2
    Cw = 2**0.25 * Aw*(Ai*r0)**(-0.75)
    # HACK - redo because getting weird answers
    #if cool:
    #    ri = r0
    #    Cw = 2.0**0.25 * (pw * X * gamma / (4.0*np.pi*mH*2.0*ci**2.0))**1.5 * (QH*3*ri/(alphaB*4*np.pi*ri))**(-0.75)
    #import pdb; pdb.set_trace()
    return Cw

def RadiationPressureCriterion(t,star,cloud,ftrap=2.0):
    OLDCODE = False
    WITHGRADIENT = False
    Sigma0 = cloud.Sigma0
    ci = cloud.ci(star,t)
    alphaB = cloud.alphaB(star,t)
    ephoton = star.EPhoton(t)
    Ft = cloud.AccretionFunc(t)
    r0 = cloud.r0
    QH = star.NPhotons(t)
    pw = star.WindMomentumRate(t)
    cifact = (ci**2 * mH/X)
    if OLDCODE:
        #Crp = Sigma0*alphaB*ephoton*gamma/(2*c) * (X / (mH*ci))**2 * Ft
        Crp = ephoton/c * (X /mH)* (ci**(-2)) * (3*QH*alphaB / (4 * np.pi * r0))**0.5
    elif WITHGRADIENT:
        # New method
        # Failure case if no emission at all
        if pw == 0.0 or QH == 0.0:
            return 0.0
        ri = r0
        rw = ri / RivsRw(t,star,cloud,ri)
        nouter = 1.0 / (cifact * (4.0*np.pi*rw**2/pw - ephoton/c * alphaB (ri-rw)))
        deltaPrp = QH * ephoton / (4 * np.pi * r0**2 * c) * ftrap
        P0 = cifact * nouter
        Crp = deltaPrp/P0
    else:
        # New method
        # Failure case if no emission at all
        if pw == 0.0 or QH == 0.0:
            return 0.0
        ri = r0
        rw = ri / RivsRw(t,star,cloud,ri,dustfrac=0.73)
        nouter = 1.0 / (cifact * (4.0*np.pi*rw**2/pw))
        deltaPrp = QH * ephoton / (4 * np.pi * r0**2 * c) * ftrap
        P0 = cifact * nouter
        Crp = deltaPrp/P0
    return Crp

def BreakoutCriterion(t,star,cloud):
    Sigma0 = cloud.Sigma0
    ci = cloud.ci(star,t)
    alphaB = cloud.alphaB(star,t)
    ephoton = star.EPhoton(t)
    Ft = cloud.AccretionFunc(t)
    r0 = cloud.r0
    M0 = cloud.M0
    QH = star.NPhotons(t)
    pw = star.WindMomentumRate(t)
    cifact = (ci**2 * mH/X)
    if QH == 0.0:
        return 0.0
    # Photon breakout condition from paper
    Cb = 0.8 * (M0/(100*Msun))**(-0.5) * (Sigma0 / (100*Msun /pc/pc))**(-1.5) * \
        (QH / 1e49) * (ci/1e6)**(-4.0)
    return Cb

def GravCriterion(t,star,cloud):
    M0 = cloud.M0/Msun
    Sigma0 = cloud.Sigma0/(Msun/pc/pc)
    ci = cloud.ci(star,t)/1e5 # km/s
    Ft = cloud.AccretionFunc(t)
    Cesc = np.sqrt(1.0 / ci**2) * 0.656 * (M0*Ft*Sigma0/1e4)**0.25
    return Cesc

def RivsRw(t,star,cloud,ri=None,dustfrac=1.0):
    # dustfract = fraction left after losses from dust (Krumholz & Matzner 09 use 0.73)
    kw = np.sqrt(7.0/25.0) * (375.0*(gamma-1.0)/(28.0*np.pi*(9.0*gamma-4.0)))
    Lw = star.WindLuminosity(t)
    pw = star.WindMomentumRate(t)
    M0 = cloud.M0
    r0 = cloud.r0
    n0 = cloud.n0
    Sigma0 = cloud.Sigma0
    ci = cloud.ci(star,t)
    Ft = cloud.AccretionFunc(t)
    QH = star.NPhotons(t)*dustfrac
    alphaB = cloud.alphaB(star,t)
    # Equation is ri**4 = const * ri * rw**4 + ri * rw**3
    # Fix ri = r0
    # So rw**4 (ri * const) + rw**3 (ri) + rw**0 (-ri**4) = 0
    const = 3.0 * QH / (4 * np.pi * alphaB) * \
        ((ci**2 * mH * 4 * np.pi) / (X * pw))**2.0
    if ri is None:
        ri = r0
    if ri <= 0:
        return 0.0
    if np.isnan(ri):
        return 0.0
    # Solve the equation
    print ri
    roots = np.roots([ri * const,ri,0.0,0.0,-ri**4])
    if roots is None:
        return 0.0
    if roots.max() <= 0.0:
        return 0.0
    # This will be a 4-element array with complex numbers in it
    # rw is the one positive, wholly real value in this array
    roots = roots[np.isreal(roots)]
    root = roots[roots > 0.0][0]
    rw = root.real
    # Return the ratio
    return ri / rw

def Findni(t,star,cloud,ri=None,nowind=False):
    kw = np.sqrt(7.0/25.0) * (375.0*(gamma-1.0)/(28.0*np.pi*(9.0*gamma-4.0)))
    Lw = star.WindLuminosity(t)
    pw = star.WindMomentumRate(t)
    M0 = cloud.M0
    r0 = cloud.r0
    n0 = cloud.n0
    Sigma0 = cloud.Sigma0
    ci = cloud.ci(star,t)
    Ft = cloud.AccretionFunc(t)
    QH = star.NPhotons(t)
    alphaB = cloud.alphaB(star,t)
    # Let ni**0.5 = nih
    # Equation is nih**4*ri**3 - constw * nih * rw**4 - consti = 0
    # Fix ri = r0
    # So rw**4 (ri * const) + rw**3 (ri) + rw**0 (-ri**4) = 0
    constw = 0.0
    if not nowind:
        constw = (pw * X / (4.0*np.pi * mH * ci**2.0))**1.5
    consti = QH * 3.0 / (alphaB * 4 * np.pi)
    if ri is None:
        ri = r0
    # Solve the equation
    roots = np.roots([ri**3.0,0.0,0.0,-constw,-consti])
    # This will be a 4-element array with complex numbers in it
    # rw is the one positive, wholly real value in this array
    roots = roots[np.isreal(roots)]
    root = roots[roots > 0.0][0]
    ni = np.sqrt(root.real)
    # Return the ratio
    return ni

def rstall_iononly(t,star,cloud):
    # Ionisation front stall radius
    r0 = cloud.r0
    n0 = cloud.n0
    ci = cloud.ci(star,t)
    Ft = cloud.AccretionFunc(t)
    QH = star.NPhotons(t)
    alphaB = cloud.alphaB(star,t)
    rstall = (8.0*np.pi*G * mH/X)**2.0 
    rstall *= (Ft*n0*r0**2/ci)**4.0
    rstall /= (3.0 * QH / 4.0 / np.pi / alphaB)
    return rstall

def rstall(t,star,cloud,cool=False):
    kw = np.sqrt(7.0/25.0) * (375.0*(gamma-1.0)/(28.0*np.pi*(9.0*gamma-4.0)))
    Lw = star.WindLuminosity(t)
    pw = star.WindMomentumRate(t)
    M0 = cloud.M0
    r0 = cloud.r0
    n0 = cloud.n0
    Sigma0 = cloud.Sigma0
    ci = cloud.ci(star,t)
    Ft = cloud.AccretionFunc(t)
    QH = star.NPhotons(t)
    alphaB = cloud.alphaB(star,t)
    v0 = cloud.v0(t)
    Cci = (ci**2)
    Cqh = (3.0*QH/(alphaB*4*np.pi)) 
    Clw = (kw*Lw*X/mH)
    Cpw = (pw*X/(4*np.pi*mH*n0*r0**2))
    Cn0 = (np.sqrt(M0*Sigma0)/(2.0*np.pi)*X/mH*Ft)
    # Based on ratio of effects
    Awadiabatic = Cci**(-0.75) * (Clw/Cn0)**(1.5)
    Awcool = (Cpw)**1.5
    if not cool:
        Aw = Awadiabatic
    else:
        Aw = Awcool
    Ai = Cqh * (Cci/Cn0)**2
    Cw = 2**0.25 * Aw*(Ai*r0)**(-0.75)
    rstall = (v0**4 - Aw * v0) / Ai
    if rstall <= 0.0:
        rstall = 1e-10
    if np.isnan(rstall):
        rstall = 1e-10
    return rstall

def deltaPrpoverP(t,star,cloud,ri,nowind=False,ftrap=2.0):
    NOGRADIENT = True
    # Radiation pressure gradient in a flat profile divided by the wind pressure = thermal pressure
    pw = star.WindMomentumRate(t)
    M0 = cloud.M0
    r0 = cloud.r0
    n0 = cloud.n0
    Sigma0 = cloud.Sigma0
    ci = cloud.ci(star,t)
    Ft = cloud.AccretionFunc(t)
    v0 = cloud.v0(t)
    QH = star.NPhotons(t)
    alphaB = cloud.alphaB(star,t)
    ephoton = star.EPhoton(t)
    ri = r0
    rw = ri / RivsRw(t,star,cloud,ri)
    cifact = (ci**2 * mH/X)
    # Failure case if no emission at all
    if pw == 0.0 or QH == 0.0:
        return 0.0
    if NOGRADIENT:
        nouter = 1.0 / (cifact * (4.0*np.pi*rw**2/pw))
    else:
        nouter = 1.0 / (cifact * (4.0*np.pi*rw**2/pw - ephoton/c * alphaB * (ri-rw)))
    deltaPrp = QH * ephoton / (4 * np.pi * ri**2 * c) * ftrap
    P0 = cifact * nouter
    dPrpoverPwind = deltaPrp/P0
    dPrpoverPnowind = (alphaB * QH / (12.0 * np.pi * ri))**0.5 * ephoton/c / cifact * ftrap
    if not nowind:
        #import pdb; pdb.set_trace()
        dPrpoverP = dPrpoverPwind
    else:
        dPrpoverP = dPrpoverPnowind
    #import pdb; pdb.set_trace()
    # Solve
    #dPrpoverP = 3.0 * ephoton / c * QH * (ri-rw) / (ri**3 - rw**3) * rw**2 / pw
    return dPrpoverP


def drdt_wind(t,ri,star,cloud,cool,nowind):
    kw = np.sqrt(7.0/25.0) * (375.0*(gamma-1.0)/(28.0*np.pi*(9.0*gamma-4.0)))
    Lw = star.WindLuminosity(t)
    pw = star.WindMomentumRate(t)
    M0 = cloud.M0
    r0 = cloud.r0
    n0 = cloud.n0
    Sigma0 = cloud.Sigma0
    ci = cloud.ci(star,t)
    Ft = cloud.AccretionFunc(t)
    v0 = cloud.v0(t)
    QH = star.NPhotons(t)
    alphaB = cloud.alphaB(star,t)
    #print kw, Lw, M0, Sigma0, ci, Ft, alphaB, QH
    f = 2.0 # CHANGE THIS!
    # Based on stalling
    #Cw = (kw*Lw*X*M0*Sigma0 * gamma/(mH*2*ci**2))**0.25 * (Ft/(16*np.pi**2*f*G))**0.5
    Cci = (ci**2)
    Cqh = (3.0*QH/(alphaB*4*np.pi)) 
    Clw = (kw*Lw*X/mH)
    Cpw = (pw*X/(4*np.pi*mH*n0*r0**2))
    Cn0 = (np.sqrt(M0*Sigma0)/(2.0*np.pi)*X/mH*Ft)
    # Based on ratio of effects
    Awadiabatic = Cci**(-0.75) * (Clw/Cn0)**(1.5)
    Awcool = (Cpw)**1.5
    if not cool:
        Aw = Awadiabatic
    else:
        Aw = Awcool
    Ai = Cqh * (Cci/Cn0)**2
    # Solve equation dri/dt**4 - Aw * dri/dt - Ai * ri = 0
    if not nowind:
        # Roots are [1,0,0,-Aw,-Airi]
        # Solve the equation
        roots = np.roots([1.0,0.0,0.0,-Aw,-Ai*ri])
        #print roots
        # This will be a 4-element array with complex numbers in it
        # rw is the one positive, wholly real value in this array
        try:
            roots = roots[np.isreal(roots)]
            root = roots[roots > 0.0][0]
            drdt = root.real
        except IndexError:
            raise IntegratorError
    else:
        drdt = (Ai*ri)**0.25
    #print t, ri, drdt, v0
    return drdt - v0

def test():
    starmass = 30.0
    sfe = 0.03
    print "Reading star mass",starmass,"Msun..."
    star = stars.Star(starmass,0.014)
    cloud = Cloud(1000.0*Msun,100.0*Msun/pc**2,1e6)
    #print "lifetime/yr", star.Lifetime()/yrins
    print "Calculating...", 
    Cw = WindCriterion(1e6*yrins,star,cloud)
    print "Cw", Cw
    print "Crp", RadiationPressureCriterion(1e6*yrins,star,cloud)
    print "Cesc", GravCriterion(t,star,cloud)

def plotvsmass(cool=False):
    Cws = []
    Crps = []
    tcools = []
    masses = stars.starmasses
    metal = 0.014
    cloud = Cloud(100.0*Msun,100.0*Msun/pc**2,metal)
    for mass in masses:
        star = stars.Star(mass,metal)
        #lifetime = star.Lifetime()
        Cws.append(WindCriterion(1e6*yrins,star,cloud,cool=cool))
        Crps.append(RadiationPressureCriterion(1e6*yrins,star,cloud))
        tcools.append(TCool(1e6*yrins,star,cloud))
    print "Cws", Cws
    print "Crps", Crps
    # Mass test
    plt.clf()
    plt.plot(masses, Cws,"k-")
    plt.plot(masses, Crps,"k--")
    plt.plot([0,120],[1,1],"k:",alpha=0.5)
    plt.ylim([1e-3,1e3])
    plt.yscale("log")
    cooltxt = ""
    if cool:
        cooltxt = "_cool"
    plt.savefig("plots/masstest"+cooltxt+".pdf")
    # Cooling time
    plt.clf()
    plt.plot(masses, tcools,"k-")
    plt.yscale("log")
    plt.savefig("plots/cooltest.pdf")

def ionisedtemperatureversusmass(soundspeed=False):
    plt.clf()
    masses = stars.starmasses
    labels = [r"$Z=0.002$",r"$Z=0.014$"]
    lines = ["--","-"]
    metals = [0.002,0.014]
    t = 1e6*yrins
    for label, line, metal in zip(labels, lines, metals):
        Tis = []
        cloud = Cloud(100.0*Msun,100.0*Msun/pc**2,metal)
        for mass in masses:
            star = stars.Star(mass,metal)
            if soundspeed:
                Ti = SoundSpeed(cloud,star,t) / 1e5 # cm/s to km/s
            else:            
                Ti = ionisedtemperatures.FindTemperature(star.Teff(t),cloud.metal)
            Tis.append(Ti)
        plt.plot(masses, Tis,"k"+line,label=label)
    plt.legend(frameon=False,loc="lower right")
    plt.xlabel("Stellar Mass / "+Msolar)
    if soundspeed:
        plt.ylabel("$c_i$ / km/s")
    else:
        plt.ylabel("Ionised Gas Temperature / K")
    filename = "ionisedtemperatureversusmass.pdf"
    if soundspeed:
        filename = "ionisedsoundspeedversusmass.pdf"
    plt.savefig("plots/"+filename)

def plotsurfaces(surfdens=1e2, metal=0.014, cool=False):
    # Surface density
    surf = surfdens*Msun/pc**2
    # Stellar mass
    starmasses = stars.starmasses
    nstarmass = len(starmasses)
    # Cloud mass
    ncloudmass = 20
    cloudmasses = np.logspace(2,5,ncloudmass)*Msun
    #surfs = np.logspace(1e1,1e4,nsurf)*Msun/pc**2
    # Other setup stuff
    #ci = 1e6
    #if metal == 0.002:
    #    ci = 2e6
    metaltxt = str(metal).replace(".","p")
    t = 1e6*yrins
    nbinary = 2
    Cws  = np.zeros((ncloudmass, nstarmass,nbinary))
    Crps = np.zeros((ncloudmass, nstarmass,nbinary))
    Cbs = np.zeros((ncloudmass, nstarmass,nbinary))
    Rstalls = np.zeros((ncloudmass, nstarmass,nbinary))
    # Escape velocity stuff
    testcloud = Cloud(100*Msun,surf,metal)
    #Cesc = GravCriterion(t,None,testcloud) # Doesn't use any of the star properties
    #Mesc = 100.0*Msun / Cesc**4
    # Iterate over cloud masses
    for icloud in range(0,ncloudmass):
        cloudmass = cloudmasses[icloud]
        cloud = Cloud(cloudmass,surf,metal)
        # Iterate over stellar masses
        for istar in range(0,nstarmass):
            starmass = starmasses[istar]
            star = stars.Star(starmass,metal)
            binary = stars.Binary(starmass,metal)
            # Single star
            Cw = WindCriterion(t,star,cloud,cool=cool)
            Crp = RadiationPressureCriterion(t,star,cloud)
            Cb = BreakoutCriterion(t,star,cloud)
            Rstall = rstall(t,star,cloud,cool=cool) / pc # MAKE IN PC, NOT IN CM
            Cws[icloud,istar,0] = Cw
            Crps[icloud,istar,0] = Crp
            Cbs[icloud,istar,0] = Cb
            Rstalls[icloud,istar,0] = Rstall
            # Binaries
            Cw = WindCriterion(t,binary,cloud,cool=cool)
            Crp = RadiationPressureCriterion(t,binary,cloud)
            Cb = BreakoutCriterion(t,binary,cloud)
            Rstall = rstall(t,binary,cloud,cool=cool) / pc # MAKE IN PC, NOT IN CM
            Cws[icloud,istar,1] = Cw
            Crps[icloud,istar,1] = Crp
            Cbs[icloud,istar,1] = Cb
            Rstalls[icloud,istar,1] = Rstall
    for data, name in zip([Cws,Crps,Cbs,Rstalls],["Cw","Crp","Cb","Rstall"]):
        fontsize = 12
        cooltxt = ""
        if cool:
            cooltxt = "_cool"
        plotname = "plots/surface_"+name+"_Sigma"+str(int(surfdens))+cooltxt+"_Z"+metaltxt+".pdf"
        print "Making figure", plotname
        plt.clf()
        X, Y = np.meshgrid(starmasses, np.log10(cloudmasses/Msun))
        Z = data[:,:,0]
        Zbin = data[:,:,1]
        print Zbin.max(), Z.max()
        #import pdb; pdb.set_trace()
        V = [0.1,1.0,10.0]
        contours = plt.contour(X,Y,Z,V,colors="black", zorder=2)
        contourbin = plt.contour(X,Y,Zbin,V,colors="black",linestyles="dashed", zorder=10)
        # Vesc line
        # (commented out because rstall is a better measure)
        #if Mesc < cloudmasses.max():
        #    plt.plot(starmasses,np.log10(Mesc/Msun)*np.ones(nstarmass),"k:")
        #    mesctxt = plt.text(starmasses[9*nstarmass//10],np.log10(Mesc/Msun)+0.05,r"$C_{esc}>1$", fontsize=fontsize)
        # Plot SFE=0.1 line
        plt.plot(starmasses,np.log10(10.0*starmasses),"k:")
        # Hack to make the SFE labels play nice
        sfetxtpos = 7*nstarmass//10
        if surfdens==1e2 and metal==0.014:
            sfetxtpos = 5*nstarmass//10
        sfetxt = plt.text(starmasses[sfetxtpos],\
                            np.log10(10.0*starmasses[sfetxtpos])+0.1,r"SFE=0.1", fontsize=fontsize)
        surfstr = "$10^{"+str(int(np.log10(surfdens)))+"}$"
        plt.title("$\Sigma_0$ = "+surfstr+" M$_{\odot}$/pc$^2$, Z="+metaltxt.replace("p","."))
        plt.legend(frameon=False)
        plt.clabel(contours, inline=True, fontsize=fontsize, fmt=clabelfmt)
        plt.clabel(contourbin, inline=True, fontsize=fontsize, fmt=clabelfmt)
        Xlow = np.log10(cloudmasses.min()/Msun)
        Xhigh = np.log10(cloudmasses.max()/Msun)
        Ylow = starmasses.min()
        Yhigh = starmasses.max()
        extents = np.array([Ylow,Yhigh,Xlow,Xhigh])
        plt.imshow(np.log(Z), extent=extents, origin='lower',cmap=cmaps[name], 
                    alpha=0.5, aspect="auto",vmin=-2,vmax=2,zorder=2)
        plt.xlabel("Stellar mass / M$_{\odot}$")
        plt.ylabel("Cloud mass / M$_{\odot}$")
        # Note: renamed stall to launch to more accurately describe the behaviour
        labels = {"Cw":"$C_{w}$","Crp":"$C_{rp}$","Cb":"$C_{B}$","Rstall":"$r_{launch}$ / pc"}
        plt.colorbar(label="log("+labels[name]+")")
        plt.savefig(plotname)

def plotvsclustermass(surfdens=1e2, sfe=0.1, metal=0.014, cool=False,legend=False):
    xmasses = []
    # Just a big ol scatter list
    Cws = []
    Crps = []
    Rstalls = []
    Cbs = []
    tcools = []
    # List of cloud masses
    logmmin = 2
    logmmax = 7
    cloudmasses = np.logspace(logmmin,logmmax,2*(logmmax - logmmin)+1)
    ncloudmasses = len(cloudmasses)
    # There and back again, used for surface plots
    taba = np.concatenate((cloudmasses,cloudmasses[::-1]))
    #ci = 1e6
    #if metal == 0.002:
    #    ci = 2e6
    metaltxt = str(metal).replace(".","p")
    t = 1e6*yrins
    # Make quartiles
    QCw = OrderedDict()
    QCrp = OrderedDict()
    QRstall = OrderedDict()
    QCb = OrderedDict()
    # Put median first so the label on the plot is the solid line
    quartiles = [49,0,24,74,99]
    lines = ["-",":","--","--",":"]
    for q in quartiles:
        QCw[q] = []
        QCrp[q] = []
        QRstall[q] = []
        QCb[q] = []
    # Every percentile
    Crpcontour = np.zeros((2*ncloudmasses,50))
    Cwcontour = np.zeros((2*ncloudmasses,50))
    Rstallcontour = np.zeros((2*ncloudmasses,50))
    Cbcontour = np.zeros((2*ncloudmasses,50))
    # Escape velocity stuff (DEPRECATED, RSTALL IS BETTER)
    testcloud = Cloud(100*Msun,surfdens*Msun/pc**2,metal)
    #Cesc = GravCriterion(t,None,testcloud) # Doesn't use any of the star properties
    #Mesc = 100.0 / Cesc**4
    # Run through cloud masses
    print "Computing clusters for cloud mass"
    for ic, cloudmass in enumerate(cloudmasses):
        print cloudmass,"..."
        cloud = Cloud(cloudmass*Msun,surfdens*Msun/pc**2,metal)
        Cwmass = []
        Crpmass = []
        Rstallmass = []
        Cbmass = []
        for seed in range(0,100):
            xmasses.append(cloudmass)
            star = stars.ClusterOnTheFly(sfe*cloudmass,metal,seed)
            Cw = WindCriterion(t,star,cloud,cool=cool)
            Crp = RadiationPressureCriterion(t,star,cloud)
            Rstall = rstall(t,star,cloud,cool=cool) / pc # MAKE IN PC, NOT IN CM
            Cb = BreakoutCriterion(t,star,cloud)
            if np.isnan(Cw):
                Cw = 0.0
            if np.isnan(Crp):
                Crp = 0.0
            if np.isnan(Rstall):
                Rstall = 0.0
            if np.isnan(Cb):
                Cb = 0.0
            Cwmass.append(Cw)
            Crpmass.append(Crp)
            Rstallmass.append(Rstall)
            Cbmass.append(Cb)
            Cws.append(Cw)
            Crps.append(Crp)
            Rstalls.append(Rstall)
            Cbs.append(Cb)
            tcools.append(TCool(1e6*yrins,star,cloud))
        # Set list of quartiles
        Cwmass.sort()
        Crpmass.sort()
        Rstallmass.sort()
        Cbmass.sort()
        for q in quartiles:
            QCw[q].append(Cwmass[q])
            QCrp[q].append(Crpmass[q])
            QRstall[q].append(Rstallmass[q])
            QCb[q].append(Cbmass[q])
        # Make contour lines for each percentile
        for p in range(0,50):
            Crpcontour[ic,p] = Crpmass[p]
            Crpcontour[2*ncloudmasses-ic-1,p] = Crpmass[99-p]
            Cwcontour[ic,p] = Cwmass[p]
            Cwcontour[2*ncloudmasses-ic-1,p] = Cwmass[99-p]
            Rstallcontour[ic,p] = Rstallmass[p]
            Rstallcontour[2*ncloudmasses-ic-1,p] = Rstallmass[99-p]
            Cbcontour[ic,p] = Cbmass[p]
            Cbcontour[2*ncloudmasses-ic-1,p] = Cbmass[99-p]
    #print "Cws", Cws
    #print "Crps", Crps
    '''
    Make Cw, Crp plot
    '''
    xmasses = np.array(xmasses)
    plt.clf()
    plotscatter = False
    ared = 0.8 # slight alpha to make difference in lightness between lines
    if plotscatter:
        # Plot as simple scatter plot
        plt.scatter(xmasses, Cws,3,"r",edgecolors='none',label="$C_{\mathrm{w}}$ (SFE="+str(int(sfe*100))+"\%)")
        plt.scatter(xmasses, Crps,3,"b",edgecolors='none',label="$C_{\mathrm{rp}}$ (SFE="+str(int(sfe*100))+"\%)")
    else:
        # Plot as surface plots with quartiles
        for p in range(0,50):
            plt.gca().fill(taba,Cwcontour[:,p],alpha=0.02,
                            edgecolor='none',facecolor="b")
            plt.gca().fill(taba,Crpcontour[:,p],alpha=0.02*ared,
                            edgecolor='none',facecolor="r")
        dolabels = True
        for line, q in zip(lines, quartiles):
            if dolabels:
                plt.plot(cloudmasses,QCw[q],linestyle=line,color="b",label="$C_{\mathrm{w}}$ (SFE="+str(int(sfe*100))+"\%)")
                plt.plot(cloudmasses,QCrp[q],linestyle=line,color="r",label="$C_{\mathrm{rp}}$ (SFE="+str(int(sfe*100))+"\%)",
                    alpha=ared)
            else:
                plt.plot(cloudmasses,QCw[q],linestyle=line,color="b")
                plt.plot(cloudmasses,QCrp[q],linestyle=line,color="r",alpha=ared)
            dolabels = False
    # Set axis limits
    ylims = [1e-3,1e3]
    # Plot line at Cw,Crp=1
    plt.plot([xmasses.min(),xmasses.max()],[1,1],"k-",alpha=0.25)
    # Vesc line
    # (commented out because rstall is a better measure)
    # if Mesc < cloudmasses.max():
    #     plt.plot(Mesc*np.ones(2),[1e-3,1e3],"k:")
    #     mesctxt = plt.text(Mesc*1.05,1e2,r"$C_{esc}>1$")
    # Make plot limits, labels, etc
    plt.ylim(ylims)
    plt.xlim([cloudmasses[0],cloudmasses[-1]])
    plt.xscale("log")
    plt.yscale("log")
    plt.ylabel("$C_{\mathrm{w}},C_{\mathrm{rp}}$")
    plt.xlabel("Cloud Mass / "+Msolar)
    surfstr = "$10^{"+str(int(np.log10(surfdens)))+"}$"
    plt.title("$\Sigma_0$ = "+surfstr+" M$_{\odot}$/pc$^2$, Z="+metaltxt.replace("p","."))
    if legend:
        plt.legend(frameon=False)
    cooltxt = ""
    if cool:
        cooltxt = "_cool"
    surftxt = ""
    if not plotscatter:
        surftxt = "_surfplot"
    plt.savefig("plots/clustermasses"+cooltxt+surftxt+"_Sigma"+str(int(surfdens))+"_Z"+metaltxt+".pdf")
    '''
    Make Rstall plot
    '''
    xmasses = np.array(xmasses)
    plt.clf()
    plotscatter = False
    pink = r"#FF69B4"
    if plotscatter:
        # Plot as simple scatter plot
        plt.scatter(xmasses, Rstalls,3,"k",edgecolors='none',label="$C_{\mathrm{w}}$ (SFE="+str(int(sfe*100))+"\%)")
    else:
        # Plot as surface plots with quartiles
        for p in range(0,50):
            plt.gca().fill(taba,Rstallcontour[:,p],alpha=0.02,
                            edgecolor='none',facecolor="k")
            plt.gca().fill(taba,Cbcontour[:,p],alpha=0.02,
                            edgecolor='none',facecolor=pink)
        dolabels = True
        for line, q in zip(lines, quartiles):
            if dolabels:
                plt.plot(cloudmasses,QRstall[q],linestyle=line,color="k",label="$r_{launch}$  (SFE="+str(int(sfe*100))+"\%)")
                plt.plot(cloudmasses,QCb[q],linestyle=line,color=pink,label="$C_{\mathrm{B}}$ (SFE="+str(int(sfe*100))+"\%)")
            else:
                plt.plot(cloudmasses,QRstall[q],linestyle=line,color="k")
                plt.plot(cloudmasses,QCb[q],linestyle=line,color=pink)
            dolabels = False
    # Set axis limits
    ylims = [1e-3,1e3]
    # Vesc line
    #if Mesc < cloudmasses.max():
    #    plt.plot(Mesc*np.ones(2),[1e-3,1e3],"k:")
    #    mesctxt = plt.text(Mesc*1.05,1e2,r"$C_{esc}>1$")
    # Plot line at rstall, Cb=1
    plt.plot([xmasses.min(),xmasses.max()],[1,1],"k-",alpha=0.25)
    # Make plot limits, labels, etc
    plt.ylim(ylims)
    plt.xlim([cloudmasses[0],cloudmasses[-1]])
    plt.xscale("log")
    plt.yscale("log")
    plt.ylabel("$r_{launch}$ / pc, $C_{B}$")
    plt.xlabel("Cloud Mass / "+Msolar)
    surfstr = "$10^{"+str(int(np.log10(surfdens)))+"}$"
    plt.title("$\Sigma_0$ = "+surfstr+" M$_{\odot}$/pc$^2$, Z="+metaltxt.replace("p","."))
    if legend:
        plt.legend(frameon=False)
    cooltxt = ""
    if cool:
        cooltxt = "_cool"
    surftxt = ""
    if not plotscatter:
        surftxt = "_surfplot"
    plt.savefig("plots/clusterrstalls_cbs"+cooltxt+surftxt+"_Sigma"+str(int(surfdens))+".pdf")
    # Cooling time
    # plt.clf()
    # plt.plot(masses, tcools,"k-")
    # plt.yscale("log")
    # plt.savefig("plots/cooltest.pdf")

def plotsurfaces_radii(metal=0.014, cool=False):
    # Stellar mass
    starmasses = stars.starmasses
    nstarmass = len(starmasses)
    # Cloud radii
    ncloudradii = 20
    cloudradii = np.logspace(-2,2,ncloudradii)*pc
    #surfs = np.logspace(1e1,1e4,nsurf)*Msun/pc**2
    # Other setup stuff
    #ci = 1e6
    #if metal == 0.002:
    #    ci = 2e6
    metaltxt = str(metal).replace(".","p")
    t = 1e6*yrins
    nbinary = 2
    Cws  = np.zeros((ncloudradii, nstarmass,nbinary))
    Crps = np.zeros((ncloudradii, nstarmass,nbinary))
    Rstalls = np.zeros((ncloudradii, nstarmass,nbinary))
    # Escape velocity stuff
    #testcloud = Cloud(100*Msun,surf,ci)
    #Cesc = GravCriterion(t,None,testcloud) # Doesn't use any of the star properties
    #Mesc = 100.0*Msun / Cesc**4
    # Iterate over cloud masses
    for icloud in range(0,ncloudradii):
        cloudradius = cloudradii[icloud]
        n0 = 1e4 # Doesn't really affect anything in Cw, Crp
        cloud = CloudFromr0n0(cloudradius,n0,metal,accreting=False,gravityOn=True)
        print cloud.M0/Msun, cloud.r0/cloudradius
        # Iterate over stellar masses
        for istar in range(0,nstarmass):
            starmass = starmasses[istar]
            star = stars.Star(starmass,metal)
            binary = stars.Binary(starmass,metal)
            # Single star
            Cw = WindCriterion(t,star,cloud,cool=cool)
            Crp = RadiationPressureCriterion(t,star,cloud)
            Rstall = rstall(t,star,cloud,cool=cool) / pc # MAKE IN PC, NOT IN CM
            Cws[icloud,istar,0] = Cw
            Crps[icloud,istar,0] = Crp
            Rstalls[icloud,istar,0] = Rstall
            # Binaries
            Cw = WindCriterion(t,binary,cloud,cool=cool)
            Crp = RadiationPressureCriterion(t,binary,cloud)
            Rstall = rstall(t,binary,cloud,cool=cool) / pc # MAKE IN PC, NOT IN CM
            Cws[icloud,istar,1] = Cw
            Crps[icloud,istar,1] = Crp
            Rstalls[icloud,istar,1] = Rstall
    for data, name in zip([Rstalls,Cws,Crps],["Rstall","Cw","Crp"]):
        fontsize = 12
        cooltxt = ""
        if cool:
            cooltxt = "_cool"
        plotname = "plots/surface_"+name+"_byradius"+cooltxt+"_Z"+metaltxt+".pdf"
        print "Making figure", plotname
        plt.clf()
        X, Y = np.meshgrid(starmasses, np.log10(cloudradii/pc))
        Z = data[:,:,0]
        Zbin = data[:,:,1]
        print Zbin.min(), Z.min()
        print Zbin.max(), Z.max()
        V = [0.1,1.0,10.0]
        contours = plt.contour(X,Y,Z,V,colors="black", zorder=2)
        contourbin = plt.contour(X,Y,Zbin,V,colors="black",linestyles="dashed", zorder=10)
        #surfstr = "$10^{"+str(int(np.log10(surfdens)))+"}$"
        #plt.title("$\Sigma_0$ = "+surfstr+" M$_{\odot}$/pc$^2$")
        plt.legend(frameon=False)
        plt.clabel(contours, inline=True, fontsize=fontsize, fmt=clabelfmt)
        plt.clabel(contourbin, inline=True, fontsize=fontsize, fmt=clabelfmt)
        plt.title("Z="+metaltxt.replace("p","."))
        Xlow = np.log10(cloudradii.min()/pc)
        Xhigh = np.log10(cloudradii.max()/pc)
        Ylow = starmasses.min()
        Yhigh = starmasses.max()
        extents = np.array([Ylow,Yhigh,Xlow,Xhigh])
        plt.imshow(np.log(Z), extent=extents, origin='lower',cmap=cmaps[name], 
                    alpha=0.5, aspect="auto",vmin=-2,vmax=2,zorder=2)
        plt.xlabel("Stellar mass / M$_{\odot}$")
        plt.ylabel("HII Region radius / pc")
        labels = {"Cw":"$C_{w}$","Crp":"$C_{rp}$","Rstall":"$r_{launch}$ / pc"}
        plt.colorbar(label="log("+labels[name]+")")
        plt.savefig(plotname)

def testvalues(mcloud=1e2, surfdens=1e2, metal=0.014, cool=False):
    # Surface density
    surf = surfdens*Msun/pc**2
    # Stellar mass
    starmasses = stars.starmasses
    nstarmass = len(starmasses)
    # Cloud mass
    cloudmass = mcloud*Msun
    #surfs = np.logspace(1e1,1e4,nsurf)*Msun/pc**2
    # Other setup stuff
    #ci = 1e6
    t = 1e6*yrins
    cloud = Cloud(cloudmass,surf,metal)
    star = stars.Star(30,metal)
    # Single star
    Cw = WindCriterion(t,star,cloud,cool=cool)
    Crp = RadiationPressureCriterion(t,star,cloud)
    Cesc = GravCriterion(t,star,cloud)
    print "Diagnostic test for models"
    print "Cloud mass:", mcloud
    print "Surface density:", surfdens
    print "Stellar mass:", star.mass
    print "---"
    print "QH", star.NPhotons(t)
    print "pw", star.WindMomentumRate(t)
    print "alphaB", cloud.alphaB(star,t)
    print "---"
    print "Cw", Cw
    print "Crp", Crp
    print "Cesc", Cesc

def testcluster(starmasses, mcloud=1e2, surfdens=1e2, metal=0.014, cool=True):
    # Surface density
    surf = surfdens*Msun/pc**2
    # Cloud mass
    cloudmass = mcloud*Msun
    #surfs = np.logspace(1e1,1e4,nsurf)*Msun/pc**2
    # Other setup stuff
    #ci = 1e6
    t = 1e6*yrins
    cloud = Cloud(cloudmass,surf,metal)
    star = stars.Cluster(starmasses,metal)
    # Single star
    Cw = WindCriterion(t,star,cloud,cool=cool)
    Crp = RadiationPressureCriterion(t,star,cloud)
    Cesc = GravCriterion(t,star,cloud)
    rivsrw = RivsRw(t,star,cloud)
    print "Diagnostic test for models"
    print "Cloud mass:", mcloud
    print "Surface density:", surfdens
    print "Stellar mass:", star.mass
    print "---"
    print "QH", star.NPhotons(t)
    print "pw", star.WindMomentumRate(t)
    print "alphaB", cloud.alphaB(star,t)
    print "---"
    print "Cw", Cw
    print "Crp", Crp
    print "Cesc", Cesc
    print "ri vs rw", rivsrw

def testamun(cloudlimit,metal=0.014):
    stars.verbose = False
    stardict = OrderedDict()
    stardict["M4_DIM"] = "30.76135824  9.91160177 10.78105538 30.88566842"
    stardict["M4_BRIGHT"] = "68.1151573  12.31290415 18.6270794"
    stardict["M5_BRIGHT"] = "68.1151573   12.31290415  18.6270794   17.99462301  54.03748279                                                         21.89887211  14.29378963   8.45167459  13.16419119  13.3033125                                                           8.1647781   13.24716712  11.69671807   8.78087569  17.76555596                                                         10.38842351  13.82479644  10.50574254  17.1176801   33.49659646                                                         10.11080321   8.06573658  40.48259705  14.47780878  13.91573894                                                         26.77688144  10.16726215   8.58661552  12.14698929   9.33191645                                                         22.45632112   8.70231881   9.65762067  14.84393557  69.3124608                                                           8.15676994  15.75376102  16.3256691  116.6760136   10.01863118                                                         14.93343147  48.65989081  18.08942282  16.49473292  52.50504075                                                          8.66518639  52.29503419   8.41391199  25.3740757   14.41284974                                                          9.24361472  10.45473048  25.13638936   8.21158929   8.5090458                                                          12.8556772   17.45313138  11.93737507  11.89786598  26.56479031                                                         17.96341229  20.03575382   8.68421145  11.93361329   9.70675954                                                          9.78582461  59.77432517  11.91998905   8.02185748   8.90316542                                                         21.35150655  20.0910768   10.38562128  19.1763621   30.18802907                                                          8.0423042    8.22455485  23.95680713   9.33382088  13.54407749                                                         13.08748301  12.06308378  16.11165484   9.41965059   8.72817561                                                          9.40438291  14.89985653   9.59336323   9.3200083   12.17783394                                                         14.86751831  13.092701    13.16279308   8.39313896  22.77406357                                                         31.52336573  47.59377135  26.6444641   13.1995285   14.86842749                                                         29.16517201  34.75631563  13.7007171   12.61265948  12.58744205                                                         26.38923677  14.32181683  22.67841826  19.76825932  11.34557163                                                         21.35464255  17.06847223   9.33015567  16.46023292  10.68777153                                                         17.81216676  53.86577454   9.44911824  11.6313694   12.6339549                                                          14.29699997   9.6830226   17.04460023  29.40376663  14.94980941                                                         20.89174197  12.00227164   9.76538001   9.06528517  10.14158163                                                         32.03191603  24.98976282  20.55123821  10.67849032  25.03152366                                                         15.25001643   8.43324663  10.4625876   14.65436319  14.06622894                                                         62.7554428    8.49982726  10.18532546   9.66235306  21.59003031                                                         71.19751123  10.34781608  18.92430094   8.07514988   8.14642606                                                          9.61999551"
    #ci = 1e6
    corem4 = CloudFromr0n0(1.0*pc,432.0,metal)
    corem5 = CloudFromr0n0(1.0*pc,390.0,metal)
    isom5 = CloudFromr0n0(10.0*pc,32.1,metal)
    surftomass_m4 = {1e2:2832.69, 1e3:348.95,"iso1pc":corem4,"isom5":isom5}
    surftomass_m5 = {1e2:104107.7, 1e3:23291.1,"iso1pc":corem5,"isom5":isom5}
    for name, starstr in stardict.iteritems():
        print "---"
        print name
        if "M4" in name:
            cloud = surftomass_m4[cloudlimit]
        if "M5" in name:
            cloud = surftomass_m5[cloudlimit]
        # Is cloud limit a surface density or a description?
        if cloudlimit == "iso1pc" or cloudlimit == "isom5":
            print "Core radius, density, mass, Sigma",
            print cloud.r0/pc, cloud.n0, cloud.M0 / Msun, cloud.Sigma0 / (Msun/pc**2)
            mcloud = cloud.M0 / Msun
            surfdensval = cloud.Sigma0  / (Msun/pc**2)
        else:
            mcloud = cloud
            surfdensval = cloudlimit
        starlist = starstr.split(" ")
        starlist = filter(None, starlist)
        starlist = [float(s) for s in starlist]
        if cloudlimit == "iso1pc" or cloudlimit == "isom5":
            starlist = [starlist[0]]
        testcluster(starlist,mcloud,surfdensval)
        print "---"

def testamun2(metal=0.014):
    t = 1e6*yrins
    cool = True
    for mstar in [30,60,120]:
        star = stars.Star(mstar,metal)
        ci = SoundSpeedFromT(1e4)#(cloud,star,t)
        cloud = CloudFromr0n0(1.0*pc,432.0,metal,ci=ci)
        rlaunch = rstall_iononly(t,star,cloud) / pc
        Cw = WindCriterion(t,star,cloud,cool=cool)
        Crp = RadiationPressureCriterion(t,star,cloud)
        Cesc = GravCriterion(t,star,cloud)
        CB = BreakoutCriterion(t,star,cloud)
        rivsrw = RivsRw(t,star,cloud)
        print "rlaunch, Cw, Crp, Cesc, CB, rwvsri, ci", rlaunch, Cw, Crp, Cesc, CB, 1/rivsrw, ci


if __name__=="__main__":
    #testvalues(cool=True)
    #testamun(cloudlimit="iso1pc")
    #testamun(cloudlimit=1e3)
    #testamun(cloudlimit="isom5")
    testamun2()
    '''
    dolegend = True
    ionisedtemperatureversusmass(soundspeed=True)
    ionisedtemperatureversusmass(soundspeed=False)
    plotsurfaces_radii(metal=0.002, cool=True)
    plotsurfaces_radii(metal=0.014, cool=True)
    for surf in [1e2,1e3,1e4]:
        plotsurfaces(surf,metal=0.014,cool=True)
        plotsurfaces(surf,metal=0.002,cool=True)
        plotvsclustermass(surf,metal=0.014,sfe=0.1,cool=True,legend=dolegend)
        plotvsclustermass(surf,metal=0.002,sfe=0.1,cool=True,legend=dolegend)
        dolegend = False
    '''
