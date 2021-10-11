"""
Analytic solutions for shell dynamics
Sam Geen, September 2020
"""

import os, sys

import numpy as np

import stars
import units

import matplotlib.pyplot as plt

import equilibriumtemperature
import ionisedtemperatures

import solvehiiprofile

import rdmfile

from matplotlib.ticker import MaxNLocator

gamma = equilibriumtemperature.gamma

Msun = "M$_{\odot}$"

NCLOUDS = np.array([1e2,1e3,1e4])

def alpha_B_HII(temperature):
    """
    Calculate the HII recombination rate
    This is the rate at which ionised hydrogen recombines 
      into neutral hydrogen
    Total recombinations per second per unit volume = alpha * ne * nH
    ne = electron number density
    nH = hydrogen number density

    Parameters
    ----------

    temperature: float
        Temperature in K

    Returns
    -------

    alpha_B_HII : float
        The recombination rate
    """
    # HII recombination rate
    # input  : T in K
    # output : HII recombination rate (in cm3 / s)
    l = 315614./temperature
    a = 2.753e-14 * l**1.5 / (1. + (l/2.74)**0.407)**2.242
    return a   

# Use a reference time to set stars' ages?
REFERENCETIME = 1e5*units.year

"""
-------------------------------------------------------------------------
GAS STUFF
-------------------------------------------------------------------------
"""
def SoundSpeed(star,t):
    #HACK - ignore change in Ti due to ionised gas density
    Ti = ionisedtemperatures.FindTemperature(star.Teff(t),star.metal)
    # Find ci
    # 2.0 factor is because of free electrons
    ci = np.sqrt(2.0 * gamma * units.kB * Ti * units.X / units.mH)
    return ci

def Tion(star,t):
    #HACK - ignore change in Ti due to ionised gas density
    Ti = ionisedtemperatures.FindTemperature(star.Teff(t),star.metal)
    return Ti

def AlphaB(star,t):
    #HACK - ignore change in Ti due to ionised gas density
    Ti = ionisedtemperatures.FindTemperature(star.Teff(t),star.metal)
    return alpha_B_HII(Ti)

def FiducialDustSigma(metal,solarsigma=1e-21):
    # Fiducial sigma to use - assume a Draine 2011 sigma scaled by metallicity
    return solarsigma * metal/0.014

class Cloud(object):
    def __init__(self,n0,r0,w,sigmaDust,Omega=4.0*np.pi):
        self._n0 = n0
        self._r0 = r0
        self._w = w
        self._sigmaDust = sigmaDust
        self._Omega = Omega

    def __str__(self):
        return "CLOUD_"+str(self._n0)+"_"+str(self._r0)+"_"+str(self._w)+"_"+str(self._sigmaDust)+"_"+str(self._Omega)

    @property
    def n0(self):
        return self._n0

    @property
    def r0(self):
        return self._r0

    @property
    def w(self):
        return self._w

    @property
    def sigmaDust(self):
        return self._sigmaDust

    @property
    def Omega(self):
        return self._Omega

    def nAtR(self, r):
        return self._n0 * (r / self._r0)**(-self._w)

    def mAtR(self, r):
        n0 = self._n0
        r0 = self._r0
        w = self._w
        threeminusw = 3.0 - w
        return 4.0 * np.pi / threeminusw * n0 * units.mH / units.X * r0 ** w * r**threeminusw

"""
-------------------------------------------------------------------------
DRAINE PARAMETERS
-------------------------------------------------------------------------
"""

def DraineOpticalDepth(star,timeins,nrms,sigmaDust):
    '''
    Optical depth reference value tau_{d,0} in Draine (2011)
    '''
    metal = star.metal
    QH = star.LIonising(timeins) / star.EPhoton(timeins)
    Tion = ionisedtemperatures.FindTemperature(star.Teff(t),star.metal)
    taud0 = 2.10 * (QH/1e49*nrms/1e3)**(1.0/3.0) * (Tion / 1e4)**0.28 * (sigmaDust/1e-21)
    return taud0

def DraineBeta(star,timeins):
    '''
    beta = Lnonionising / Lionising in Draine (2011)
    '''
    beta = star.LNonIonising(timeins) / star.LIonising(timeins)
    return beta

def DraineGamma(star,timeins,sigmaDust):
    '''
    gamma = dimensionless parameter for effectiveness of radiation pressure in Draine (2011)
    '''
    metal = star.metal
    Eion = star.EPhoton(timeins)
    Tion = ionisedtemperatures.FindTemperature(star.Teff(t),star.metal)
    gamma = 11.2 * (Tion / 1e4)**1.83 * (18 * units.eV / Eion) * (sigmaDust/1e-21)
    return gamma

"""
-------------------------------------------------------------------------
ANALYTIC SOLUTIONS
-------------------------------------------------------------------------
"""

def WindPressure(star,cloud,r,t,cooled=False,tref=None):
    if cooled:
        raise NotImplementedError
    Omega = cloud.Omega
    vol = (Omega / 3.0) * r**3
    if tref is not None:
        Eemitted = star.WindLuminosity(tref)*tref
    else:
        Eemitted = star.WindEnergy(t)
    w = cloud.w
    efact = (5.0-w)/(11.0-w)
    E = Eemitted * efact # From Weaver
    #print ("E R VOL", E, r, vol)
    return E / vol / 1.5 # Factor of 3/2 for 3 degrees of freedom, i.e. E = 3/2 N kB T = 3/2 P V

def RadiationPressure(star,cloud,r,t):
    L = star.LIonising(t) + star.LNonIonising(t)
    # NOTE: Omega isn't used here because light travels in straight lines
    return L / (4.0 * np.pi * r**2 * units.c)

def dRdTforw2Old(star,cloud,radiation=True):
    # Expansion rate of wind bubble for w=2
    Omega = cloud.Omega
    t= 1e5*units.year
    Lw = star.WindLuminosity(t)
    Lrad = 0.0
    if radiation:
        Lrad = star.LIonising(t) + star.LNonIonising(t)
    rho0 = cloud.n0 * units.mH /  units.X
    r0 = cloud.r0
    w = cloud.w
    # Solve Av^3 + Bv^2 + Cv + D = 0
    A = Omega * rho0 * r0*r0
    B = 0.0
    C = -Omega / (4.0 * np.pi * units.c)*Lrad
    D = -2.0*(5.0-w)/(11.0-w) * Lw
    roots = np.roots([A,B,C,D])
    root = roots[np.where(roots > 0.0)][0]
    #print("ROOTS", roots, root)
    return np.real(root)

def dRdT(star,cloud,t):
    # Expansion rate of wind bubble
    w = cloud.w
    return (3.0 / (5.0-w)) * WindRadius(star,cloud,t) / t

def Ri(star,cloud,r,t):
    shellrs,shellns,nmean,nrms,shellmass = solvehiiprofile.findprofile(star,cloud,r,t)
    ri = shellrs[-1]
    return ri

def RwfromRiforw2(star,cloud,riIn):
    rw = riIn+0.0
    ri = riIn+0.0
    v2 = dRdTforw2(star,cloud)
    error = 1.0
    while np.abs(error) > 0.0001:
        t = rw / v2
        ri = Ri(star,cloud,rw,t)
        error = (ri/riIn)-1
        rw *= 1.0-error
    return rw

    # r0 = 0.01 * units.pc
    # r1 = 100000*units.pc
    # # Do a binary partition search for the overflow radius
    # while (r1 / r0) > 1.001:
    #     rm = 10.0**(0.5*(np.log10(r0) + np.log10(r1)))
    #     ro0, to0 = findroto(r0)
    #     rom, tom = findroto(rm)
    #     ro1, to1 = findroto(r1)
    #     if rom > 0:
    #         r1 = rm
    #     else:
    #         r0 = rm


def dRidT(star,cloud,r,t):
    r2 = r*1.001
    drdt = dRdT(star,cloud,t)
    t2 = r2 / drdt
    shellrs,shellns,nmean,nrms,shellmass = solvehiiprofile.findprofile(star,cloud,r,t)
    ri = shellrs[-1]
    shellrs,shellns,nmean,nrms,shellmass = solvehiiprofile.findprofile(star,cloud,r2,t2)
    ri2 = shellrs[-1]
    dridt = (ri2 - ri) / (t2 - t)
    return dridt

def dRdTforw2(star,cloud):
    # Expansion rate of wind bubble, special case for w=2 where it is independent of time
    n0 = cloud.n0
    r0 = cloud.r0
    w = 2
    tref = 1e5*units.year
    Lw = star.WindLuminosity(tref)
    rho0 = n0 * units.mH / units.X
    rw = (3.0 / (5.0-w)) * (WindRadiusConstant(star,cloud)*Lw/(rho0*r0**w))**(1.0/3.0)
    return rw

def WindRadiusConstant(star,cloud):
    n0 = cloud.n0
    r0 = cloud.r0
    w = cloud.w
    Omega = cloud.Omega
    Aw = (1-w/3.0)*(1-w/5.0)**3.0 / (1-2*w/7.0) / (1 - w/11.0)
    return 4.0*np.pi / Omega * Aw * (250.0/308.0 / np.pi)

def WindRadius(star,cloud,t):
    n0 = cloud.n0
    r0 = cloud.r0
    w = cloud.w
    Lw = star.WindLuminosity(t)
    rho0 = n0 * units.mH / units.X
    rw = (WindRadiusConstant(star,cloud)*Lw*t**3.0/(rho0*r0**w))**(1.0/(5.0-w))
    return rw

def WindTime(star,cloud,r):
    n0 = cloud.n0
    r0 = cloud.r0
    w = cloud.w
    t = 1e5 * units.year # fiducial wind luminosity reference time
    Lw = star.WindLuminosity(t)
    rho0 = n0 * units.mH / units.X
    tw = (r**(5.0-w)*(rho0*r0**w)/(WindRadiusConstant(star,cloud)*Lw))**(1.0/3.0)
    return tw

def GravityPressure(star,cloud,r,t):
    # Find edge of ionised shell
    shellrs,shellns,nmean,nrms,shellmass = solvehiiprofile.findprofile(star,cloud,r,t)
    ri = shellrs[-1]
    n0 = cloud.n0
    r0 = cloud.r0
    w = cloud.w
    M = cloud.mAtR(ri)
    Omega = cloud.Omega
    Pgrav = units.G*M*M / Omega / ri**4
    return Pgrav

def OverflowParameter(star,cloud,r,t):
    n0 = cloud.n0
    r0 = cloud.r0
    w = cloud.w
    LHS = n0 * r0**w * r ** (1.0-w) / (3.0 - w)
    niiVSnrms = 1.4
    niavVSnrms = 1.0/1.4
    QH = star.NIonising(t)
    ci = SoundSpeed(star,t)
    alphaB = AlphaB(star,t)
    Pw = WindPressure(star,cloud,r,t)
    Prad = RadiationPressure(star,cloud,r,t)
    RHS = QH * units.mH * ci ** 2 * niiVSnrms * niavVSnrms / (4.0 * np.pi * r**2 * alphaB * units.X * (Pw + Prad))
    #print("OVERFLOWPARAMS", LHS, RHS, QH, ci, r, alphaB, Pw, Prad)
    return RHS / LHS # if > 1, outflow happens
    
def IonisedShellThicknessOld(star,cloud,r,t,thickshell=False,uniform=False):
    if thickshell:
        ri, rw = IonisedRadiusThickShell_w2(star,cloud,r,t)
        deltar = ri - rw
    else:
        n0 = cloud.n0
        r0 = cloud.r0
        w = cloud.w
        tref = 1e5*units.year
        if uniform:
            niiVSnrms = 1.0
            niavVSnrms = 1.0
            Prad = 0.0
        else:
            # Estimate from Draine+ 2011
            niiVSnrms = 1.4
            niavVSnrms = 1.0/1.4
            Prad = RadiationPressure(star,cloud,r,tref)
        Omega = cloud.Omega # solid angle
        QH = star.NIonising(t)
        ci = SoundSpeed(star,t)
        alphaB = AlphaB(star,t)
        Pw = WindPressure(star,cloud,r,t,tref=tref)
        deltar = QH / (alphaB * Omega) * niiVSnrms**2.0 * (units.mH / units.X)**2.0 * ci**4.0 * (r * (Pw + Prad))**(-2.0)
    return deltar

def IonisedRadiusThickShell_w2Old(star,cloud,r,t):
    '''
    Ionised radius assuming a uniform photoionised shell
    '''
    Lw = star.WindLuminosity(t)
    n0 = cloud.n0
    r0 = cloud.r0
    w = cloud.w
    Omega = cloud.Omega # solid angle
    # Ionising stuff
    QH = star.NIonising(t)
    ci = SoundSpeed(star,t)
    alphaB = AlphaB(star,t)
    # Wind stuff
    Pw = WindPressure(star,cloud,r,t)
    Lw = star.WindLuminosity(t)
    v2 = dRdTforw2(star,cloud)
    # Calculate radii
    rw = v2*t
    A = 3.63 * (Omega**2 / 4.0 / np.pi) * (QH/alphaB) * (ci**2 * units.mH * v2 / Lw / units.X)**2
    ri = rw * (1 + A * rw)**(1.0/3.0)
    return ri, rw

def ShellThicknessConstant(star,cloud):
    n0 = cloud.n0
    r0 = cloud.r0
    w = cloud.w
    r0w = r0**w
    Omega = cloud.Omega # solid angle
    t = 1e5 * units.year # fiducial time
    QH = star.NIonising(t)
    ci = SoundSpeed(star,t)
    alphaB = AlphaB(star,t)
    Lw = star.WindLuminosity(t)
    powerindex = 3.0/(5.0-4.0*w)
    mH = units.mH
    X = units.X
    rho0 = n0 * mH/X
    Bw = Omega**2 / 16.0 / np.pi
    Bw *= ((11-w)/(5-w))**2
    Bw *= (mH/X)**2
    Bw *= QH/alphaB
    Bw *= ci**4
    Bw *= (WindRadiusConstant(star,cloud)/Lw**2/rho0/r0w)**(2.0/3.0)
    return Bw

def IonisedShellThickness(star,cloud,rw,t):
    '''
    Ionised radius assuming a uniform photoionised shell
    '''
    w = cloud.w
    powerindex = (2.0*w-1.0)/3.0
    Bw = ShellThicknessConstant(star,cloud)
    ri = rw * (1.0 + 3.0 * Bw * rw**powerindex)**(1.0/3.0)
    dr = ri - rw
    return dr

def IonisedRadiusShell_w2(star,cloud,r,t):
    '''
    Ionised radius assuming a uniform photoionised shell, force w=2
    '''
    cloud2 = Cloud(cloud.n0,cloud.r0,2.0,cloud.sigmaDust,cloud.Omega)
    return IonisedShellThickness(star,cloud2,r,t)

    
def IonisedRadiusThickShell_w2(star,cloud,rw,t):
    '''
    Legacy function name, kept to keep code working, we assume a thick shell always for now
    '''
    cloud2 = Cloud(cloud.n0,cloud.r0,2.0,cloud.sigmaDust,cloud.Omega)
    w = cloud2.w
    powerindex = (2.0*w-1.0)/3.0
    Bw = ShellThicknessConstant(star,cloud2)
    ri = rw * (1.0 + 3.0 * Bw * rw**powerindex)**(1.0/3.0)
    dr = ri - rw
    return ri, rw

def NeutralShellSoundSpeed(star,cloud,r,t):
    Pw = WindPressure(star,cloud,r,t)
    Prad = RadiationPressure(star,cloud,r,t)
    Ptot = Pw + Prad
    Tshell = equilibriumtemperature.NeutralTemperatureFromPressure(Ptot,star.metal)
    cshell = np.sqrt(gamma * units.kB * Tshell * units.X / units.mH)
    return cshell

def NeutralShellThickness(star,cloud,r,t):
    n0 = cloud.n0
    r0 = cloud.r0
    w = cloud.w
    Omega = cloud.Omega # solid angle
    drdt = ExpansionRate(star,cloud,r,t)
    deltar = (NeutralShellSoundSpeed(star,cloud,r,t) / drdt)**2.0 * r / (3.0 - w)
    return deltar

def CalculateOutflowRadiusTimeOld(star,cloud,rcutoff = None,tcutoff = None):
    '''
    # FIRST VERSION: ITERATING OVER PRESSURE EQUATION
    r = 1e-2 * units.pc
    dt = 1e2 * units.year
    t = 0.0
    Cout = 0.0
    niter = 0
    while Cout < 1.0:
        t += dt
        #print(t)
        Cout = OverflowParameter(star,cloud,r,t)
        #print(Cout, r / units.pc, t / units.year / 1e3)
        Pw = WindPressure(star,cloud,r,t)
        Prad = RadiationPressure(star,cloud,r,t)
        Ptot = Pw + Prad
        deltar = IonisedShellThickness(star,cloud,r,t)
        drdt = np.sqrt(Ptot * units.X / (cloud.nAtR(r) * units.mH))
        if rcutoff is not None:
            if r > rcutoff:
                return cutoff, t
        if tcutoff is not None:
            if t > tcutoff:
                return cutoff, t
        r += drdt * dt
        niter += 1
        '''
    # SECOND VERSION: FINDING POINT ON R = v2 . t LINE WHERE OUTFLOW OCCURS
    start = 1e0*units.year
    end = 1e9*units.year
    steps = 10.0
    drdt = dRdTforw2(star,cloud)
    rout = None
    tout = None
    for niter in range(0,10):
        ts = np.logspace(np.log10(start),np.log10(end),10)
        rs = drdt*ts
        for r, t in zip(rs,ts):
            Cout = OverflowParameter(star,cloud,r,t)
            #print(niter,r,t,Cout)
            if Cout >= 1.0 or t > star.Lifetime():
                #print("BREAK")
                rout = r
                tout = t
                start = tprev
                end = t
                break
            tprev = t
    if tout > star.Lifetime():
        return None, None, None
    deltar = IonisedShellThickness(star,cloud,rout,tout)
    return rout, tout, deltar

def CalculateOutflowRadiusTimew2Old(star,cloud):
    n0 = cloud.n0
    r0 = cloud.r0
    w = cloud.w
    Omega = cloud.Omega # solid angle
    t = 1e5 * units.year # fiducial time
    QH = star.NIonising(t)
    ci = SoundSpeed(star,t)
    alphaB = AlphaB(star,t)
    niiVSnrms = 1.4
    niavVSnrms = 1.0/1.4
    Lw = star.WindLuminosity(t)
    Lrad = star.LIonising(t) + star.LNonIonising(t)
    v2 = dRdTforw2(star,cloud)
    # Calculate the overflow radius and time if w=2 (exact solution)
    ro2 = alphaB * units.X * n0 * r0*r0 / (QH * units.mH * ci*ci * niavVSnrms * niiVSnrms)
    ro2 *= (2.0*(5.0-w)/(11.0-w) * Lw / v2 + Lrad * Omega / (4.0 * np.pi * units.c))
    to2 = ro2 / v2
    deltar = IonisedShellThickness(star,cloud,ro2,to2)
    return ro2, to2, deltar

def CalculateOutflowRadiusTimew2(star,cloud):
    cloud2 = Cloud(cloud.n0,cloud.r0,2.0,cloud.sigmaDust,cloud.Omega)
    return CalculateOutflowRadiusTime(star,cloud2)

def CalculateOutflowRadiusTime(star,cloud):
    n0 = cloud.n0
    r0 = cloud.r0
    w = cloud.w
    r0w = r0**w
    Omega = cloud.Omega # solid angle
    t = 1e5 * units.year # fiducial time
    QH = star.NIonising(t)
    ci = SoundSpeed(star,t)
    alphaB = AlphaB(star,t)
    Lw = star.WindLuminosity(t)
    powerindex = 3.0/(5.0-4.0*w)
    mH = units.mH
    X = units.X
    Cw = WindRadiusConstant(star,cloud)
    roverflow = Omega/4.0/np.pi * QH/alphaB * (3.0-w)/2.0 * (11.0-w)/(5.0-w) * (mH/X)**(2.0/3.0) * ci**2.0
    roverflow *= Lw**(-2.0/3.0) * (n0*r0w)**(-4.0/3.0) * Cw**(1.0/3.0)
    roverflow = roverflow**powerindex
    toverflow = WindTime(star,cloud,roverflow)
    droverflow = IonisedShellThickness(star,cloud,roverflow,toverflow)
    return roverflow, toverflow, droverflow

def CalculateOutflowRadiusTime2(star,cloud):
    '''
    This version includes a factor that accounts for mass enclosed out to ri
    '''

    n0 = cloud.n0
    r0 = cloud.r0
    w = cloud.w
    if w != 2:
        print("Doesn't work if w!=2, sorry")
        raise ValueError
    r0w = r0**w
    Omega = cloud.Omega # solid angle
    t = 1e5 * units.year # fiducial time
    QH = star.NIonising(t)
    ci = SoundSpeed(star,t)
    alphaB = AlphaB(star,t)
    Lw = star.WindLuminosity(t)
    powerindex = 3.0/(5.0-4.0*w)
    mH = units.mH
    X = units.X
    rho0 = n0 * mH/X
    Aw = WindRadiusConstant(star,cloud)
    Bw = ShellThicknessConstant(star,cloud)
    Dw = (1.5 * QH/alphaB * mH/X * ci**2 * Omega/(4.0*np.pi)  / (n0*r0w))**3.0 * (Aw / Lw **2 / (rho0 * r0w))
    A = Dw
    B = 0.0
    C = -3.0*Bw
    D = -1.0
    # Solve Av^3 + Bv^2 + Cv + D = 0
    roots = np.roots([A,B,C,D])
    root = roots[np.where(roots > 0.0)][0]
    root = np.real(root)
    roverflow = root
    #roverflow = Omega/4.0/np.pi * QH/alphaB * (3.0-w)/2.0 * (11.0-w)/(5.0-w) * (mH/X)**(2.0/3.0) * ci**2.0
    #roverflow *= Lw**(-2.0/3.0) * (n0*r0w)**(-4.0/3.0) * Cw**(1.0/3.0)
    #roverflow = roverflow**powerindex
    toverflow = WindTime(star,cloud,roverflow)
    droverflow = IonisedShellThickness(star,cloud,roverflow,toverflow)
    return roverflow, toverflow, droverflow

def ExpansionRate(star,cloud,r,t):
    Pw = WindPressure(star,cloud,r,t)
    Prad = RadiationPressure(star,cloud,r,t)
    dr = IonisedShellThickness(star,cloud,r,t)
    Ptot = Pw + Prad
    drdt = np.sqrt(Ptot * units.X / (cloud.nAtR(r) * units.mH))
    return drdt

def MinResolutionNeeded(star,cloud,r,t):
    '''
    Calculation of the minimum grid resolution needed to trap photons
    '''
    n0 = cloud.n0
    r0 = cloud.r0
    w = cloud.w
    Omega = cloud.Omega # solid angle
    QH = star.NIonising(t)
    alphaB = AlphaB(star,t)
    res = Omega * alphaB * (n0 * r0**w) ** 2 * r**(4 - 2 * w) / (3.0 - w)**2 / QH
    return res

def CoolingRate(star,cloud,r,t):
    '''
    Cooling rate dE/dt of the wind bubble
    '''
    # Calculate the core temperature and density
    Lw = star.WindLuminosity(t)
    Zsolar = star.metal/0.014
    Omega = cloud.Omega # solid angle
    Cconduct = 6e-7 # erg / s / cm / K^{-7/2} from Mac Low & McCray 1988
    Ccool = 1e-22 # erg cm^3 / s
    Tcut = 1e5 # K
    w = cloud.w
    efact = (5.0-w) / (11.0-w)
    Tc = ((2.0 * efact * Lw) / (Omega * Cconduct * r))**(2.0/7.0)
    Eb = 2.0 * efact * Lw * t
    pc = Eb / (Omega * r**3)
    nc = pc / (units.kB * Tc)
    xcut = 1 - (Tcut / Tc)**2.5
    # Calculate the integral for the interior of the bubble
    # (Requires some integration solving, I used Wolfram Alpha)
    structureIntegral = -23 * xcut*xcut - 50 * xcut + 625
    structureIntegral /= (1-xcut)**(2.0/25.0)
    structureIntegral -= 625.0
    structureIntegral *= 25.0/1104.0
    dEdt = structureIntegral * Omega * r**3 * nc*nc * Ccool * Zsolar * Tc**(-0.7)
    return dEdt, Eb


"""
-------------------------------------------------------------------------
PLOTTING ROUTINES
-------------------------------------------------------------------------
"""

def CloudColour(ncloud,metal):
    #nclouds = np.log10(np.array([3e1,1e2,3e2,1e3,3e3]))
    nclouds = np.log10(NCLOUDS)
    col = np.log10(ncloud) - nclouds[0]
    col /= (nclouds[-1] - nclouds[0])
    cmap = plt.get_cmap("copper")
    #if metal == 0.014:
    #    cmap = plt.get_cmap("copper")
    #elif metal == 0.002:
    #    cmap = plt.get_cmap("winter")
    #else:
    #    print("Metallicity",metal,"not implemented here")
    #    raise NotImplementedError
    return cmap(col)

def CloudLine(metal):
    if metal == 0.014:
        return (8, 0)
    elif metal == 0.002:
        return (6, 2)
    else:
        print("Metallicity",metal,"not implemented here")
        raise NotImplementedError

def PlotdRdTforw2(metals=[0.014],rotating=True):
    plt.clf()
    #nclouds = [3e1,3e2,3e3]
    nclouds = NCLOUDS
    masses = np.arange(20,121,5)
    leglines = None
    rdm = rdmfile.RDMFile(__file__)
    lines = {0.014:"-",0.002:"--"}
    for metal in metals:
        for ncloud in nclouds:
            cloud = Cloud(ncloud,1.0*units.pc,2,FiducialDustSigma(metal))
            drdts = []
            drdtsnp = []
            print ("NCLOUD",ncloud)
            print(".......................")
            for starmass in masses:
                print("STARMASS",starmass)
                star = stars.Star(starmass, metal, rotating=rotating)
                star.ForceAge(REFERENCETIME)
                drdt = dRdTforw2(star,cloud)
                drdts.append(drdt)
                #drdt = dRdTforw2(star,cloud)
                #drdtsnp.append(drdt)
            col = CloudColour(ncloud,metal)
            dashes = CloudLine(metal)
            label = None
            if metal == metals[0]:
                label = str(ncloud)+" cm$^{-3}$"
            l1, = plt.plot(masses,np.array(drdts)/1e5,color=col,label=label,linestyle=lines[metal])
            rdm.AddPoints(masses,np.array(drdts)/1e5,label=label)
            #l2, = plt.plot(masses,np.array(drdtsnp)/1e5,color=col,linestyle=":")
            #rdm.AddPoints(masses,np.array(drdtsnp)/1e5,label=label+"_np")
            #if leglines is None:
            #    leglines = [l1,l2]
    leg1 = plt.legend(ncol=2,fontsize="small",frameon=False,
                    title="$n_0$, where $n(r) = n_0 (r / 1~\mathrm{pc})^{-2}$",title_fontsize="small")     
    leg1._legend_box.align = "right"
    plt.text(120,15,"$Z="+str(metals[0])+"$",fontsize="small",horizontalalignment="right",verticalalignment="top")
    #leg2 = plt.legend(leglines,["Winds + Radiation Pressure","Winds Only"],
    #        fontsize="small",frameon=False,loc="upper left")
    plt.gca().add_artist(leg1)
    plt.gca().yaxis.set_ticks_position('both')
    plt.xlabel("Initial Stellar Mass / "+Msun)
    plt.ylabel("$\mathrm{d}r_{w,2}/\mathrm{d}r \equiv v_2$ / km/s")
    plt.yscale("log")
    plt.ylim([1,400])
    filename = "../plots/drdtforw2_metals"+str(metals[0])+".pdf"
    plt.savefig(filename)
    rdm.Write(filename)

def PlotdRdTforPowerIndex(metal=0.014,rotating=True):
    plt.clf()
    #nclouds = [3e1,3e2,3e3]
    ncloud = 4000.0
    masses = np.arange(20,121,5)
    starmass = 35.0
    star = stars.Star(starmass, metal, rotating=rotating)
    star.ForceAge(REFERENCETIME)
    leglines = None
    rdm = rdmfile.RDMFile(__file__)
    ws = [0.0,1.0,2.0]
    lines = ["--",":","-"]
    linestyles = {w:line for w, line in zip(ws, lines)}
    kyr = 1e3*units.year
    times = np.linspace(0.0,1.0,100)*200*kyr
    for w in ws:
        cloud = Cloud(ncloud,1.0*units.pc,w,FiducialDustSigma(metal))
        drdts = []
        drdtsnp = []
        drdts = dRdT(star,cloud,times)
        col = "k" # CloudColour(ncloud,metal)
        dashes = CloudLine(metal)
        label = "$\omega="+str(int(w))+"$"
        l1, = plt.plot(times/kyr,np.array(drdts)/1e5,color=col,label=label,linestyle=linestyles[w])
        rdm.AddPoints(times/kyr,np.array(drdts)/1e5,label=label)
        #import pdb; pdb.set_trace()

    leg1 = plt.legend(fontsize="small",frameon=False)#,title="$(\Omega = 4 \pi)$",title_fontsize="small")     
    leg1._legend_box.align = "right"
    #plt.text(120,15,"$Z="+str(metals[0])+"$",fontsize="small",horizontalalignment="right",verticalalignment="top")
    #leg2 = plt.legend(leglines,["Winds + Radiation Pressure","Winds Only"],
    #        fontsize="small",frameon=False,loc="upper left")
    #plt.gca().add_artist(leg1)
    plt.gca().yaxis.set_ticks_position('both')
    plt.xlabel("Time / kyr")
    plt.ylabel("$\mathrm{d}r_w/\mathrm{d}r$ / km/s")
    plt.xscale("log")
    plt.yscale("log")
    plt.ylim([7,30])
    plt.xlim([5,200])
    plt.yticks(ticks=[10,20,30], labels=["10","20","30","40"])

    #ax = plt.figure().gca()
    #ax.yaxis.set_major_locator(MaxNLocator(integer=True))   
    filename = "../plots/drdtforpowerindex"+str(metal)+".pdf"
    plt.savefig(filename)
    rdm.Write(filename)

def PlotdRdTforSolidAngle(metal=0.014,rotating=True):
    plt.clf()
    #nclouds = [3e1,3e2,3e3]
    ncloud = 4000.0
    masses = np.arange(20,121,5)
    starmass = 35.0
    star = stars.Star(starmass, metal, rotating=rotating)
    star.ForceAge(REFERENCETIME)
    leglines = None
    rdm = rdmfile.RDMFile(__file__)
    Omegas = np.array([1.0,2.0,4.0])
    lines = ["--",":","-"]
    linestyles = {Omega:line for Omega, line in zip(Omegas, lines)}
    kyr = 1e3*units.year
    times = np.linspace(0.0,1.0,100)*200*kyr
    for Omega in Omegas:
        cloud = Cloud(ncloud,1.0*units.pc,2,FiducialDustSigma(metal),Omega=Omega*np.pi)
        drdts = []
        drdtsnp = []
        drdts = dRdT(star,cloud,times)
        col = "k" # CloudColour(ncloud,metal)
        dashes = CloudLine(metal)
        label = "$\Omega="+str(int(Omega))+"\pi$"
        l1, = plt.plot(times/kyr,np.array(drdts)/1e5,color=col,label=label,linestyle=linestyles[Omega])
        rdm.AddPoints(times/kyr,np.array(drdts)/1e5,label=label)
        #import pdb; pdb.set_trace()

    leg1 = plt.legend(fontsize="small",frameon=False)     
    leg1._legend_box.align = "right"
    #plt.text(120,15,"$Z="+str(metals[0])+"$",fontsize="small",horizontalalignment="right",verticalalignment="top")
    #leg2 = plt.legend(leglines,["Winds + Radiation Pressure","Winds Only"],
    #        fontsize="small",frameon=False,loc="upper left")
    #plt.gca().add_artist(leg1)
    plt.gca().yaxis.set_ticks_position('both')
    plt.xlabel("Time / kyr")
    plt.ylabel("$\mathrm{d}r_w/\mathrm{d}r$ / km/s")
    plt.xscale("log")
    #plt.yscale("log")
    plt.ylim([5,20])
    plt.xlim([5,200])
    #plt.yticks(ticks=[10,20,30,40], labels=["10","20","30","40"])

    #ax = plt.figure().gca()
    plt.gca().yaxis.set_major_locator(MaxNLocator(integer=True))   
    filename = "../plots/drdtforsolidangle"+str(metal)+".pdf"
    plt.savefig(filename)
    rdm.Write(filename)

def PlotThicknessRatioNew(metal=0.014,rotating=True):
    plt.clf()
    rdm = rdmfile.RDMFile(__file__)
    #nclouds = [3e1,3e2,3e3]
    mstars = [30,60,120]
    nclouds = NCLOUDS
    masses = np.arange(20,121,5)
    leglines = []
    Omegas = np.array([1.0,2.0,4.0])
    lines = ["-","--","-."]
    linestyles = {mstar:line for mstar, line in zip(mstars, lines)}
    kyr = 1e3*units.year
    for starmass in mstars:
        star = stars.Star(starmass, metal, rotating=rotating)
        star.ForceAge(REFERENCETIME)
        for ncloud in nclouds:    
            times = np.linspace(0.0,1.0,100)*1000*kyr
            cloud = Cloud(ncloud,1.0*units.pc,2,FiducialDustSigma(metal),Omega=4.0*np.pi)
            rw = WindRadius(star,cloud,times)
            dr = IonisedShellThickness(star,cloud,rw,times)
            ro, to, dro = CalculateOutflowRadiusTime(star,cloud)
            ri = rw + dr
            ratios = dr / rw
            overflowed = False
            if rw.max() > ro:
                overflowed = True
            ratios = ratios[rw <= ro]
            times = times[rw < ro]
            col = CloudColour(ncloud,metal)
            dashes = CloudLine(metal)
            label = ""
            if starmass == mstars[0]:
                label = str(int(ncloud))+" cm$^{-3}$"
            l1, = plt.plot(times/kyr,np.array(ratios),color=col,label=label,linestyle=linestyles[starmass])
            rdm.AddPoints(times/kyr,np.array(ratios),label=label)
            if overflowed:
                plt.scatter([times[-1]/kyr],ratios[-1],color=col)
            if ncloud == nclouds[0]:
                leglines += [l1]
            #import pdb; pdb.set_trace()
    leg1 = plt.legend(ncol=2,fontsize="small",frameon=False, loc = "lower right",
        title="$n_0$, where $n(r) = n_0 (r / 1~\mathrm{pc})^{-2}$",title_fontsize="small") 

    leg2 = plt.legend(leglines,[str(starmass)+" M$_{\odot}$" for starmass in mstars],ncol=3,
        fontsize="small",frameon=False,loc="upper right")
    plt.gca().add_artist(leg1)
    plt.text(600,3e-3,"$Z="+str(metal)+"$",fontsize="small",horizontalalignment="center",verticalalignment="bottom") 

    #plt.text(120,15,"$Z="+str(metals[0])+"$",fontsize="small",horizontalalignment="right",verticalalignment="top")
    #leg2 = plt.legend(leglines,["Winds + Radiation Pressure","Winds Only"],
    #        fontsize="small",frameon=False,loc="upper left")
    #plt.gca().add_artist(leg1)
    plt.plot([0,1000],[0.2599,0.2599],"k-",alpha=0.3)
    plt.gca().yaxis.set_ticks_position('both')
    plt.xlabel("Time / kyr")
    plt.ylabel("$(r_i - r_w) / r_w$")
    #plt.xscale("log")
    plt.yscale("log")
    plt.ylim([0.0001,10])
    plt.xlim([5,1000])
    #plt.yticks(ticks=[10,20,30,40], labels=["10","20","30","40"])

    #ax = plt.figure().gca()
    #plt.gca().yaxis.set_major_locator(MaxNLocator(integer=True))   
    filename = "../plots/thicknessratio_new_"+str(metal)+".pdf"
    plt.savefig(filename)
    rdm.Write(filename)

def PlotMinResolutionNeeded(metals=[0.014],rotating=True):
    plt.clf()
    #nclouds = [3e1,3e2,3e3]
    nclouds = NCLOUDS
    masses = np.arange(20,121,5)
    metal = 0.014
    # NOTE: As long as w = 2, radius is unimportant
    r = 0.1*units.pc
    rdm = rdmfile.RDMFile(__file__)
    for ncloud in nclouds[::-1]:
        ls = []
        for metal in metals:
            cloud = Cloud(ncloud,1.0*units.pc,2,FiducialDustSigma(metal))
            minreses = []
            print ("NCLOUD",ncloud)
            print(".......................")
            for starmass in masses:
                print("STARMASS",starmass)
                star = stars.Star(starmass, metal, rotating=rotating)
                star.ForceAge(REFERENCETIME)
                ro2, to2, deltar = CalculateOutflowRadiusTimew2(star,cloud)
                minres = MinResolutionNeeded(star,cloud,ro2,to2)
                minreses.append(minres)
            col = CloudColour(ncloud,metal)
            dashes = CloudLine(metal)
            label = None
            if metal == metals[0]:
                label = str(ncloud)+" cm$^{-3}$"
            l, = plt.plot(masses,np.array(minreses) / units.pc,color=col,label=label,dashes=dashes)
            rdm.AddPoints(masses,np.array(minreses) / units.pc,label=label)
            ls.append(l)
    leg1 = plt.legend(ncol=2,fontsize="x-small",frameon=False,
                    title="$n_0$, where $n(r) = n_0 (r / 1~\mathrm{pc})^{-2}$",title_fontsize="small")     
    leg1._legend_box.align = "right"
    leg2 = plt.legend(ls,["$Z="+str(metal)+"$" for metal in metals],
            fontsize="small",frameon=False,loc="upper right",ncol=2)
    plt.gca().add_artist(leg1)
    plt.xlabel("Initial Stellar Mass / "+Msun)
    plt.ylabel("Minimum Resolution Needed / pc")
    plt.yscale("log")
    plt.ylim([1e-6,1e3])
    filename = "../plots/minimumresolution.pdf"
    plt.savefig(filename)
    rdm.Write(filename)

def PlotThicknessRatiow2(metal=0.014,rotating=True,vsR = False,justthickness=False,thickshell=False):
    plt.clf()
    #nclouds = [3e1,3e2,3e3]
    nclouds = NCLOUDS
    #masses = np.arange(20,121,5)
    masses = [30,60,120]
    lines = ["-","--","-."]
    iline = 0
    leglines = []
    rdm = rdmfile.RDMFile(__file__)
    for ncloud in nclouds:
        cloud = Cloud(ncloud,1.0*units.pc,2,FiducialDustSigma(metal))
        drdts = []
        drdtsnp = []
        print ("NCLOUD",ncloud)
        print(".......................")
        for starmass in masses:
            print("STARMASS",starmass)
            star = stars.Star(starmass, metal, rotating=rotating)
            star.ForceAge(REFERENCETIME)
            rs = np.logspace(-2,2,100)*units.pc
            drdt = dRdTforw2(star,cloud)
            ts = rs / drdt
            ro2, to2, deltaro2 = CalculateOutflowRadiusTimew2(star,cloud)
            overflowed = False
            if rs[-1] > ro2:
                overflowed = True
            rs = rs[rs < ro2]
            ts = ts[ts < to2]
            ys = []
            for r, t in zip(rs, ts):
                rthick = IonisedShellThickness(star,cloud,r,t,True)
                rthin = IonisedShellThickness(star,cloud,r,t,False)
                # VALUES INDICATE A MAX +/- 20% CORRECTION TO THE THIN SHELL APPROXIMATION
                #print ("THICKNESS CORRECTION:", rthick / rthin)
                deltar = IonisedShellThickness(star,cloud,r,t)
                if not justthickness:
                    ys.append(deltar/r)
                else:
                    ys.append(deltar / units.pc)
            col = CloudColour(ncloud, metal)
            label = None
            if starmass == masses[0]:
                label = str(ncloud)+" cm$^{-3}$"
            line = lines[iline]
            xs = ts/units.year/1e6
            if vsR:
                xs = rs / units.pc
            l, = plt.plot(xs,np.array(ys),color=col,linestyle=line,label=label)
            rdm.AddPoints(xs,np.array(ys),label=label)
            if overflowed:
                plt.scatter(xs[-1],ys[-1],color=col)
            if ncloud == nclouds[0]:
                leglines += [l]
            iline = (iline+1) % len(lines)
    # Line at 0.01, 0.1, 1
    if not justthickness:
        plt.plot([0,100],[0.01,0.01],color="k",linestyle="-",alpha=0.2)
        plt.plot([0,100],[0.1,0.1],color="k",linestyle="-",alpha=0.2)
        plt.plot([0,100],[1,1],color="k",linestyle="-",alpha=0.2)
    else:
        plt.plot([1e-2,1e2],[1e-2,1e2],color="k",alpha=0.2,linestyle=":")
        plt.plot([1e-2,1e2],[1e-3,1e1],color="k",alpha=0.2,linestyle=":")
        plt.plot([1e-2,1e2],[1e-4,1e0],color="k",alpha=0.2,linestyle=":")
    if justthickness and metal==0.002:
        loc = "lower right"
        leg1 = plt.legend(ncol=1,fontsize="small",frameon=False,
                        title="$n_0$, where $n(r) = n_0 (r / 1~\mathrm{pc})^{-2}$",title_fontsize="small",loc=loc)     
        leg1._legend_box.align = "right"
    if justthickness and metal==0.014:
        loc = "lower right"
        leg2 = plt.legend(leglines,[str(starmass)+" M$_{\odot}$" for starmass in masses],
                fontsize="small",frameon=False,loc=loc)
        #plt.gca().add_artist(leg1)
    plt.xlabel("Time / Myr")
    if vsR:
        plt.xlabel("$r_w$ / pc")
    if not justthickness:
        plt.ylabel("$\Delta r_i / r_w$")
    else:
        plt.ylabel("$\Delta r_i$ / pc")
    plt.xscale("log")
    plt.yscale("log")
    xlims = plt.gca().get_xlim()
    ylims = plt.gca().get_ylim()
    plt.gca().yaxis.set_ticks_position('both')
    if vsR:
        plt.xlim([1e-2,30.0])
    #plt.ylim([1e-5,ylims[1]])
    ymax = 1e1
    plt.ylim([1e-5,ymax])
    xtxtloc = 20
    xtxtal = "right"
    if justthickness:
        xtxtloc = 0.011
        xtxtal = "left"
    plt.text(xtxtloc,ymax*0.7,"$Z="+str(metal)+"$",fontsize="small",horizontalalignment=xtxtal,verticalalignment="top")
    suffix = ""
    if vsR:
        suffix = "_vsR"
    suffix += "_Z"+str(metal)
    if thickshell:
        suffix += "_thickshell"
    if not justthickness:
        filename = "../plots/shellthicknessratio_w2"+suffix+".pdf"
    else:
        filename = "../plots/shellthickness_w2"+suffix+".pdf"
    plt.savefig(filename)
    rdm.Write(filename)

def PlotMeanShellDensity(metal=0.014,rotating=True,vsR=False,plotnrmsvsnmean=False):
    plt.clf()
    #nclouds = [3e1,3e2,3e3]
    nclouds = NCLOUDS
    #masses = np.arange(20,121,5)
    masses = [30,60,120]
    lines = ["-","--","-."]
    iline = 0
    leglines = []
    rdm = rdmfile.RDMFile(__file__)
    for ncloud in nclouds:
        cloud = Cloud(ncloud,1.0*units.pc,2,FiducialDustSigma(metal))
        drdts = []
        drdtsnp = []
        print ("NCLOUD",ncloud)
        print(".......................")
        for starmass in masses:
            print("STARMASS",starmass)
            star = stars.Star(starmass, metal, rotating=rotating)
            star.ForceAge(REFERENCETIME)
            rs = np.logspace(-2,2,100)*units.pc
            drdt = dRdTforw2(star,cloud)
            ts = rs / drdt
            ys = []
            for r, t in zip(rs, ts):
                shellrs,shellns,nmean,nrms,shellmass = solvehiiprofile.findprofile(star,cloud,r,t)
                if plotnrmsvsnmean:
                    ys.append(nrms/nmean)
                else:
                    ys.append(shellns[-1]/nrms)
            # Calculate overflow radius
            roana, toana, deltaroana = CalculateOutflowRadiusTimew2(star,cloud)
            ro2 = (shellmass * (3.0-cloud.w)/cloud.Omega / (cloud.n0*cloud.r0**cloud.w) * units.X/units.mH)**(1.0/(3.0-cloud.w))
            print ("RO / pc: numerical, analytic: ", ro2/units.pc, roana/units.pc)
            overflowed = False
            if rs[-1] > ro2:
                overflowed = True
            ts = ts[rs < ro2]
            ys = np.array(ys)
            ys = ys[rs < ro2]
            rs = rs[rs < ro2]
            # Make plot line
            col = CloudColour(ncloud, metal)
            label = None
            if starmass == masses[0]:
                label = str(ncloud)+" cm$^{-3}$"
            line = lines[iline]
            xs = ts/units.year/1e6
            if vsR:
                xs = rs / units.pc
            l, = plt.plot(xs,np.array(ys),color=col,linestyle=line,label=label)
            rdm.AddPoints(xs,np.array(ys),label=label)
            if len(xs) > 0:
                if overflowed:
                    plt.scatter(xs[-1],ys[-1],color=col)
            if ncloud == nclouds[0]:
                leglines += [l]
            iline = (iline+1) % len(lines)
    # Line at 0.01, 0.1, 1
    # if not justthickness:
    #     plt.plot([0,100],[0.01,0.01],color="k",linestyle="-",alpha=0.2)
    #     plt.plot([0,100],[0.1,0.1],color="k",linestyle="-",alpha=0.2)
    #     plt.plot([0,100],[1,1],color="k",linestyle="-",alpha=0.2)
    # else:
    #     plt.plot([1e-2,1e2],[1e-2,1e2],color="k",alpha=0.2,linestyle=":")
    #     plt.plot([1e-2,1e2],[1e-3,1e1],color="k",alpha=0.2,linestyle=":")
    #     plt.plot([1e-2,1e2],[1e-4,1e0],color="k",alpha=0.2,linestyle=":")
    if plotnrmsvsnmean and metal==0.002:
        loc = "upper right"
        leg1 = plt.legend(ncol=1,fontsize="small",frameon=False,
                        title="$n_0$, where $n(r) = n_0 (r / 1~\mathrm{pc})^{-2}$",title_fontsize="small",loc=loc)     
        leg1._legend_box.align = "right"
    if plotnrmsvsnmean and metal==0.014:
        loc = "upper right"
        leg2 = plt.legend(leglines,[str(starmass)+" M$_{\odot}$" for starmass in masses],
                fontsize="small",frameon=False,loc=loc)
        #plt.gca().add_artist(leg1)
    plt.xlabel("Time / Myr")
    if vsR:
        plt.xlabel("$r_w$ / pc")
    if plotnrmsvsnmean:
        plt.ylabel("$n_{i,rms} / \\bar{n_i}$")
    else:
        plt.ylabel("$n_{i}(r_i) / n_{i,rms}$")
    plt.xscale("log")
    plt.yscale("linear")
    xlims = plt.gca().get_xlim()
    ylims = plt.gca().get_ylim()
    plt.gca().yaxis.set_ticks_position('both')
    if vsR:
        plt.xlim([1e-2,30.0])
    #plt.ylim([1e-5,ylims[1]])
    #ymax = 1e1
    #plt.ylim([1e-5,ymax])
    xtxtloc = 20
    xtxtal = "right"
    #if plotnrmsvsnmean:
    #    xtxtloc = 0.011
    #    xtxtal = "left"
    #plt.text(xtxtloc,ymax*0.7,"$Z="+str(metal)+"$",fontsize="small",horizontalalignment=xtxtal,verticalalignment="top")
    suffix = ""
    if vsR:
        suffix = "_vsR"
    suffix += "_Z"+str(metal)
    if plotnrmsvsnmean:
        suffix += "_nrmsvsnmean"
    else:
        suffix += "_nedgevsnrms"
    filename = "../plots/shellmeandensity_w2"+suffix+".pdf"
    plt.savefig(filename)
    rdm.Write(filename)

def PlotShellThicknessNumerical(metal=0.014,rotating=True,plotdrdt=False):
    plt.clf()
    #nclouds = [3e1,3e2,3e3]
    nclouds = NCLOUDS
    #masses = np.arange(20,121,5)
    masses = [30,120]
    lines = ["-","--"]
    iline = 0
    leglines = []
    rdm = rdmfile.RDMFile(__file__)
    for ncloud in nclouds:
        cloud = Cloud(ncloud,1.0*units.pc,2,FiducialDustSigma(metal))
        drdts = []
        drdtsnp = []
        print ("NCLOUD",ncloud)
        print(".......................")
        for starmass in masses:
            print("STARMASS",starmass)
            star = stars.Star(starmass, metal, rotating=rotating)
            star.ForceAge(REFERENCETIME)
            rs = np.logspace(-2,3,1000)*units.pc
            drdt = dRdTforw2(star,cloud)
            ts = rs / drdt
            drnums = []
            ro2 = None
            for r, t in zip(rs, ts):
                shellrs,shellns,nmean,nrms,shellmass = solvehiiprofile.findprofile(star,cloud,r,t)
                drnum = (shellrs[-1]-shellrs[0])/r
                #rotmp = (shellmass * (3.0-cloud.w)/cloud.Omega / (cloud.n0*cloud.r0**cloud.w) * units.X/units.mH)**(1.0/(3.0-cloud.w))
                ri = shellrs[-1]
                sweptupmass = cloud.mAtR(ri)
                if sweptupmass <= shellmass and ro2 is None:
                    ro2 = r
                # Plot shell velocity instead?
                if plotdrdt:
                    r2 = r*1.001
                    t2 = r2 / drdt
                    shellrs,shellns,nmean,nrms,shellmass = solvehiiprofile.findprofile(star,cloud,r2,t2)
                    ri2 = shellrs[-1]
                    drdti = (ri2 - ri) / (t2 - t)
                    drnum = (drdti / drdt - 1) / drnum
                drnums.append(drnum)

            #dranas = IonisedShellThickness(star,cloud,rs,ts)/rs
            # Calculate overflow radius
            #roana, toana, deltaroana = CalculateOutflowRadiusTime2(star,cloud)
            #ro2 = (shellmass * (3.0-cloud.w)/cloud.Omega / (cloud.n0*cloud.r0**cloud.w) * units.X/units.mH)**(1.0/(3.0-cloud.w))
            if ro2 is None:
                ro2 = 1e6 * rs.max()
            #if roana is None:
            #    roana = 1e6 * rs.max()
            #print ("RO / pc: numerical, analytic: ", ro2/units.pc, roana/units.pc)
            overflowed = False
            if rs[-1] > ro2:
                overflowed = True
            #ts = ts[rs < ro2]
            #tanas = ts[rs < roana]
            drnums = np.array(drnums)
            drnums = drnums[rs < ro2]
            #dranas = dranas[rs < roana]
            #ranas = rs[rs < roana]
            rs = rs[rs < ro2]
            # Make plot line
            col = CloudColour(ncloud, metal)
            label = ""
            if starmass == masses[0]:
                label = "$10^{"+str(int(np.log10(ncloud)))+"}$ cm$^{-3}$"
            line = lines[iline]
            xs = rs / units.pc
            #xanas = ranas / units.pc
            l, = plt.plot(xs,np.array(drnums),color=col,linestyle=line,label=label)
            #l2, = plt.plot(xanas,np.array(dranas),color=col,linestyle=line,alpha=0.3)
            rdm.AddPoints(xs,np.array(drnums),label=label+" numerical")
            #rdm.AddPoints(xanas,np.array(dranas),label=label+" analytic")
            if len(xs) > 0:
                if overflowed:
                    plt.scatter(xs[-1],drnums[-1],color=col)
                    #plt.scatter(xanas[-1],dranas[-1],color=col,alpha=0.3)
            if ncloud == nclouds[0]:
                leglines += [l]
            iline = (iline+1) % len(lines)
    if metal==0.002:
        loc = "lower right"
        leg1 = plt.legend(ncol=1,fontsize="small",frameon=False,
                        title="$n_0$, where $n(r) = n_0 (r / 1~\mathrm{pc})^{-2}$",title_fontsize="small",loc=loc)     
        leg1._legend_box.align = "right"
    if metal==0.014:
        loc = "lower right"
        leg2 = plt.legend(leglines,[str(starmass)+" M$_{\odot}$" for starmass in masses],
                fontsize="small",frameon=False,loc=loc)
        #plt.gca().add_artist(leg1)
    plt.xlabel("Time / Myr")
    plt.xlabel("$r_w$ / pc")
    plt.xscale("log")
    plt.yscale("log")
    xlims = plt.gca().get_xlim()
    ylims = plt.gca().get_ylim()
    plt.gca().yaxis.set_ticks_position('both')
    plt.xlim([3e-2,300.0])
    if not plotdrdt:
        plt.ylim([1e-4,10])
        plt.ylabel("$(r_i - r_w)/r_w$")
    else:
        plt.ylabel("d$r_i$/d$t$ / d$r_w$/d$t - 1$")
        plt.ylabel("(d$r_i$/d$t$ / d$r_w$/d$t - 1) / ((r_i - r_w)/r_w)$")
    xlims = plt.gca().get_xlim()
    ylims = plt.gca().get_ylim()
    plt.text(xlims[0]*2,ylims[-1]*0.7,"$Z="+str(metal)+"$",fontsize="small",horizontalalignment="left",verticalalignment="top")
    xtxtloc = 20
    xtxtal = "right"
    suffix = ""
    suffix += "_Z"+str(metal)
    if not plotdrdt:
        filename = "../plots/num_shellthickness"+suffix+".pdf"
    else:
        filename = "../plots/num_relativevelocity"+suffix+".pdf"
    plt.savefig(filename)
    rdm.Write(filename)

def PlotShellThicknessAnalyticvsNumerical(metal=0.014,rotating=True):
    plt.clf()
    #nclouds = [3e1,3e2,3e3]
    nclouds = NCLOUDS
    #masses = np.arange(20,121,5)
    masses = [30,120]
    lines = ["-","--"]
    iline = 0
    leglines = []
    rdm = rdmfile.RDMFile(__file__)
    for ncloud in nclouds:
        cloud = Cloud(ncloud,1.0*units.pc,2,FiducialDustSigma(metal))
        drdts = []
        drdtsnp = []
        print ("NCLOUD",ncloud)
        print(".......................")
        for starmass in masses:
            print("STARMASS",starmass)
            star = stars.Star(starmass, metal, rotating=rotating)
            star.ForceAge(REFERENCETIME)
            rs = np.logspace(-2,3,1000)*units.pc
            drdt = dRdTforw2(star,cloud)
            ts = rs / drdt
            drnums = []
            ro2 = None
            for r, t in zip(rs, ts):
                shellrs,shellns,nmean,nrms,shellmass = solvehiiprofile.findprofile(star,cloud,r,t)
                drnums.append((shellrs[-1]-shellrs[0])/r)
                #rotmp = (shellmass * (3.0-cloud.w)/cloud.Omega / (cloud.n0*cloud.r0**cloud.w) * units.X/units.mH)**(1.0/(3.0-cloud.w))
                ri = shellrs[-1]
                sweptupmass = cloud.mAtR(ri)
                if sweptupmass <= shellmass and ro2 is None:
                    ro2 = r
            dranas = IonisedShellThickness(star,cloud,rs,ts)/rs
            # Calculate overflow radius
            roana, toana, deltaroana = CalculateOutflowRadiusTime2(star,cloud)
            #ro2 = (shellmass * (3.0-cloud.w)/cloud.Omega / (cloud.n0*cloud.r0**cloud.w) * units.X/units.mH)**(1.0/(3.0-cloud.w))
            if ro2 is None:
                ro2 = 1e6 * rs.max()
            if roana is None:
                roana = 1e6 * rs.max()
            print ("RO / pc: numerical, analytic: ", ro2/units.pc, roana/units.pc)
            overflowed = False
            if rs[-1] > ro2:
                overflowed = True
            #ts = ts[rs < ro2]
            tanas = ts[rs < roana]
            drnums = np.array(drnums)
            drnums = drnums[rs < ro2]
            dranas = dranas[rs < roana]
            ranas = rs[rs < roana]
            rs = rs[rs < ro2]
            # Make plot line
            col = CloudColour(ncloud, metal)
            label = ""
            if starmass == masses[0]:
                label = str(ncloud)+" cm$^{-3}$"
            line = lines[iline]
            xs = rs / units.pc
            xanas = ranas / units.pc
            l, = plt.plot(xs,np.array(drnums),color=col,linestyle=line,label=label)
            l2, = plt.plot(xanas,np.array(dranas),color=col,linestyle=line,alpha=0.3)
            rdm.AddPoints(xs,np.array(drnums),label=label+" numerical")
            rdm.AddPoints(xanas,np.array(dranas),label=label+" analytic")
            if len(xs) > 0:
                if overflowed:
                    plt.scatter(xs[-1],drnums[-1],color=col)
                    plt.scatter(xanas[-1],dranas[-1],color=col,alpha=0.3)
            if ncloud == nclouds[0]:
                leglines += [l]
            iline = (iline+1) % len(lines)
    if metal==0.002:
        loc = "lower right"
        leg1 = plt.legend(ncol=1,fontsize="small",frameon=False,
                        title="$n_0$, where $n(r) = n_0 (r / 1~\mathrm{pc})^{-2}$",title_fontsize="small",loc=loc)     
        leg1._legend_box.align = "right"
    if metal==0.014:
        loc = "lower right"
        leg2 = plt.legend(leglines,[str(starmass)+" M$_{\odot}$" for starmass in masses],
                fontsize="small",frameon=False,loc=loc)
        #plt.gca().add_artist(leg1)
    plt.xlabel("Time / Myr")
    plt.xlabel("$r_w$ / pc")
    plt.ylabel("$(r_i - r_w)/r_w$")
    plt.xscale("log")
    plt.yscale("log")
    xlims = plt.gca().get_xlim()
    ylims = plt.gca().get_ylim()
    plt.gca().yaxis.set_ticks_position('both')
    plt.xlim([3e-2,300.0])
    xtxtloc = 20
    xtxtal = "right"
    suffix = ""
    suffix += "_Z"+str(metal)
    filename = "../plots/shellthickness_anavsnum"+suffix+".pdf"
    plt.savefig(filename)
    rdm.Write(filename)

def PlotOverflowRadiusTimeNumerical(rotating=True,plotradius=True):
    plt.clf()
    #nclouds = [3e1,3e2,3e3]
    nclouds = NCLOUDS
    masses = np.arange(20,121,5)
    #masses = [30,60,120]
    metals = [0.002,0.014]
    lines = ["--","-"]
    metallines = {metal:line for metal, line in zip(metals,lines)}
    leglines = []
    rdm = rdmfile.RDMFile(__file__)
    for metal in metals:
        line = metallines[metal]
        # Plot star lifetime as a dashed line
        #if not plotradius:
        #    lifetimes = [stars.Star(starmass, metal, rotating=rotating).Lifetime() / units.year / 1e6 for starmass in masses]
        #    plt.plot(masses, lifetimes,linestyle=line,alpha=0.3,color="k")
        for ncloud in nclouds[::-1]:
            if metallines[metal] == "-":
                label = "$10^{"+str(int(np.log10(ncloud)))+"}$ cm$^{-3}$"
            else:
                label = ""
            ros = []
            tos = []
            cloud = Cloud(ncloud,1.0*units.pc,2,FiducialDustSigma(metal))
            drdts = []
            drdtsnp = []
            print ("NCLOUD",ncloud)
            print(".......................")
            for starmass in masses:
                print("STARMASS",starmass)
                star = stars.Star(starmass, metal, rotating=rotating)
                star.ForceAge(REFERENCETIME)
                rs = np.logspace(-2,3,100)*units.pc
                drdt = dRdTforw2(star,cloud)
                ts = rs / drdt
                #drnums = []
                def findroto(rin,forceTrue=False):
                    tin = rin / dRdTforw2(star,cloud)
                    shellrs,shellns,nmean,nrms,shellmass = solvehiiprofile.findprofile(star,cloud,rin,tin)
                    #drnums.append((shellrs[-1]-shellrs[0])/rin)
                    ri = shellrs[-1]
                    sweptupmass = cloud.mAtR(ri)
                    if sweptupmass <= shellmass or forceTrue:
                        return ri, tin
                    else:
                        return -200, -200
                r0 = 0.01 * units.pc
                r1 = 100000*units.pc
                # Do a binary partition search for the overflow radius
                while (r1 / r0) > 1.001:
                    rm = 10.0**(0.5*(np.log10(r0) + np.log10(r1)))
                    ro0, to0 = findroto(r0)
                    rom, tom = findroto(rm)
                    ro1, to1 = findroto(r1)
                    if rom > 0:
                        r1 = rm
                    else:
                        r0 = rm
                # Get a numerical result to good approximation
                ro2, to2 = findroto(r0,forceTrue=True)
                ros.append(ro2)
                tos.append(to2)

            ros = np.array(ros)
            tos = np.array(tos)
            xs = np.array(masses)
            xs = xs[ros > 0.0]
            tos = tos[ros > 0.0]
            ros = ros[ros > 0.0]
                
            # Make plot line
            col = CloudColour(ncloud, metal)
            #xs = masses
            ys = np.array(ros) / units.pc
            if not plotradius:
                ys = np.array(tos) / units.year / 1e6
            l, = plt.plot(xs,np.array(ys),color=col,linestyle=line,label=label)
            if ncloud == nclouds[0]:
                leglines.append(l)

            rdm.AddPoints(xs,np.array(ys),label=label+" numerical")
    # Do legends
    if plotradius:
        loc = "upper right"
        leg1 = plt.legend(ncol=1,fontsize="small",frameon=False,
                        #title="$n_0$, where $n(r) = n_0 (r / 1~\mathrm{pc})^{-2}$",
                        title_fontsize="small",loc=loc,
                        columnspacing = 0.8)     
        leg1._legend_box.align = "right"
    #leg1.set_label_position('bottom')
    else:
        leg2 = plt.legend(leglines,["$Z="+str(metal)+"$" for metal in metals],
                fontsize="small",frameon=False,loc="upper right",ncol=2)
    #plt.gca().add_artist(leg1)
    plt.xlabel("Stellar Mass / M$_{\odot}$")
    plt.ylabel("Overflow radius / pc")
    if not plotradius:
        plt.ylabel("Overflow time / Myr")
    #plt.xscale("log")
    plt.yscale("log")
    xlims = plt.gca().get_xlim()
    ylims = plt.gca().get_ylim()
    plt.gca().yaxis.set_ticks_position('both')
    #plt.xlim([3e-2,300.0])
    xtxtloc = 20
    xtxtal = "right"
    suffix = ""
    if plotradius:
        suffix = "radius"+suffix
    else:
        suffix = "time"+suffix
    filename = "../plots/num_overflow"+suffix+".pdf"
    plt.savefig(filename)
    rdm.Write(filename)

def PlotThicknessCorrection(metal=0.014,rotating=True,vsR = False,uniform=False):
    plt.clf()
    #nclouds = [3e1,3e2,3e3]
    nclouds = NCLOUDS
    #masses = np.arange(20,121,5)
    masses = [30,60,120]
    lines = ["-","--","-."]
    iline = 0
    leglines = []
    rdm = rdmfile.RDMFile(__file__)
    for ncloud in nclouds:
        cloud = Cloud(ncloud,1.0*units.pc,2,FiducialDustSigma(metal))
        drdts = []
        drdtsnp = []
        print ("NCLOUD",ncloud)
        print(".......................")
        for starmass in masses:
            print("STARMASS",starmass)
            star = stars.Star(starmass, metal, rotating=rotating)
            star.ForceAge(REFERENCETIME)
            rs = np.logspace(-2,2,100)*units.pc
            drdt = dRdTforw2(star,cloud)
            ts = rs / drdt
            ro2, to2, deltaro2 = CalculateOutflowRadiusTimew2(star,cloud)
            overflowed = False
            if rs[-1] > ro2:
                overflowed = True
            rs = rs[rs < ro2]
            ts = ts[ts < to2]
            ys = []
            for r, t in zip(rs, ts):
                rthick = IonisedShellThickness(star,cloud,r,t,True)
                rthin = IonisedShellThickness(star,cloud,r,t,False,uniform=uniform)
                ys.append(rthin / rthick)
            col = CloudColour(ncloud, metal)
            label = None
            if starmass == masses[0]:
                label = str(ncloud)+" cm$^{-3}$"
            line = lines[iline]
            xs = ts/units.year/1e6
            if vsR:
                xs = rs / units.pc
            l, = plt.plot(xs,np.array(ys),color=col,linestyle=line,label=label)
            rdm.AddPoints(xs,np.array(ys),label=label)
            if overflowed:
                plt.scatter(xs[-1],ys[-1],color=col)
            if ncloud == nclouds[0]:
                leglines += [l]
            iline = (iline+1) % len(lines)
    # Line at 0.01, 0.1, 1
    #plt.plot([0,100],[0.01,0.01],color="k",linestyle="-",alpha=0.2)
    #plt.plot([0,100],[0.1,0.1],color="k",linestyle="-",alpha=0.2)
    plt.plot([0,100],[1,1],color="k",linestyle="-",alpha=0.2)
    if metal==0.014:
        leg1 = plt.legend(ncol=1,fontsize="small",frameon=False,
                        title="$n_0$, where $n(r) = n_0 (r / 1~\mathrm{pc})^{-2}$",title_fontsize="small")     
        leg1._legend_box.align = "right"
        loc = "lower right"
        leg2 = plt.legend(leglines,[str(starmass)+" M$_{\odot}$" for starmass in masses],
                fontsize="small",frameon=False,loc=loc)
        plt.gca().add_artist(leg1)
    plt.xlabel("Time / Myr")
    if vsR:
        plt.xlabel("$r_w$ / pc")
    unitxt = "radiation-shaped $n_i(r)$"
    if uniform:
        unitxt = "uniform $n_i(r)$"
    plt.ylabel("$\Delta r_i / \Delta r_{i,thick}$ ("+unitxt+")")
    plt.xscale("log")
    plt.yscale("log")
    xlims = plt.gca().get_xlim()
    ylims = plt.gca().get_ylim()
    plt.gca().yaxis.set_ticks_position('both')
    if vsR:
        plt.xlim([1e-2,30.0])
    #plt.ylim([1e-5,ylims[1]])
    ymax = 4.0
    plt.ylim([0.25,ymax])
    plt.text(0.02,ymax*0.7,"$Z="+str(metal)+"$",fontsize="small",horizontalalignment="left",verticalalignment="top")
    suffix = ""
    if vsR:
        suffix = "_vsR"
    if uniform:
        suffix += "_uniformthinsoln"
    suffix += "_Z"+str(metal)
    filename = "../plots/shellthickness_correction_w2"+suffix+".pdf"
    plt.savefig(filename) 
    rdm.Write(filename)

def PlotShellExpansionCorrection(metal=0.014,rotating=True,vsR = False,ratio=False):
    plt.clf()
    #nclouds = [3e1,3e2,3e3]
    nclouds = NCLOUDS
    #masses = np.arange(20,121,5)
    masses = [30,60,120]
    lines = ["-","--","-."]
    iline = 0
    leglines = []
    rdm = rdmfile.RDMFile(__file__)
    for ncloud in nclouds:
        cloud = Cloud(ncloud,1.0*units.pc,2,FiducialDustSigma(metal))
        drdts = []
        drdtsnp = []
        print ("NCLOUD",ncloud)
        print(".......................")
        for starmass in masses:
            print("STARMASS",starmass)
            star = stars.Star(starmass, metal, rotating=rotating)
            star.ForceAge(REFERENCETIME)
            rs = np.logspace(-2,2,100)*units.pc
            drdt = dRdTforw2(star,cloud)
            ts = rs / drdt
            ro2, to2, deltaro2 = CalculateOutflowRadiusTimew2(star,cloud)
            overflowed = False
            if rs[-1] > ro2:
                overflowed = True
            rs = rs[rs <= ro2]
            ts = ts[ts <= to2]
            ys = []
            for r, t in zip(rs, ts):
                rthick, rw = IonisedRadiusThickShell_w2(star,cloud,r,t)
                ys.append(rthick)
            print(len(ys), len(ts))
            drdtthick = np.diff(ys)/np.diff(ts) / 1e5 # in km/s
            drdtthin = np.diff(rs)/np.diff(ts) / 1e5 # in km/s
             # To match positions of differenced velocities
            rnews = 0.5*(rs[1:] + rs[:-1])
            tnews = 0.5*(ts[1:] + ts[:-1])
            col = CloudColour(ncloud, metal)
            label = None
            if starmass == masses[0]:
                label = str(ncloud)+" cm$^{-3}$"
            line = lines[iline]
            xs = tnews/units.year/1e6
            if vsR:
                xs = rnews / units.pc
            if not ratio:  
                l, = plt.plot(xs,np.array(drdtthick),color=col,linestyle=line,label=label)
                plt.plot(xs,np.array(drdtthin),color=col,linestyle=line,label=label,alpha=0.5)
                rdm.AddPoints(xs,np.array(drdtthin),label=str(ncloud)+" cm$^{-3}$"+" THICK")
                rdm.AddPoints(xs,np.array(drdtthin),label=str(ncloud)+" cm$^{-3}$"+" THIN")
                if overflowed:
                    plt.scatter(xs[-1],(drdtthick)[-1],color=col)
                    plt.scatter(xs[-1],(drdtthin)[-1],color=col,alpha=0.5)
            else:
                l, = plt.plot(xs,np.array(drdtthick/drdtthin),color=col,linestyle=line,label=label)
                rdm.AddPoints(xs,np.array(drdtthick/drdtthin),label=label)
                if overflowed:
                    plt.scatter(xs[-1],(drdtthick/drdtthin)[-1],color=col)
            if ncloud == nclouds[0]:
                leglines += [l]
            iline = (iline+1) % len(lines)
    # Line at 0.01, 0.1, 1
    #plt.plot([0,100],[0.01,0.01],color="k",linestyle="-",alpha=0.2)
    #plt.plot([0,100],[0.1,0.1],color="k",linestyle="-",alpha=0.2)
    if ratio:
        plt.plot([0,100],[1,1],color="k",linestyle="-",alpha=0.2)
    if metal==0.014:
        leg1 = plt.legend(ncol=1,fontsize="small",frameon=False,
                        title="$n_0$, where $n(r) = n_0 (r / 1~\mathrm{pc})^{-2}$",title_fontsize="small")     
        leg1._legend_box.align = "right"
        loc = "lower right"
        leg2 = plt.legend(leglines,[str(starmass)+" M$_{\odot}$" for starmass in masses],
                fontsize="small",frameon=False,loc=loc)
        plt.gca().add_artist(leg1)
    plt.xlabel("Time / Myr")
    if vsR:
        plt.xlabel("$r_w$ / pc")
    if ratio:
        plt.ylabel("d$r$/d$t_{thick}$ / d$r$/d$t_{thin}$")    
    else:
        plt.ylabel("d$r$/d$t$")
    plt.xscale("log")
    #plt.yscale("log")
    xlims = plt.gca().get_xlim()
    ylims = plt.gca().get_ylim()
    plt.gca().yaxis.set_ticks_position('both')
    if vsR:
        plt.xlim([1e-2,30.0])
    #plt.ylim([1e-5,ylims[1]])
    ymax = ylims[1]
    if ratio:
        ymin = 0.4
    else:
        ymin = ylims[0]
    plt.ylim([ymin,ymax])
    plt.text(0.02,ymax*0.7,"$Z="+str(metal)+"$",fontsize="small",horizontalalignment="left",verticalalignment="top")
    suffix = ""
    if vsR:
        suffix += "_vsR"
    if ratio:
        suffix += "_ratio"
    suffix += "_Z"+str(metal)
    filename = "../plots/shellexpansion_correction_w2"+suffix+".pdf"
    plt.savefig(filename) 
    rdm.Write(filename)

def PlotOutflowRadii(metals=[0.014],rotating=True):
    plt.clf()
    #nclouds = [3e1,3e2,3e3]
    nclouds = NCLOUDS
    masses = np.arange(20,121,5)
    rdm = rdmfile.RDMFile(__file__)
    for ncloud in nclouds[::-1]:
        ls = []
        for metal in metals:
            cloud = Cloud(ncloud,1.0*units.pc,2,FiducialDustSigma(metal))
            radii = []
            print ("NCLOUD",ncloud)
            print(".......................")
            for starmass in masses:
                print("STARMASS",starmass)
                star = stars.Star(starmass, metal, rotating=rotating)
                star.ForceAge(REFERENCETIME)
                #r, t, dr = CalculateOutflowRadiusTime(star,cloud,rcutoff=1000*units.pc,tcutoff=1e7*units.year)
                r, t, dr = CalculateOutflowRadiusTimew2(star,cloud)
                radii.append(r)
            col = CloudColour(ncloud,metal)
            dashes = CloudLine(metal)
            label = None
            if metal == metals[0]:
                label=str(ncloud)+" cm$^{-3}$"
            l, = plt.plot(masses,np.array(radii) / units.pc,color=col,label=label,dashes=dashes)
            rdm.AddPoints(masses,np.array(radii) / units.pc,label=label)
            ls.append(l)
    leg1 = plt.legend(ncol=2,fontsize="small",frameon=False,
                    title="$n_0$, where $n(r) = n_0 (r / 1~\mathrm{pc})^{-2}$",title_fontsize="small")     
    leg1._legend_box.align = "left"
    #leg2 = plt.legend(ls,["$Z="+str(metal)+"$" for metal in metals],
    #        fontsize="small",frameon=False,loc="lower right",ncol=2)
    plt.gca().add_artist(leg1)
    plt.gca().yaxis.set_ticks_position('both')
    plt.xlabel("Initial Stellar Mass / "+Msun)
    plt.ylabel("Radius of Overflow $r_{o,2}$ / pc")
    plt.yscale("log")
    plt.ylim([1e-4,100])
    filename = "../plots/outflowradii.pdf"
    plt.savefig(filename)
    rdm.Write(filename)

def PlotOutflowTimes(metals=[0.014],rotating=True):
    plt.clf()
    #nclouds = [3e1,3e2,3e3]
    nclouds = NCLOUDS
    masses = np.arange(20,121,5)
    rdm = rdmfile.RDMFile(__file__)
    for ncloud in nclouds[::-1]:
        ls = []
        for metal in metals:
            cloud = Cloud(ncloud,1.0*units.pc,2,FiducialDustSigma(metal))
            times = []
            print ("NCLOUD",ncloud)
            print(".......................")
            for starmass in masses:
                print("STARMASS",starmass)
                star = stars.Star(starmass, metal, rotating=rotating)
                star.ForceAge(REFERENCETIME)
                #r, t, dr = CalculateOutflowRadiusTime(star,cloud,rcutoff=1000*units.pc,tcutoff=1e7*units.year)
                r, t, dr = CalculateOutflowRadiusTimew2(star,cloud)
                times.append(t)
            col = CloudColour(ncloud,metal)
            dashes = CloudLine(metal)
            label = None
            if metal == metals[0]:
                label=str(ncloud)+" cm$^{-3}$"
            l, = plt.plot(masses,np.array(times) / units.year / 1e6,color=col,label=label,dashes=dashes)
            rdm.AddPoints(masses,np.array(times) / units.year / 1e6,label=label)
            ls.append(l)
    #leg1 = plt.legend(ncol=2,fontsize="small",frameon=False,
    #                title="$n_0$, where $n(r) = n_0 (r / 1~\mathrm{pc})^2$",title_fontsize="small")     
    #leg1._legend_box.align = "right"
    leg2 = plt.legend(ls,["$Z="+str(metal)+"$" for metal in metals],
            fontsize="small",frameon=False,loc="upper right",ncol=2)
    #plt.gca().add_artist(leg1)
    plt.gca().yaxis.set_ticks_position('both')
    plt.xlabel("Initial Stellar Mass / "+Msun)
    plt.ylabel("Time of Overflow  $t_{o,2}$ / Myr")
    plt.yscale("log")
    plt.ylim([1e-4,30])
    filename = "../plots/outflowtimes.pdf"
    plt.savefig(filename) 
    rdm.Write(filename)

def PlotGravityVsBubblePressure(metals=[0.014],rotating=True):
    plt.clf()
    #nclouds = [3e1,3e2,3e3]
    nclouds = NCLOUDS
    masses = np.arange(20,121,5)
    rdm = rdmfile.RDMFile(__file__)
    plt.plot(masses,masses*0.0+1.0,color="k",alpha=0.3)
    for ncloud in nclouds[::-1]:
        ls = []
        for metal in metals:
            cloud = Cloud(ncloud,1.0*units.pc,2,FiducialDustSigma(metal))
            ratios = []
            print ("NCLOUD",ncloud)
            print(".......................")
            for starmass in masses:
                print("STARMASS",starmass)
                star = stars.Star(starmass, metal, rotating=rotating)
                star.ForceAge(REFERENCETIME)
                #r, t, dr = CalculateOutflowRadiusTime(star,cloud,rcutoff=1000*units.pc,tcutoff=1e7*units.year)
                t = 1e5*units.year
                r = dRdTforw2(star,cloud)*t         
                Prad = RadiationPressure(star,cloud,r,t)
                Pwind = WindPressure(star,cloud,r,t)
                Lwind = star.WindLuminosity(t) # From Weaver
                Lrad = star.LIonising(t) + star.LNonIonising(t)
                ci = SoundSpeed(star,t)
                #ratio = GravityPressure(star,cloud,r,t) / (Pwind + Prad)
                ratio = GravityPressure(star,cloud,r,t) / Pwind
                ratios.append(ratio)
            col = CloudColour(ncloud,metal)
            dashes = CloudLine(metal)
            label = None
            if metal == metals[0]:
                label = "$10^{"+str(int(np.log10(ncloud)))+"}$ cm$^{-3}$"
            l, = plt.plot(masses,np.array(ratios),color=col,label=label,dashes=dashes)
            rdm.AddPoints(masses,np.array(ratios),label=label)
            ls.append(l)
    leg1 = plt.legend(ncol=3,fontsize="small",frameon=False,
                    title="$n_0$, where $n(r) = n_0 (r / 1~\mathrm{pc})^2$",title_fontsize="small",
                    loc = "lower left",columnspacing=0.8)     
    leg1._legend_box.align = "left"
    leg2 = plt.legend(ls,["$Z="+str(metal)+"$" for metal in metals],
            fontsize="small",frameon=False,loc="upper right",ncol=2)
    plt.gca().add_artist(leg1)
    plt.gca().yaxis.set_ticks_position('both')
    plt.xlabel("Initial Stellar Mass / "+Msun)
    plt.yscale("log")
    plt.xlim([masses[0],masses[-1]])
    plt.ylim([4e-7,10.0])
    #plt.ylabel("$P_{grav} / (P_{w} + P_{rad})$")
    plt.ylabel("$P_{grav} / P_{w}$")
    filename = "../plots/PgravvsPbubble.pdf"
    plt.savefig(filename)
    rdm.Write(filename)

def PlotRadiationVsWindPressure(metals=[0.014],rotating=True):
    plt.clf()
    #nclouds = [3e1,3e2,3e3]
    nclouds = NCLOUDS
    masses = np.arange(20,121,5)
    rdm = rdmfile.RDMFile(__file__)
    for ncloud in nclouds[::-1]:
        ls = []
        for metal in metals:
            cloud = Cloud(ncloud,1.0*units.pc,2,FiducialDustSigma(metal))
            ratios = []
            print ("NCLOUD",ncloud)
            print(".......................")
            for starmass in masses:
                print("STARMASS",starmass)
                star = stars.Star(starmass, metal, rotating=rotating)
                star.ForceAge(REFERENCETIME)
                #r, t, dr = CalculateOutflowRadiusTime(star,cloud,rcutoff=1000*units.pc,tcutoff=1e7*units.year)
                t = 1e5*units.year
                r = dRdTforw2(star,cloud)*t         
                Prad = RadiationPressure(star,cloud,r,t)
                Pwind = WindPressure(star,cloud,r,t)
                Lwind = star.WindLuminosity(t) # From Weaver
                Lrad = star.LIonising(t) + star.LNonIonising(t)
                ci = SoundSpeed(star,t)
                ratio = Prad / Pwind
                ratios.append(ratio)
            col = CloudColour(ncloud,metal)
            dashes = CloudLine(metal)
            label = None
            if metal == metals[0]:
                label=str(ncloud)+" cm$^{-3}$"
            l, = plt.plot(masses,np.array(ratios),color=col,label=label,dashes=dashes)
            rdm.AddPoints(masses,np.array(ratios),label=label)
            ls.append(l)
    leg1 = plt.legend(ncol=2,fontsize="small",frameon=False,
                    title="$n_0$, where $n(r) = n_0 (r / 1~\mathrm{pc})^{-2}$",title_fontsize="small")     
    leg1._legend_box.align = "right"
    leg2 = plt.legend(ls,["$Z="+str(metal)+"$" for metal in metals],
            fontsize="small",frameon=False,loc="upper left",ncol=2)
    plt.gca().add_artist(leg1)
    plt.gca().yaxis.set_ticks_position('both')
    plt.xlabel("Initial Stellar Mass / "+Msun)
    plt.ylabel("$P_{rad} / P_{w} $")
    plt.yscale("log")
    plt.ylim([0.001,1.0])
    filename = "../plots/PradvsPwind.pdf"
    plt.savefig(filename)
    rdm.Write(filename)

def PlotShellThickness(ncloud,plotneutral=False):
    plt.clf()
    masses = [30,60,120]
    rdm = rdmfile.RDMFile(__file__)
    metal = 0.014
    for starmass in masses:
        rs = []
        drs = []
        ts = []
        star = stars.Star(starmass, 0.014, rotating=True)
        star.ForceAge(REFERENCETIME)
        cloud = Cloud(ncloud,1.0*units.pc,2,FiducialDustSigma(metal))
        r = 1e-3 * units.pc
        dt = 1e2 * units.year
        t = 0.0
        Cout = 0.0
        niter = 0
        while Cout < 1.0:
            t += dt
            #print(t)
            Cout = OverflowParameter(star,cloud,r,t)
            #print(Cout, r / units.pc, t / units.year / 1e3)
            Pw = WindPressure(star,cloud,r,t)
            Prad = RadiationPressure(star,cloud,r,t)
            if plotneutral:
                dr = NeutralShellThickness(star,cloud,r,t)
            else:
                dr = IonisedShellThickness(star,cloud,r,t)
            Ptot = Pw + Prad
            drdt = np.sqrt(Ptot * units.X / (cloud.nAtR(r) * units.mH))
            r += drdt * dt
            rs.append(r)
            ts.append(t)
            drs.append(dr)
            niter += 1
        #plt.plot(np.array(ts)/units.year/1e6,np.array(drs) / units.pc,label=str(ncloud)+" cm$^{-3}$")
        plt.plot(np.array(rs[1:]) / units.pc, np.array(drs[1:]) / units.pc,label=str(starmass)+" Msun")
        rdm.AddPoints(np.array(rs[1:]) / units.pc, np.array(drs[1:]) / units.pc,label=str(starmass)+" Msun")
    plt.legend()
    plt.xlabel("Radius / pc")
    plt.xscale("log")
    plt.yscale("log")
    plt.text(0.6,0.1,"$n_H = "+str(ncloud)+"$ cm$^{-3}$",transform=plt.gca().transAxes)
    if plotneutral:
        state = "neutral"
    else:
        state = "ionised"
    plt.ylabel("Thickness of "+state+" shell / pc")
    filename = "../plots/"+state+"shellthickness_ncloud"+str(ncloud)+".pdf"
    plt.savefig(filename)
    rdm.Write(filename)

def PlotShellThicknessRatio():
    plt.clf()
    masses = [30,60,120]
    ncloud = 100.0
    rdm = rdmfile.RDMFile(__file__)
    metal = 0.014
    for starmass in masses:
        rs = []
        drs = []
        ts = []
        star = stars.Star(starmass, 0.014, rotating=True)
        star.ForceAge(REFERENCETIME)
        cloud = Cloud(ncloud,1.0*units.pc,2,FiducialDustSigma(metal))
        r = 1e-2 * units.pc
        dt = 1e2 * units.year
        t = 0.0
        dr = r
        Cout = 0.0
        niter = 0
        while Cout < 1.0:
            t += dt
            #print(t)
            Cout = OverflowParameter(star,cloud,r,t)
            #print(Cout, r / units.pc, t / units.year / 1e3)
            Pw = WindPressure(star,cloud,r,t)
            Prad = RadiationPressure(star,cloud,r,t)
            dr = ShellThickness(star,cloud,r,t)
            Ptot = Pw + Prad
            drdt = np.sqrt(Ptot * units.X / (cloud.nAtR(r) * units.mH))
            r += drdt * dt
            rs.append(r)
            ts.append(t)
            drs.append(dr)
            niter += 1
        #plt.plot(np.array(ts)/units.year/1e6,np.array(drs) / units.pc,label=str(ncloud)+" cm$^{-3}$")
        plt.plot(np.array(rs) / units.pc, np.array(drs) / np.array(rs),label=str(starmass)+" Msun")
        rdm.AddPoints(np.array(rs) / units.pc, np.array(drs) / np.array(rs),label=str(starmass)+" Msun")
    plt.legend()
    plt.xlabel("Radius / pc")
    plt.ylabel("$\Delta r / r$")
    plt.xscale("log")
    plt.yscale("log")
    filename = "../plots/shellthicknessratio.pdf"
    plt.savefig(filename)
    rdm.Write(filename)

def PlotCoolingRate(metal=0.014,rotating=True):
    plt.clf()
    rdm = rdmfile.RDMFile(__file__)
    #nclouds = [3e1,3e2,3e3]
    nclouds = NCLOUDS
    #masses = np.arange(20,121,5)
    masses = [30,120]
    rs = np.logspace(-2,2,100)*units.pc
    lines = ["-","--"]
    iline = 0
    leglines = []
    for ncloud in nclouds:
        cloud = Cloud(ncloud,1.0*units.pc,2,FiducialDustSigma(metal))
        drdts = []
        drdtsnp = []
        print ("NCLOUD",ncloud)
        print(".......................")
        for starmass in masses:
            print("STARMASS",starmass)
            star = stars.Star(starmass, metal, rotating=rotating)
            star.ForceAge(REFERENCETIME)
            Lw = star.WindLuminosity(REFERENCETIME)
            drdt = dRdTforw2(star,cloud)
            ts = rs / drdt
            ys = []
            for r, t in zip(rs, ts):
                dEdt, Eb = CoolingRate(star,cloud,r,t)
                ys.append(dEdt / Lw)
            col = CloudColour(ncloud,metal)
            label = None
            if starmass == masses[0]:
                label = str(ncloud)+" cm$^{-3}$"
            line = lines[iline]
            l, = plt.plot(ts/units.year/1e6,np.array(ys),color=col,linestyle=line,label=label)
            rdm.AddPoints(ts/units.year/1e6,np.array(ys),label=label)
            if ncloud == nclouds[0]:
                leglines += [l]
            iline = (iline+1) % len(lines)
    if metal == 0.002:
        leg1 = plt.legend(ncol=1,fontsize="small",frameon=False,
                        title="$n_0$, where $n(r) = n_0 (r / 1~\mathrm{pc})^{-2}$",title_fontsize="small",loc="upper right")     
        leg1._legend_box.align = "right"
    else:
        leg2 = plt.legend(leglines,[str(starmass)+" M$_{\odot}$" for starmass in masses],
                fontsize="small",frameon=False,loc="upper right")
    #plt.gca().add_artist(leg1)
    
    plt.xlabel("Time / Myr")
    plt.ylabel("$dE_{cool}/dt~/~L_{w}$")
    plt.xscale("log")
    plt.yscale("log")
    xlims = plt.gca().get_xlim()
    ylims = plt.gca().get_ylim()
    plt.xlim([xlims[0],20.0])
    plt.ylim([1e-6,6e-1])
    plt.text(xlims[0]*1.5,3e-6,"$Z="+str(metal)+"$",fontsize="small",horizontalalignment="left",verticalalignment="top")
    filename = "../plots/coolingplot_Z"+str(metal)+".pdf"
    plt.savefig(filename)
    rdm.Write(filename)

def PlotdRdTforw2forMetals(rotating=True):
    plt.clf()
    rdm = rdmfile.RDMFile(__file__)
    #nclouds = [3e1,3e2,3e3]
    nclouds = NCLOUDS
    masses = np.arange(20,121,5)
    leglines = []
    for ncloud in nclouds:
        windratios = []
        drdts = []
        drdtslowZ = []
        print ("NCLOUD",ncloud)
        print(".......................")
        for starmass in masses:
            print("STARMASS",starmass)
            cloud = Cloud(ncloud,1.0*units.pc,2,FiducialDustSigma(0.014))
            star = stars.Star(starmass, 0.014, rotating=rotating)
            star.ForceAge(REFERENCETIME)
            wlHiZ = star.WindLuminosity(1e5 * units.year)
            drdt = dRdTforw2(star,cloud)
            drdts.append(drdt)
            cloud = Cloud(ncloud,1.0*units.pc,2,FiducialDustSigma(0.002))
            star = stars.Star(starmass, 0.002, rotating=rotating)
            star.ForceAge(REFERENCETIME)
            wlLowZ = star.WindLuminosity(1e5 * units.year)
            windratios.append(wlHiZ / wlLowZ)
            drdt = dRdTforw2(star,cloud)
            drdtslowZ.append(drdt)
        col = CloudColour(ncloud,0.014)
        l1, = plt.plot(masses,np.array(drdts)/1e5,color=col,label=str(ncloud)+" cm$^{-3}$")
        rdm.AddPoints(masses,np.array(drdts)/1e5,label=str(ncloud)+" cm$^{-3}$ Z=0.014")
        col = CloudColour(ncloud,0.002)
        l2, = plt.plot(masses,np.array(drdtslowZ)/1e5,color=col,linestyle="--")
        rdm.AddPoints(masses,np.array(drdtslowZ)/1e5,label=str(ncloud)+" cm$^{-3}$ Z=0.002")
        if len(leglines) == 0:
            leglines = [l1,l2]
        #print(np.array(drdts) / np.array(drdtslowZ))
        #print(windratios)
    leg1 = plt.legend(ncol=2,fontsize="small",frameon=False,
                    title="$n_0$, where $n(r) = n_0 (r / 1~\mathrm{pc})^{-2}$",title_fontsize="small")     
    leg1._legend_box.align = "right"
    leg2 = plt.legend(leglines,["$Z=0.014$","$Z=0.002$"],
            fontsize="small",frameon=False,loc="upper left")
    plt.gca().add_artist(leg1)
    plt.xlabel("Initial Stellar Mass / "+Msun)
    plt.ylabel("$v_2$ / km/s")
    plt.yscale("log")
    filename = "../plots/drdtforw2Metal.pdf"
    plt.savefig(filename)
    rdm.Write(filename)


def PrintStarPropertiesTable(rotating):
    # Make table intro
    rotatingtxt = "nonrotating"
    if rotating:
        rotatingtxt = "rotating"
    f = open("../plots/starpropstable_bothmetal_"+rotatingtxt+".tex","w")
    intro = r'''
    \begin{center}
    \begin{tabular}{lllllllllllllllllllll}
    \multicolumn{1}{|c|}{\textbf{log$(M_{ini}$ / M$_{\odot})$}} & 
    \multicolumn{2}{|c|}{\textbf{log$(L_{w}$ / erg/s$)$}} & 
    \multicolumn{2}{|c|}{\textbf{log$(\dot{M}_w$ / M$_{\odot}$ / yr$)$}} &
    \multicolumn{2}{|c|}{\textbf{log$(L_{n}$ / erg/s$)$}} &
    \multicolumn{2}{|c|}{\textbf{log$(L_{i}$ / erg/s$)$}} &
    \multicolumn{2}{|c|}{\textbf{log$(Q_H / $s$^{-1})$}}  &
    \multicolumn{2}{|c|}{\textbf{log$(T_i / $K$)$}} \\
    '''
    intro += r"\hline "+os.linesep+" "
    outro = r'''
    \end{tabular}
    \end{center}
    '''
    table = intro+""    
    # Fill table
    masses = np.arange(20,121,5)
    # Put the highest masses on top
    for starmass in masses[::-1]:
        print("STARMASS",starmass)
        # Solar metallicity (S)
        star = stars.Star(starmass, 0.014, rotating=rotating)
        starage = 1e5 * units.year
        star.ForceAge(starage)
        LwS = star.WindLuminosity(starage)
        MdotS = star.WindMassLoss(starage) / units.Msun * units.year
        LiS = star.LIonising(starage)
        LnS = star.LNonIonising(starage)
        QHS = star.NIonising(starage)
        TiS = Tion(star,starage)
        # Low metallicty (L)
        star = stars.Star(starmass, 0.002, rotating=rotating)
        star.ForceAge(REFERENCETIME)
        LwL = star.WindLuminosity(starage)
        MdotL = star.WindMassLoss(starage) / units.Msun * units.year
        LiL = star.LIonising(starage)
        LnL = star.LNonIonising(starage)
        QHL = star.NIonising(starage)
        TiL = Tion(star,starage)
        def fstr(val):
            v = np.log10(val)
            #return f'{np.log(val):.3}'
            return f'{v:.4}'
        line = ""
        line += str(starmass)+" & "
        line += fstr(LwS)+" & ("+fstr(LwL)+") & "
        line += fstr(MdotS)+" & ("+fstr(MdotS)+") & "
        line += fstr(LnS)+" & ("+fstr(LnL)+") & "
        line += fstr(LiS)+" & ("+fstr(LiL)+") & "
        line += fstr(QHS)+" & ("+fstr(QHL)+") & "
        line += fstr(TiS)+" & ("+fstr(TiL)+") \\\\ "+os.linesep
        table += line

    # Finish and save table
    table += outro
    f.write(table)
    f.close()

"""
-------------------------------------------------------------------------
TEST FUNCTIONS
-------------------------------------------------------------------------
"""

def OrionTest():
    # Parameters for Orion
    Mave = 1600.0*units.Msun
    Mmin = 600.0*units.Msun
    Mmax = 2600.0*units.Msun
    mstar = 30
    rotating = False
    metal = 0.014
    def zerodp(input):
        return "{:.0f}".format(input)
    def onedp(input):
        return "{:.1f}".format(input)
    def threedp(input):
        return "{:.3f}".format(input)
    def plusminus(arr,formfunc):
        vav = arr[1]
        minus = vav - np.min(arr)
        plus = np.max(arr) - vav
        return "$"+formfunc(vav)+"^{+"+formfunc(plus)+"}_{-"+formfunc(minus)+"}$"
    for mstar in [30,35]:
        print("MSTAR:", mstar)
        for ri, Omega in zip([2*units.pc,4.0*units.pc],[4.0*np.pi,1.0*np.pi]):
            vs = []
            ages = []
            ros = []
            n0s = []
            for M in [Mmin,Mave,Mmax]:
                #print("M:",M/units.Msun,"Msun")
                #n0 = 6244.0 # From Pabst+ 2019's mass of 2600 Msun
                #n0 = 4083.0 # From Pabst+ 2019's mass of 1400 Msun
                Mref = Omega*ri*(units.pc)**2*(units.mH / units.X)
                n0 = M / Mref
                n0s.append(n0)
                #print(n0,"cm^-3")
                cloud = Cloud(n0,1.0*units.pc,2,FiducialDustSigma(metal),Omega)
                star = stars.Star(mstar, 0.014, rotating=rotating)
                star.ForceAge(REFERENCETIME)
                rw = RwfromRiforw2(star,cloud,ri)
                v2 = dRdTforw2(star,cloud)
                t = rw / v2
                age = t / units.year/1e6
                drdt = dRidT(star,cloud,rw,t)
                ro, to, dro = CalculateOutflowRadiusTimew2(star,cloud)
                vs.append(drdt/1e5)
                ages.append(age)
                ros.append(ro/units.pc)
            print(onedp(ri/units.pc)+" & $"+str(int(Omega/np.pi))+" \pi$ & "+plusminus(vs,onedp)+" & "+
            plusminus(ages,threedp)+" & "+plusminus(ros,onedp)+" & "+plusminus(n0s,zerodp)+"//")
            #print(onedp(rw/units.pc),"&",
            #onedp(drdt/1e5),"km/s &",
            #threedp(age), "Myr &", 
            #onedp(ro/units.pc),"pc")

def test():
    '''
    Basic test script
    '''
    metal = 0.014
    for starmass in [30,60,120]:
        star = stars.Star(starmass, 0.014, rotating=True)
        star.ForceAge(REFERENCETIME)
        cloud = Cloud(1e3,1.0*units.pc,2,FiducialDustSigma(metal))
        drdt = dRdTforw2(star,cloud)
        tau = DraineOpticalDepth(star,1e6*units.year,1e3,FiducialDustSigma(metal))
        beta = DraineBeta(star,1e6*units.year)
        gamma = DraineGamma(star,1e6*units.year,FiducialDustSigma(metal))
        dedt = CoolingRate(star,cloud,0.1*units.pc,0.1*units.pc/drdt)
        print(drdt, tau, beta, gamma, dedt)
        #r, t = CalculateOutflowRadiusTime(star,cloud)
        #print(starmass, "Msun", r / units.pc, "pc", t / units.year / 1e6, "Myr")


if __name__=="__main__":



    # THIS IS SLOW
    for truefalse in [True, False]:
        PlotOverflowRadiusTimeNumerical(rotating=True,plotradius=truefalse)

    for metal in [0.014,0.002]:
        for truefalse in [True, False]:
            PlotMeanShellDensity(metal=metal,rotating=True,vsR = True,plotnrmsvsnmean=truefalse)
        PlotShellThicknessAnalyticvsNumerical(metal=metal,rotating=True)
        PlotShellThicknessNumerical(metal=metal,rotating=True)

    OrionTest()

    for metal in [0.014,0.002]:
        PlotShellThicknessNumerical(metal=metal,rotating=True,plotdrdt=True)

    for truefalse in [True, False]:
        PrintStarPropertiesTable(rotating=truefalse)

    PlotGravityVsBubblePressure(metals=[0.014,0.002],rotating=True)
    for metal in [0.014,0.002]:
        PlotdRdTforw2(metals=[metal],rotating=True)
        PlotCoolingRate(metal=metal,rotating=True)
    PlotOutflowTimes(metals=[0.014,0.002],rotating=True)
    PlotOutflowRadii(metals=[0.014,0.002],rotating=True)
    PlotdRdTforw2forMetals(rotating=True)
    
    PlotThicknessRatioNew(metal=0.014,rotating=True)
    PlotThicknessRatioNew(metal=0.002,rotating=True)

    PlotdRdTforSolidAngle(metal=0.014,rotating=True)
    PlotdRdTforPowerIndex(metal=0.014,rotating=True)

    OrionTest()



        #PlotShellExpansionCorrection(metal=metal,rotating=True,vsR=True)
        #PlotShellExpansionCorrection(metal=metal,rotating=True,vsR=True,ratio=True)
        #PlotThicknessCorrection(metal=metal,rotating=True,vsR=True)
        #PlotThicknessCorrection(metal=metal,rotating=True,vsR=True,uniform=True)
        #for thickshell in [True, False]:
        #    for justthickness in [True, False]:
        #        PlotThicknessRatiow2(metal=metal,rotating=True,vsR = True, justthickness=justthickness,
        #            thickshell=thickshell)
    #PlotShellThickness(10)
    #PlotThicknessRatiow2(metal=0.014,rotating=True,vsR = True)
    #PlotThicknessRatiow2(metal=0.002,rotating=True,vsR = True)
    #PlotThicknessRatiow2(metal=0.014,rotating=True,vsR = True, justthickness=True)
    #PlotThicknessRatiow2(metal=0.002,rotating=True,vsR = True, justthickness=True)
    
    #PlotRadiationVsWindPressure(metals=[0.014,0.002],rotating=True)
    #PlotMinResolutionNeeded(metals=[0.014,0.002],rotating=True)
    #PlotThicknessRatiow2(rotating=True)
    #PlotCoolingRate(metal=0.002,rotating=True)
    #PlotOutflowTimes(metals=[0.014,0.002],rotating=True)
    #PlotOutflowRadii(metals=[0.014,0.002],rotating=True)
    #PlotGravityVsBubblePressure(metals=[0.014,0.002],rotating=True)
    #PlotRadiationVsWindPressure(metals=[0.014,0.002],rotating=True)
    #PlotdRdTforw2(metals=[0.014],rotating=True)
    #PlotdRdTforw2(metals=[0.002],rotating=True)
    #PlotdRdTforw2forMetals(rotating=True)
    # THIS IS SLOW
    #for metal in [0.014,0.002]:
    #    for truefalse in [True, False]:
    #        PlotMeanShellDensity(metal=metal,rotating=True,vsR = True,plotnrmsvsnmean=truefalse)
    '''
    PlotOutflowTimes()
    PlotOutflowRadii()
    PlotdRdTforw2()
    PlotMinResolutionNeeded()
    '''
    #for neutral in [False]:
    #    PlotShellThickness(100,plotneutral=neutral)
    #    PlotShellThickness(300,plotneutral=neutral)
    #    PlotShellThickness(1000,plotneutral=neutral)