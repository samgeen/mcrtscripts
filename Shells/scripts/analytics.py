"""
Analytic solutions for shell dynamics
Sam Geen, September 2020
"""

import os, sys
sys.path.append("/home/samgeen/Programming/Astro/Weltgeist")
sys.path.append("/home/samgeen/Programming/Astro/WindInUV")

import numpy as np

import weltgeist
import weltgeist.units as wunits # make this easier to type
import stars

import matplotlib.pyplot as plt

import equilibriumtemperature

gamma = equilibriumtemperature.gamma

"""
-------------------------------------------------------------------------
GAS STUFF
-------------------------------------------------------------------------
"""
def SoundSpeed(star,t):
    #HACK - ignore change in Ti due to ionised gas density
    Ti = weltgeist.ionisedtemperatures.FindTemperature(star.Teff(t),star.metal)
    # Find ci
    # 2.0 factor is because of free electrons
    ci = np.sqrt(2.0 * gamma * wunits.kB * Ti * wunits.X / wunits.mH)
    return ci

def AlphaB(star,t):
    #HACK - ignore change in Ti due to ionised gas density
    Ti = weltgeist.ionisedtemperatures.FindTemperature(star.Teff(t),star.metal)
    return weltgeist.radiation.alpha_B_HII(Ti)

class Cloud(object):
    def __init__(self,n0,r0,w,sigmaDust,Omega=4.0*np.pi):
        self._n0 = n0
        self._r0 = r0
        self._w = w
        self._sigmaDust = sigmaDust
        self._Omega = Omega

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
        return 4.0 * np.pi / threeminusw * n0 * wunits.mH / wunits.X * r0 ** w * r**threeminusw

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
    Teff = star.Teff(timeins)
    Tion = weltgeist.radiation.IonisedGasTemperature(Teff, metal)
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
    Teff = star.Teff(timeins)
    Tion = weltgeist.radiation.IonisedGasTemperature(Teff, metal)
    gamma = 11.2 * (Tion / 1e4)**1.83 * (18 * wunits.eV / Eion) * (sigmaDust/1e-21)
    return gamma

"""
-------------------------------------------------------------------------
ANALYTIC SOLUTIONS
-------------------------------------------------------------------------
"""

def WindPressure(star,cloud,r,t,cooled=False):
    if cooled:
        raise NotImplementedError
    Omega = cloud.Omega
    vol = (Omega / 3.0) * r**3
    E = star.WindEnergy(t) * 5.0/11.0 # From Weaver
    #print ("E R VOL", E, r, vol)
    return E / vol

def RadiationPressure(star,cloud,r,t):
    L = star.LIonising(t) + star.LNonIonising(t)
    # NOTE: Omega isn't used here because light travels in straight lines
    return L / (4.0 * np.pi * r**2 * wunits.c)

def dRdTforw2(star,cloud,radiation=True):
    # Expansion rate of wind bubble for w=2
    Omega = cloud.Omega
    t= 1e5*wunits.year
    Lw = star.WindLuminosity(t)
    Lrad = 0.0
    if radiation:
        Lrad = star.LIonising(t) + star.LNonIonising(t)
    rho0 = cloud.n0 * wunits.mH /  wunits.X
    r0 = cloud.r0
    # Solve Av^3 + Bv^2 + Cv + D = 0
    A = 0.6 * Omega * rho0 * r0*r0
    B = 0.0
    C = -Omega / (4.0 * np.pi * wunits.c)*Lrad
    D = -10.0/11.0 * Lw
    roots = np.roots([A,B,C,D])
    root = roots[np.where(roots > 0.0)][0]
    #print("ROOTS", roots, root)
    return np.real(root)

def GravityPressure(star,cloud,r,t):
    n0 = cloud.n0
    r0 = cloud.r0
    w = cloud.w
    M = cloud.mAtR(r)
    Omega = cloud.Omega
    Pgrav = wunits.G*M*M / Omega / r**4
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
    RHS = QH * wunits.mH * ci ** 2 * niiVSnrms * niavVSnrms / (4.0 * np.pi * r**2 * alphaB * wunits.X * (Pw + Prad))
    #print("OVERFLOWPARAMS", LHS, RHS, QH, ci, r, alphaB, Pw, Prad)
    return RHS / LHS # if > 1, outflow happens
    
def IonisedShellThickness(star,cloud,r,t):
    n0 = cloud.n0
    r0 = cloud.r0
    w = cloud.w
    niiVSnrms = 1.4
    niavVSnrms = 1.0/1.4
    Omega = cloud.Omega # solid angle
    QH = star.NIonising(t)
    ci = SoundSpeed(star,t)
    alphaB = AlphaB(star,t)
    Pw = WindPressure(star,cloud,r,t)
    Prad = RadiationPressure(star,cloud,r,t)
    deltar = QH / (alphaB * Omega) * niiVSnrms**2.0 * (wunits.mH / wunits.X)**2.0 * ci**4.0 * (r * (Pw + Prad))**(-2.0)
    return deltar

def NeutralShellSoundSpeed(star,cloud,r,t):
    Pw = WindPressure(star,cloud,r,t)
    Prad = RadiationPressure(star,cloud,r,t)
    Ptot = Pw + Prad
    Tshell = equilibriumtemperature.NeutralTemperatureFromPressure(Ptot,star.metal)
    cshell = np.sqrt(gamma * wunits.kB * Tshell * wunits.X / wunits.mH)
    return cshell

def NeutralShellThickness(star,cloud,r,t):
    n0 = cloud.n0
    r0 = cloud.r0
    w = cloud.w
    Omega = cloud.Omega # solid angle
    drdt = ExpansionRate(star,cloud,r,t)
    deltar = (NeutralShellSoundSpeed(star,cloud,r,t) / drdt)**2.0 * r / (3.0 - w)
    return deltar

def CalculateOutflowRadiusTime(star,cloud,rcutoff = None,tcutoff = None):
    '''
    # FIRST VERSION: ITERATING OVER PRESSURE EQUATION
    r = 1e-2 * wunits.pc
    dt = 1e2 * wunits.year
    t = 0.0
    Cout = 0.0
    niter = 0
    while Cout < 1.0:
        t += dt
        #print(t)
        Cout = OverflowParameter(star,cloud,r,t)
        #print(Cout, r / wunits.pc, t / wunits.year / 1e3)
        Pw = WindPressure(star,cloud,r,t)
        Prad = RadiationPressure(star,cloud,r,t)
        Ptot = Pw + Prad
        deltar = IonisedShellThickness(star,cloud,r,t)
        drdt = np.sqrt(Ptot * wunits.X / (cloud.nAtR(r) * wunits.mH))
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
    start = 1e0*wunits.year
    end = 1e9*wunits.year
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

def CalculateOutflowRadiusTimew2(star,cloud):
    n0 = cloud.n0
    r0 = cloud.r0
    w = cloud.w
    Omega = cloud.Omega # solid angle
    t = 1e5 * wunits.year # fiducial time
    QH = star.NIonising(t)
    ci = SoundSpeed(star,t)
    alphaB = AlphaB(star,t)
    niiVSnrms = 1.4
    niavVSnrms = 1.0/1.4
    Lw = star.WindLuminosity(t)
    Lrad = star.LIonising(t) + star.LNonIonising(t)
    v2 = dRdTforw2(star,cloud)
    # Calculate the overflow radius and time if w=2 (exact solution)
    ro2 = alphaB * wunits.X * n0 * r0*r0 / (QH * wunits.mH * ci*ci * niavVSnrms * niiVSnrms)
    ro2 *= (10.0/11.0 * Lw / v2 + Lrad * Omega / (4.0 * np.pi * wunits.c))
    to2 = ro2 / v2
    deltar = IonisedShellThickness(star,cloud,ro2,to2)
    return ro2, to2, deltar

def ExpansionRate(star,cloud,r,t):
    Pw = WindPressure(star,cloud,r,t)
    Prad = RadiationPressure(star,cloud,r,t)
    dr = IonisedShellThickness(star,cloud,r,t)
    Ptot = Pw + Prad
    drdt = np.sqrt(Ptot * wunits.X / (cloud.nAtR(r) * wunits.mH))
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
    Tc = ((10.0 * Lw) / (11.0 * Omega * Cconduct * r))**(2.0/7.0)
    Eb = 10.0 * Lw * t / 11.0
    pc = Eb / (Omega * r**3)
    nc = pc / (wunits.kB * Tc)
    xcut = 1 - (Tcut / Tc)**2.5
    # Calculate the integral for the interior of the bubble
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
    nclouds = np.log10(np.array([3e1,1e2,3e2,1e3,3e3]))
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
    nclouds = [3e1,1e2,3e2,1e3,3e3]
    masses = np.arange(20,121,5)
    leglines = None
    for metal in metals:
        for ncloud in nclouds:
            cloud = Cloud(ncloud,1.0*wunits.pc,2,1e-21)
            drdts = []
            drdtsnp = []
            print ("NCLOUD",ncloud)
            print(".......................")
            for starmass in masses:
                print("STARMASS",starmass)
                star = stars.Star(starmass, metal, rotating=rotating)
                drdt = dRdTforw2(star,cloud)
                drdts.append(drdt)
                drdt = dRdTforw2(star,cloud, radiation=False)
                drdtsnp.append(drdt)
            col = CloudColour(ncloud,metal)
            dashes = CloudLine(metal)
            label = None
            if metal == metals[0]:
                label = str(ncloud)+" cm$^{-3}$"
            l1, = plt.plot(masses,np.array(drdts)/1e5,color=col,label=label,dashes=dashes)
            l2, = plt.plot(masses,np.array(drdtsnp)/1e5,color=col,linestyle=":")
            if leglines is None:
                leglines = [l1,l2]
    leg1 = plt.legend(ncol=2,fontsize="small",frameon=False,
                    title="$n_0$, where $n(r) = n_0 (r / 1~\mathrm{pc})^2$",title_fontsize="small")     
    leg1._legend_box.align = "right"
    plt.text(120,15,"$Z="+str(metals[0])+"$",fontsize="small",horizontalalignment="right",verticalalignment="top")
    leg2 = plt.legend(leglines,["Winds + Radiation Pressure","Winds Only"],
            fontsize="small",frameon=False,loc="upper left")
    plt.gca().add_artist(leg1)
    plt.gca().yaxis.set_ticks_position('both')
    plt.xlabel("Mass / Msun")
    plt.ylabel("$v_2$ / km/s")
    plt.yscale("log")
    plt.ylim([1,400])
    plt.savefig("../plots/drdtforw2_metals"+str(metals[0])+".pdf")

def PlotMinResolutionNeeded(metals=[0.014],rotating=True):
    plt.clf()
    nclouds = [3e1,1e2,3e2,1e3,3e3]
    masses = np.arange(20,121,5)
    metal = 0.014
    # NOTE: As long as w = 2, radius is unimportant
    r = 0.1*wunits.pc
    for ncloud in nclouds[::-1]:
        ls = []
        for metal in metals:
            cloud = Cloud(ncloud,1.0*wunits.pc,2,1e-21)
            minreses = []
            print ("NCLOUD",ncloud)
            print(".......................")
            for starmass in masses:
                print("STARMASS",starmass)
                star = stars.Star(starmass, metal, rotating=rotating)
                ro2, to2, deltar = CalculateOutflowRadiusTimew2(star,cloud)
                minres = MinResolutionNeeded(star,cloud,ro2,to2)
                minreses.append(minres)
            col = CloudColour(ncloud,metal)
            dashes = CloudLine(metal)
            label = None
            if metal == metals[0]:
                label = str(ncloud)+" cm$^{-3}$"
            l, = plt.plot(masses,np.array(minreses) / wunits.pc,color=col,label=label,dashes=dashes)
            ls.append(l)
    leg1 = plt.legend(ncol=2,fontsize="x-small",frameon=False,
                    title="$n_0$, where $n(r) = n_0 (r / 1~\mathrm{pc})^2$",title_fontsize="small")     
    leg1._legend_box.align = "right"
    leg2 = plt.legend(ls,["$Z="+str(metal)+"$" for metal in metals],
            fontsize="small",frameon=False,loc="upper right",ncol=2)
    plt.gca().add_artist(leg1)
    plt.xlabel("Mass / Msun")
    plt.ylabel("Minimum Resolution Needed / pc")
    plt.yscale("log")
    plt.ylim([1e-6,1e3])
    plt.savefig("../plots/minimumresolution.pdf")

def PlotThicknessRatiow2(metal=0.014,rotating=True,vsR = False,justthickness=False):
    plt.clf()
    nclouds = [3e1,3e2,3e3]
    #masses = np.arange(20,121,5)
    masses = [30,60,120]
    lines = ["-","--",":"]
    iline = 0
    leglines = []
    for ncloud in nclouds:
        cloud = Cloud(ncloud,1.0*wunits.pc,2,1e-21)
        drdts = []
        drdtsnp = []
        print ("NCLOUD",ncloud)
        print(".......................")
        for starmass in masses:
            print("STARMASS",starmass)
            star = stars.Star(starmass, metal, rotating=rotating)
            rs = np.logspace(-2,2,100)*wunits.pc
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
                deltar = IonisedShellThickness(star,cloud,r,t)
                if not justthickness:
                    ys.append(deltar/r)
                else:
                    ys.append(deltar / wunits.pc)
            col = CloudColour(ncloud, metal)
            label = None
            if starmass == masses[0]:
                label = str(ncloud)+" cm$^{-3}$"
            line = lines[iline]
            xs = ts/wunits.year/1e6
            if vsR:
                xs = rs / wunits.pc
            l, = plt.plot(xs,np.array(ys),color=col,linestyle=line,label=label)
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
    if justthickness and metal==0.014:
        leg1 = plt.legend(ncol=1,fontsize="small",frameon=False,
                        title="$n_0$, where $n(r) = n_0 (r / 1~\mathrm{pc})^2$",title_fontsize="small")     
        leg1._legend_box.align = "right"
        loc = "lower right"
        leg2 = plt.legend(leglines,[str(starmass)+" M$_{\odot}$" for starmass in masses],
                fontsize="small",frameon=False,loc=loc)
        plt.gca().add_artist(leg1)
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
    plt.ylim([1e-5,ylims[1]])
    plt.text(20,ylims[1]*0.7,"$Z="+str(metal)+"$",fontsize="small",horizontalalignment="right",verticalalignment="top")
    suffix = ""
    if vsR:
        suffix = "_vsR"
    suffix += "_Z"+str(metal)
    if not justthickness:
        plt.savefig("../plots/shellthicknessratio_w2"+suffix+".pdf")
    else:
        plt.savefig("../plots/shellthickness_w2"+suffix+".pdf") 

def PlotOutflowRadii(metals=[0.014],rotating=True):
    plt.clf()
    nclouds = [3e1,1e2,3e2,1e3,3e3]
    masses = np.arange(20,121,5)
    for ncloud in nclouds[::-1]:
        ls = []
        for metal in metals:
            cloud = Cloud(ncloud,1.0*wunits.pc,2,1e-21)
            radii = []
            print ("NCLOUD",ncloud)
            print(".......................")
            for starmass in masses:
                print("STARMASS",starmass)
                star = stars.Star(starmass, metal, rotating=rotating)
                #r, t, dr = CalculateOutflowRadiusTime(star,cloud,rcutoff=1000*wunits.pc,tcutoff=1e7*wunits.year)
                r, t, dr = CalculateOutflowRadiusTimew2(star,cloud)
                radii.append(r)
            col = CloudColour(ncloud,metal)
            dashes = CloudLine(metal)
            label = None
            if metal == metals[0]:
                label=str(ncloud)+" cm$^{-3}$"
            l, = plt.plot(masses,np.array(radii) / wunits.pc,color=col,label=label,dashes=dashes)
            ls.append(l)
    leg1 = plt.legend(ncol=2,fontsize="small",frameon=False,
                    title="$n_0$, where $n(r) = n_0 (r / 1~\mathrm{pc})^2$",title_fontsize="small")     
    leg1._legend_box.align = "left"
    #leg2 = plt.legend(ls,["$Z="+str(metal)+"$" for metal in metals],
    #        fontsize="small",frameon=False,loc="lower right",ncol=2)
    plt.gca().add_artist(leg1)
    plt.gca().yaxis.set_ticks_position('both')
    plt.xlabel("Mass / Msun")
    plt.ylabel("Radius of Overflow $r_{o,2}$ / pc")
    plt.yscale("log")
    plt.ylim([1e-4,100])
    plt.savefig("../plots/outflowradii.pdf")

def PlotOutflowTimes(metals=[0.014],rotating=True):
    plt.clf()
    nclouds = [3e1,1e2,3e2,1e3,3e3]
    masses = np.arange(20,121,5)
    for ncloud in nclouds[::-1]:
        ls = []
        for metal in metals:
            cloud = Cloud(ncloud,1.0*wunits.pc,2,1e-21)
            times = []
            print ("NCLOUD",ncloud)
            print(".......................")
            for starmass in masses:
                print("STARMASS",starmass)
                star = stars.Star(starmass, metal, rotating=rotating)
                #r, t, dr = CalculateOutflowRadiusTime(star,cloud,rcutoff=1000*wunits.pc,tcutoff=1e7*wunits.year)
                r, t, dr = CalculateOutflowRadiusTimew2(star,cloud)
                times.append(t)
            col = CloudColour(ncloud,metal)
            dashes = CloudLine(metal)
            label = None
            if metal == metals[0]:
                label=str(ncloud)+" cm$^{-3}$"
            l, = plt.plot(masses,np.array(times) / wunits.year / 1e6,color=col,label=label,dashes=dashes)
            ls.append(l)
    #leg1 = plt.legend(ncol=2,fontsize="small",frameon=False,
    #                title="$n_0$, where $n(r) = n_0 (r / 1~\mathrm{pc})^2$",title_fontsize="small")     
    #leg1._legend_box.align = "right"
    leg2 = plt.legend(ls,["$Z="+str(metal)+"$" for metal in metals],
            fontsize="small",frameon=False,loc="upper right",ncol=2)
    #plt.gca().add_artist(leg1)
    plt.gca().yaxis.set_ticks_position('both')
    plt.xlabel("Mass / Msun")
    plt.ylabel("Time of Overflow  $t_{o,2}$ / Myr")
    plt.yscale("log")
    plt.ylim([1e-4,30])
    plt.savefig("../plots/outflowtimes.pdf") 

def PlotGravityVsBubblePressure(metals=[0.014],rotating=True):
    plt.clf()
    nclouds = [3e1,1e2,3e2,1e3,3e3]
    masses = np.arange(20,121,5)
    for ncloud in nclouds[::-1]:
        ls = []
        for metal in metals:
            cloud = Cloud(ncloud,1.0*wunits.pc,2,1e-21)
            ratios = []
            print ("NCLOUD",ncloud)
            print(".......................")
            for starmass in masses:
                print("STARMASS",starmass)
                star = stars.Star(starmass, metal, rotating=rotating)
                #r, t, dr = CalculateOutflowRadiusTime(star,cloud,rcutoff=1000*wunits.pc,tcutoff=1e7*wunits.year)
                t = 1e5*wunits.year
                r = dRdTforw2(star,cloud)*t         
                Prad = RadiationPressure(star,cloud,r,t)
                Pwind = WindPressure(star,cloud,r,t)
                Lwind = star.WindLuminosity(t) # From Weaver
                Lrad = star.LIonising(t) + star.LNonIonising(t)
                ci = SoundSpeed(star,t)
                ratio = GravityPressure(star,cloud,r,t) / (WindPressure(star,cloud,r,t) + RadiationPressure(star,cloud,r,t))
                ratios.append(ratio)
            col = CloudColour(ncloud,metal)
            dashes = CloudLine(metal)
            label = None
            if metal == metals[0]:
                label=str(ncloud)+" cm$^{-3}$"
            l, = plt.plot(masses,np.array(ratios),color=col,label=label,dashes=dashes)
            ls.append(l)
    #leg1 = plt.legend(ncol=2,fontsize="small",frameon=False,
    #                title="$n_0$, where $n(r) = n_0 (r / 1~\mathrm{pc})^2$",title_fontsize="small")     
    #leg1._legend_box.align = "right"
    leg2 = plt.legend(ls,["$Z="+str(metal)+"$" for metal in metals],
            fontsize="small",frameon=False,loc="upper right",ncol=2)
    #plt.gca().add_artist(leg1)
    plt.gca().yaxis.set_ticks_position('both')
    plt.xlabel("Mass / Msun")
    plt.ylabel("$P_{rad} / P_{w} $")
    plt.yscale("log")
    plt.ylim([1e-6,1.0])
    plt.ylabel("$P_{grav} / (P_{w} + P_{rad})$")
    plt.savefig("../plots/PgravvsPbubble.pdf")

def PlotRadiationVsWindPressure(metals=[0.014],rotating=True):
    plt.clf()
    nclouds = [3e1,1e2,3e2,1e3,3e3]
    masses = np.arange(20,121,5)
    for ncloud in nclouds[::-1]:
        ls = []
        for metal in metals:
            cloud = Cloud(ncloud,1.0*wunits.pc,2,1e-21)
            ratios = []
            print ("NCLOUD",ncloud)
            print(".......................")
            for starmass in masses:
                print("STARMASS",starmass)
                star = stars.Star(starmass, metal, rotating=rotating)
                #r, t, dr = CalculateOutflowRadiusTime(star,cloud,rcutoff=1000*wunits.pc,tcutoff=1e7*wunits.year)
                t = 1e5*wunits.year
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
            ls.append(l)
    leg1 = plt.legend(ncol=2,fontsize="small",frameon=False,
                    title="$n_0$, where $n(r) = n_0 (r / 1~\mathrm{pc})^2$",title_fontsize="small")     
    leg1._legend_box.align = "right"
    leg2 = plt.legend(ls,["$Z="+str(metal)+"$" for metal in metals],
            fontsize="small",frameon=False,loc="upper left",ncol=2)
    plt.gca().add_artist(leg1)
    plt.gca().yaxis.set_ticks_position('both')
    plt.xlabel("Mass / Msun")
    plt.ylabel("$P_{rad} / P_{w} $")
    plt.yscale("log")
    plt.ylim([0.001,1.0])
    plt.savefig("../plots/PradvsPwind.pdf")

def PlotShellThickness(ncloud,plotneutral=False):
    plt.clf()
    masses = [30,60,120]
    for starmass in masses:
        rs = []
        drs = []
        ts = []
        star = stars.Star(starmass, 0.014, rotating=True)
        cloud = Cloud(ncloud,1.0*wunits.pc,2,1e-21)
        r = 1e-3 * wunits.pc
        dt = 1e2 * wunits.year
        t = 0.0
        Cout = 0.0
        niter = 0
        while Cout < 1.0:
            t += dt
            #print(t)
            Cout = OverflowParameter(star,cloud,r,t)
            #print(Cout, r / wunits.pc, t / wunits.year / 1e3)
            Pw = WindPressure(star,cloud,r,t)
            Prad = RadiationPressure(star,cloud,r,t)
            if plotneutral:
                dr = NeutralShellThickness(star,cloud,r,t)
            else:
                dr = IonisedShellThickness(star,cloud,r,t)
            Ptot = Pw + Prad
            drdt = np.sqrt(Ptot * wunits.X / (cloud.nAtR(r) * wunits.mH))
            r += drdt * dt
            rs.append(r)
            ts.append(t)
            drs.append(dr)
            niter += 1
        #plt.plot(np.array(ts)/wunits.year/1e6,np.array(drs) / wunits.pc,label=str(ncloud)+" cm$^{-3}$")
        plt.plot(np.array(rs[1:]) / wunits.pc, np.array(drs[1:]) / wunits.pc,label=str(starmass)+" Msun")
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
    plt.savefig("../plots/"+state+"shellthickness_ncloud"+str(ncloud)+".pdf")

def PlotShellThicknessRatio():
    plt.clf()
    masses = [30,60,120]
    ncloud = 100.0
    for starmass in masses:
        rs = []
        drs = []
        ts = []
        star = stars.Star(starmass, 0.014, rotating=True)
        cloud = Cloud(ncloud,1.0*wunits.pc,2,1e-21)
        r = 1e-2 * wunits.pc
        dt = 1e2 * wunits.year
        t = 0.0
        dr = r
        Cout = 0.0
        niter = 0
        while Cout < 1.0:
            t += dt
            #print(t)
            Cout = OverflowParameter(star,cloud,r,t)
            #print(Cout, r / wunits.pc, t / wunits.year / 1e3)
            Pw = WindPressure(star,cloud,r,t)
            Prad = RadiationPressure(star,cloud,r,t)
            dr = ShellThickness(star,cloud,r,t)
            Ptot = Pw + Prad
            drdt = np.sqrt(Ptot * wunits.X / (cloud.nAtR(r) * wunits.mH))
            r += drdt * dt
            rs.append(r)
            ts.append(t)
            drs.append(dr)
            niter += 1
        #plt.plot(np.array(ts)/wunits.year/1e6,np.array(drs) / wunits.pc,label=str(ncloud)+" cm$^{-3}$")
        plt.plot(np.array(rs) / wunits.pc, np.array(drs) / np.array(rs),label=str(starmass)+" Msun")
    plt.legend()
    plt.xlabel("Radius / pc")
    plt.ylabel("$\Delta r / r$")
    plt.xscale("log")
    plt.yscale("log")
    plt.savefig("../plots/shellthicknessratio.pdf")

def PlotCoolingRate(metal=0.014,rotating=True):
    plt.clf()
    nclouds = [1e1,1e2,1e3]
    #masses = np.arange(20,121,5)
    masses = [30,60,120]
    rs = np.logspace(-2,2,100)*wunits.pc
    lines = ["-","--",":"]
    iline = 0
    leglines = []
    for ncloud in nclouds:
        cloud = Cloud(ncloud,1.0*wunits.pc,2,1e-21)
        drdts = []
        drdtsnp = []
        print ("NCLOUD",ncloud)
        print(".......................")
        for starmass in masses:
            print("STARMASS",starmass)
            star = stars.Star(starmass, metal, rotating=rotating)
            drdt = dRdTforw2(star,cloud)
            ts = rs / drdt
            ys = []
            for r, t in zip(rs, ts):
                dEdt, Eb = CoolingRate(star,cloud,r,t)
                print(Eb, dEdt)
                Lw = star.WindLuminosity(t)
                ys.append(dEdt / Lw)
            col = CloudColour(ncloud,metal)
            label = None
            if starmass == masses[0]:
                label = str(ncloud)+" cm$^{-3}$"
            line = lines[iline]
            l, = plt.plot(ts/wunits.year/1e6,np.array(ys),color=col,linestyle=line,label=label)
            if ncloud == nclouds[0]:
                leglines += [l]
            iline = (iline+1) % len(lines)
    leg1 = plt.legend(ncol=1,fontsize="small",frameon=False,
                    title="$n_0$, where $n(r) = n_0 (r / 1~\mathrm{pc})^2$",title_fontsize="small",loc="upper right")     
    leg1._legend_box.align = "right"
    leg2 = plt.legend(leglines,[str(starmass)+" M$_{\odot}$" for starmass in masses],
            fontsize="small",frameon=False,loc="upper left")
    plt.gca().add_artist(leg1)
    
    plt.xlabel("Time / Myr")
    plt.ylabel("$dE_{cool}/dt~/~L_{w}$")
    plt.xscale("log")
    plt.yscale("log")
    xlims = plt.gca().get_xlim()
    ylims = plt.gca().get_ylim()
    plt.xlim([xlims[0],20.0])
    plt.ylim([2e-7,3e-1])
    plt.text(xlims[0]*1.5,1e-6,"$Z="+str(metal)+"$",fontsize="small",horizontalalignment="left",verticalalignment="top")
    plt.savefig("../plots/coolingplot_Z"+str(metal)+".pdf")

def PlotdRdTforw2forMetals(rotating=True):
    plt.clf()
    nclouds = [1e1,3e1,1e2,3e2,1e3]
    masses = np.arange(20,121,5)
    for ncloud in nclouds:
        windratios = []
        drdts = []
        drdtslowZ = []
        print ("NCLOUD",ncloud)
        print(".......................")
        for starmass in masses:
            print("STARMASS",starmass)
            cloud = Cloud(ncloud,1.0*wunits.pc,2,1e-21)
            star = stars.Star(starmass, 0.014, rotating=rotating)
            wlHiZ = star.WindLuminosity(1e5 * wunits.year)
            drdt = dRdTforw2(star,cloud)
            drdts.append(drdt)
            cloud = Cloud(ncloud,1.0*wunits.pc,2,1e-21*0.002/0.014)
            star = stars.Star(starmass, 0.002, rotating=rotating)
            wlLowZ = star.WindLuminosity(1e5 * wunits.year)
            windratios.append(wlHiZ / wlLowZ)
            drdt = dRdTforw2(star,cloud)
            drdtslowZ.append(drdt)
        col = CloudColour(ncloud,0.014)
        plt.plot(masses,np.array(drdts)/1e5,color=col,label=str(ncloud)+" cm$^{-3}$")
        col = CloudColour(ncloud,0.002)
        plt.plot(masses,np.array(drdtslowZ)/1e5,color=col)
        print(np.array(drdts) / np.array(drdtslowZ))
        print(windratios)
    leg = plt.legend(ncol=2,fontsize="small",frameon=False,
                    title="$n_0$, where $n(r) = n_0 (r / 1~\mathrm{pc})^2$",title_fontsize="small")     
    leg._legend_box.align = "right"
    plt.xlabel("Mass / Msun")
    plt.ylabel("$v_2$ / km/s")
    plt.yscale("log")
    plt.savefig("../plots/drdtforw2Metal.pdf")


"""
-------------------------------------------------------------------------
TEST FUNCTIONS
-------------------------------------------------------------------------
"""

def OrionTest():
    # Parameters for Orion
    Mave = 1600.0*wunits.Msun
    Mmin = 600.0*wunits.Msun
    Mmax = 2600.0*wunits.Msun
    mstar = 30
    rotating = False
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
        for rw, Omega in zip([2*wunits.pc,4.0*wunits.pc],[4.0*np.pi,1.0*np.pi]):
            vs = []
            ages = []
            ros = []
            n0s = []
            for M in [Mmin,Mave,Mmax]:
                #print("M:",M/wunits.Msun,"Msun")
                #n0 = 6244.0 # From Pabst+ 2019's mass of 2600 Msun
                #n0 = 4083.0 # From Pabst+ 2019's mass of 1400 Msun
                Mref = Omega*rw*(wunits.pc)**2*(wunits.mH / wunits.X)
                n0 = M / Mref
                n0s.append(n0)
                #print(n0,"cm^-3")
                cloud = Cloud(n0,1.0*wunits.pc,2,1e-21,Omega)
                star = stars.Star(mstar, 0.014, rotating=rotating)
                drdt = dRdTforw2(star,cloud)
                ro, to, dro = CalculateOutflowRadiusTimew2(star,cloud)
                age = rw / drdt / wunits.year/1e6
                vs.append(drdt/1e5)
                ages.append(age)
                ros.append(ro/wunits.pc)
            print(onedp(rw/wunits.pc)+" pc & $"+str(int(Omega/np.pi))+" \pi$ & "+plusminus(vs,onedp)+" km/s & "+
            plusminus(ages,threedp)+" Myr & "+plusminus(ros,onedp)+" pc & "+plusminus(n0s,zerodp)+"//")
            #print(onedp(rw/wunits.pc),"&",
            #onedp(drdt/1e5),"km/s &",
            #threedp(age), "Myr &", 
            #onedp(ro/wunits.pc),"pc")

def test():
    '''
    Basic test script
    '''
    for starmass in [30,60,120]:
        star = stars.Star(starmass, 0.014, rotating=True)
        cloud = Cloud(1e3,1.0*wunits.pc,2,1e-21)
        drdt = dRdTforw2(star,cloud)
        tau = DraineOpticalDepth(star,1e6*wunits.year,1e3,1e-21)
        beta = DraineBeta(star,1e6*wunits.year)
        gamma = DraineGamma(star,1e6*wunits.year,1e-21)
        dedt = CoolingRate(star,cloud,0.1*wunits.pc,0.1*wunits.pc/drdt)
        print(drdt, tau, beta, gamma, dedt)
        #r, t = CalculateOutflowRadiusTime(star,cloud)
        #print(starmass, "Msun", r / wunits.pc, "pc", t / wunits.year / 1e6, "Myr")


if __name__=="__main__":
    #PlotShellThickness(10)
    PlotMinResolutionNeeded(metals=[0.014,0.002],rotating=True)
    PlotThicknessRatiow2(metal=0.014,rotating=True,vsR = True)
    PlotThicknessRatiow2(metal=0.002,rotating=True,vsR = True)
    PlotThicknessRatiow2(metal=0.014,rotating=True,vsR = True, justthickness=True)
    PlotThicknessRatiow2(metal=0.002,rotating=True,vsR = True, justthickness=True)
    PlotThicknessRatiow2(rotating=True)
    PlotCoolingRate(metal=0.014,rotating=True)
    PlotCoolingRate(metal=0.002,rotating=True)
    PlotOutflowTimes(metals=[0.014,0.002],rotating=True)
    PlotOutflowRadii(metals=[0.014,0.002],rotating=True)
    PlotGravityVsBubblePressure(metals=[0.014,0.002],rotating=True)
    PlotRadiationVsWindPressure(metals=[0.014,0.002],rotating=True)
    PlotdRdTforw2(metals=[0.014],rotating=True)
    PlotdRdTforw2(metals=[0.002],rotating=True)
    PlotCoolingRate(rotating=True)
    PlotdRdTforw2forMetals(rotating=True)
    OrionTest()
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