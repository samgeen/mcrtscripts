'''
Compare to the Raga solution (c.f. Pascal Tremblin's work)
Sam Geen, April 2015
'''

import os

import customplot
import matplotlib.pyplot as plt
import matplotlib.lines as mlines

import scipy.interpolate

import Hamu

import radii as radiimodule # oops bad naming

import numpy as np

import profilesphere, outflowmodel, plottimeseries, timeplot, rayprof
import testprofilefit

from pymses.utils import constants as C

import linestyles

# Some global properties
# BE VERY CAREFUL TO MAKE SURE THAT THESE ARE GOOOOOOOOD!
# Embedded -> star at centre of cloud (embedded = True)
# Blister -> star at edge of cloud (embedded = False)
Te = 8400.0 # K
Text = 50.0 # K, I dunno
kB = 1.3806e-16 # erg / K
gamma = 1.4
X = 0.74
#mu = X*2.0 + 0.25*(1-X)*2.0 # Ionised hydrogen plus once-ionised He
mu = 0.61 # From Matzner 2002
mH = 1.67e-24 # g
cs = np.sqrt(gamma * kB * Te / (mH*mu))
c0 = np.sqrt(gamma * kB * Text / (mH*mu))
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

class ProfAna(object):
    def __init__(self,sim):
        self._n0 = outflowmodel.FindnH(sim)
        rhofit,w = testprofilefit.Run(sim.Name(),outflowmodel.tstart)
        r0 = np.exp((np.log(rhofit) - np.log(self._n0))/w)
        r0 *= outflowmodel.pcincm
        self._r0w = r0**w
        self._minusw = -w

    def rho(self,r):
        #if r < self._r0w**(-1.0/self._minusw):
        #    return self._n0
        #else:
        if r == 0:
            return self._n0
        return self._n0 * r**self._minusw * self._r0w

    def vesc(self,r):
        #if self._minusw < -3.0:
        #    import pdb; pdb.set_trace()
        mp = outflowmodel.mH/outflowmodel.X
        vesc = np.sqrt(8.0*np.pi * G * self._n0*mp * self._r0w * \
                           r**(2.0+self._minusw) / (3.0+self._minusw))
        vesc /= 1e5 # expects km/s
        return -vesc

class DensAna(object):
    def __init__(self,profobj):
        self._prof = profobj

    def __call__(self,r,t):
        # HACK
        #if r < 3*outflowmodel.pcincm:
        #    return [500.0]
        #else:
        #    return [1.0]
        return [self._prof.rho(r)]

    def LastTime(self):
        return 5.0*Myrins

class VelAna(object):
    def __init__(self,profobj):
        self._prof = profobj

    def __call__(self,r,t):
        return [self._prof.vesc(r)]

    def LastTime(self):
        return 5.0*Myrins

def smooth(x,window_len=11,window='hanning'):
 
    if x.ndim != 1:
        raise ValueError, "smooth only accepts 1 dimension arrays."
        
    if x.size < window_len:
        raise ValueError, "Input vector needs to be bigger than window size."
  
    if window_len<3:
        return x
        
  
    if not window in ['flat', 'hanning', 'hamming', 'bartlett', 'blackman']:
        raise ValueError, "Window is on of 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'"
  
   
    s=np.r_[x[window_len-1:0:-1],x,x[-1:-window_len:-1]]
    #print(len(s))
    if window == 'flat': #moving average
        w=np.ones(window_len,'d')
    else:
        w=eval('np.'+window+'(window_len)')
  
    y=np.convolve(w/w.sum(),s,mode='valid')
    return y

class Profiles(object):
    def __init__(self, sim,hydro,profmodule=profilemodule):
        self._sim = sim
        self._proffunc = None
        self._profmodule = profmodule
        self._lasttime = 0.0
        self._hydro = hydro
        self._Setup()
        
    def _Setup(self):
        # Gather profiles
        profs = []
        r = None
        t = []
        myr = None
        profs = None
        rs = None
        ts = []
        nt = len(self._sim.Times())
        times = np.array(self._sim.Times())
        itime = 0
        for snap in self._sim.Snapshots():
            if myr is None:
                myr = snap.RawData().info["unit_time"].express(C.Myr)
                times *= myr*Myrins
            if self._profmodule.__name__=="rayprof":
                r, prof = self._profmodule.medianprofileHamu(snap,self._hydro)
            else:
                r, prof = self._profmodule.profileHamu(snap,self._hydro,1e6)
            #r, prof = rayprof.medianprofileHamu(snap,self._hydro) - bad
            for i in range(0,len(prof)):
                if prof[i] == 0.0:
                    try:
                        prof[i] = prof[i+1]
                    except:
                        prof[i] = prof[i-1]
            if profs is None:
                nr = len(r)
                profs = np.zeros((nr,nt))
            profs[:,itime] = prof
            itime += 1
        rs = r
        ts = times
        #import pdb; pdb.set_trace()
        #nt = len(ts)
        #for i in range(0,nr):
        #    profs[i,:] = smooth(profs[i,:],window_len=10)[0:nt]
        rs *= pcincm
        self._lasttime = times[-1]
        # Make interpolator
        # TODO: CHECK CORRECT x, y AXIS ORIENTATION
        self._proffunc = scipy.interpolate.interp2d(rs,ts,profs.T,kind="linear")
        #import pdb; pdb.set_trace()
        #self._proffunc = scipy.interpolate.RectBivariateSpline(rs,ts,profs)

    def LastTime(self):
        return self._lasttime

    def __call__(self, radii, times):
        return self._proffunc(radii, times)

class SolutionFull(object):
    '''
    Full Raga + infall solution
    '''
    def __init__(self,n0,rs,rcloud):
        self._n0 = n0
        self._rs = rs
        self._rcloud = rcloud

    def fracs(self,r,n,vr,turb):
        '''
        Compute fractions from inputs
        '''
        vdispsq = (turb*1e5)**2 + c0**2
        return r/self._rs, n/self._n0, vr*1e5/cs, vdispsq/cs**2

    def __call__(self,r,n,vr,turb):
        rs = self._rs
        n0 = self._n0
        #import pdb; pdb.set_trace()
        vdispsq = (turb*1e5)**2 + c0**2
        #print "vs:", turb*1e5, vr*1e5,c0,cs
        #print "rs:", r, rs
        #print "ns:", n, n0
        #vdispsq = c0**2
        # OUTFLOW HACK
        #if r > self._rcloud:
        #    return cs
        drdt = cs*((rs/r)**(3./4.)*np.sqrt(n0/n)-vdispsq/cs**2/((rs/r)**(3./4.)*np.sqrt(n0/n))+vr*1e5/cs)
        return drdt

class SolutionSpitzerLike(object):
    '''
    Just the Spitzer-like parts
    '''
    def __init__(self,n0,rs,rcloud):
        self._n0 = n0
        self._rs = rs
        self._rcloud = rcloud
        print "n0, rs", self._n0, self._rs

    def fracs(self,r,n,vr,turb):
        '''
        Compute fractions from inputs
        '''
        vdispsq = (turb*1e5)**2 + c0**2
        return r/self._rs, n/self._n0, vr*1e5/cs, vdispsq/cs**2

    def __call__(self,r,n,vr,turb):
        rs = self._rs
        n0 = self._n0
        #if r > self._rcloud:
        #    return cs
        return cs*((rs/r)**(3./4.)*np.sqrt(n0/n))

def ComputeExpansion(sim,solntype="full",proftype="sim"):
    '''
    Compute the solution by integrating in time
    solntype - do "spitzlike" (Spitzer 1978 soln) or "full" (raga + infall)
    proftype - "sim" (sampled from simulation) or "ana" (simple analytic)
    '''
    simname = sim.Name()
    fluxstr = simname[1:3]
    tend = 5.0
    tini = 1.25
    if not "noturb" in simname:
        nofluxsim = Hamu.Simulation(simname.replace(fluxstr,"00"))
    else:
        nofluxsim = Hamu.Simulation("noturb00")
        tini = 0.0
    rcloud = 3.0*pcincm
    if "F2" in simname:
        tini = 1.25*2.0
        nofluxsim = Hamu.Simulation("N00_M4_B02")
    elif "F3" in simname:
        tini = 1.25*3.0
        nofluxsim = Hamu.Simulation("N00_M4_B02")
        tend = 7.0
    elif "_C2" in simname:
        rcloud *= 0.75**2
        tini = 1.25*0.75**3
    elif "_C" in simname:
        rcloud *= 0.5**2
        tini = 1.25*0.5**3
        nofluxsim = Hamu.Simulation("N00_M4_B02_C8")
    outflowmodel.tstart = tini
    timeplot.starts[simname] = tini

    flux = outflowmodel.FindFlux(sim)
    
    # Get the simulation profiles
    nprof = Profiles(nofluxsim,"rho")#,profmodule=rayprof)
    vrprof = Profiles(nofluxsim,"vrad")
    spdprof = Profiles(nofluxsim,"spd")

    # Compute the Stromgren radius in the density field
    numr = 100000
    dr = 1e-5*pcincm # some small value
    rcurr = 0.0
    rs = 0.0
    nrecomb = 0.0
    beta2 = outflowmodel.beta2
    for i in range(0,numr):
        ncurr = nprof(rs,0)
        nrecomb += beta2 * ncurr ** 2 * rs ** 2 * 4 * np.pi * dr
        if nrecomb < flux:
            rs += dr
        else:
            break
    n0 = np.sqrt(flux / (beta2 * 4.0/3.0 * np.pi * (rs)**3))
    print rs / pcincm
    print n0

    # Use analytic profiles instead?
    if proftype == "ana":
        #r0 = outflowmodel.FindScaleRadius(sim,n0)
        #w = -outflowmodel.FindProfilePower(sim)
        profana = ProfAna(sim)
        nprof = DensAna(profana)
        vrprof = VelAna(profana)
        spdprof = lambda r,t: [0.0] # Just set to zero for now

    # HACK - CHECK THAT THIS WORKS OK!!!
    #n0 = outflowmodel.FindnH(sim)
    #rs = outflowmodel.Findrstromgren(sim,dens=n0)
    # Set function to use
    if solntype == "full":
        soln = SolutionFull(n0,rs,rcloud)
    elif solntype == "spitzlike":
        soln = SolutionSpitzerLike(n0,rs,rcloud)
    else:
        print "Oops, we don't have this type of solution:", solntype
        raise ValueError
    
    # Run integrator
        
    dt = 1e-4
    numr = int((tend-tini)/dt+1)
    rii = np.zeros((numr))
    times = np.zeros((numr))
    rfracs = np.zeros((numr))
    nfracs = np.zeros((numr))
    vrfracs = np.zeros((numr))
    vtfracs = np.zeros((numr))
    denss = np.zeros((numr))
    vrs = np.zeros((numr))
    rii[0] = rs
    times[0] = tini*Myrins
    tend *= Myrins
    dt *= Myrins
    itime = 0
    print numr
    while itime < numr-1: #(times[itime] < tend):
        # Sample simulation values
        nnext = nprof(rii[itime],times[itime])[0]
        if np.isnan(nnext) or nnext < 0.0:
            nnext = ncurr
        ncurr = nnext
        vrcurr = vrprof(rii[itime],times[itime])[0]
        denss[itime] = ncurr
        vrs[itime] = vrcurr
        spdcurr = spdprof(rii[itime],times[itime])[0]
        spdcent = spdprof(0,times[itime])[0]
        #turbcurr = spdcent
        #turbcurr = spdcurr
        #turbcurr = np.sqrt(np.max([(spdcurr**2 - vrcurr**2,0.0)]))
        turbcurr = 0.0
        #if times[itime] > 3.0*tini*Myrins and solntype == "full":
        #    if rii[itime] < rcloud:
        #turbcurr = spdcent
                #soln = SolutionSpitzerLike(n0,rs,rcloud)
                #solntype = "outflow"
        # HACK
        #turbcurr = spdcent
        #vrcurr = 0.0
        # Compute solutions
        drdt = soln(rii[itime],ncurr,vrcurr,turbcurr)
        if np.isnan(drdt) or np.isinf(drdt):
           drdt = 0.0 
        rii[itime+1] = rii[itime]+drdt*dt
        if rii[itime+1] < 0:
            rii[itime+1] = 0.0
        # Compute fractional values from model inputs
        rfracs[itime],nfracs[itime],vrfracs[itime],vtfracs[itime] = \
            soln.fracs(rii[itime],ncurr,vrcurr,turbcurr)
        # Update time
        times[itime+1] = times[itime]+dt
        itime += 1

    # Cut solution after last time in profile data
    intime = times < nprof.LastTime()
    rii = rii[intime]
    rfracs = rfracs[intime]
    nfracs = nfracs[intime]
    vrfracs = vrfracs[intime]
    vtfracs = vtfracs[intime]
    denss = denss[intime]
    vrs = vrs[intime]
    times = times[intime]
    # Scale outputs
    rii /= pcincm
    times /= Myrins
    times -= tini
    return times, rii, rfracs, nfracs, vrfracs, vtfracs, denss, vrs

def PlotForSims(sims,name,proftype="ana"):
    plt.clf()
    #cols = ["r","m","b"]
    cols = []
    for sim in sims:
        cols.append(linestyles.col(sim.Name()))
    try:
        os.mkdir("../plots/raga/"+name)
    except:
        pass
    solnsfull = {}
    solnsspitz = {}
    icol = 0 
    for sim in sims:
        # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        # Run this first to set the tstart values in timeplot
        # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        solnsfull[sim.Name()]  = ComputeExpansion(sim,"full",proftype)
        solnsspitz[sim.Name()] = ComputeExpansion(sim,"spitzlike",proftype)
    icol = 0 
    for sim in sims:
        col = cols[icol]
        #tsim, rsim = timeplot.funcovertime(sim,plottimeseries.meanradius)
        # HACK? - use median radius not mass-weighted mean
        tsim, rsim = timeplot.funcovertime(sim,radiimodule.MedianRadiusProfile)
        plt.plot(tsim,rsim,color=col,
                 linestyle="-",label=sim.Label())#+" sim")
        icol += 1
    icol = 0
    for sim in sims:
        col = cols[icol]
        # d[1-4] are dummy variables we don't need here
        traga, rraga, d1, d2, d3, d4, d5, d6 = solnsspitz[sim.Name()]
        #plt.plot(traga,rraga,color=col,
        #         linestyle=":")#,label="FE") # Free Expansion
        icol += 1
    icol = 0
    for sim in sims:
        col = cols[icol]
        # d[1-4] are dummy variables we don't need here
        traga, rraga, d1, d2, d3, d4, d5, d6 = solnsfull[sim.Name()]
        plt.plot(traga,rraga,color=col,
                 linestyle="--")#,label="IL") # Infall Limited
        icol += 1
    if name != "freefall":
        plt.xlim([0,4])
    else:
        plt.xlim([0,1])
    if name == "compact":
        plt.ylim([0.01,14])
    elif name == "freefall":
        plt.ylim([0.2,5])
    else:
        plt.ylim([0.6,14])
    anatxt = "Sampled"
    if proftype == "ana":
        anatxt = "Power Law"
    if name != "freefall":
        plt.text(0.1,10,anatxt)
    plt.yscale("log")
    plt.xlabel("Time / Myr")
    plt.ylabel("Radius / pc")
    ncol = 1
    loc = (0.14,0.05)
    if name == "compact":
        loc = "lower right"
        ncol = 1
    leg = None
    leg2 = None
    if proftype == "ana" or name == "freefall":
        leg = plt.legend(fontsize="small",loc=loc,
                         ncol=ncol,frameon=False)
    if name == "compact" and proftype != "ana" or name == "freefall":
        sline = mlines.Line2D([], [], color='k', label='Simulation')
        #fline = mlines.Line2D([], [], color='k', linestyle=":",
        #                      label='Free Expansion')
        #iline = mlines.Line2D([], [], color='k', linestyle="--",
        #                      label='Infall-Limited')
        iline = mlines.Line2D([], [], color='k', linestyle="--",
                              label='Non-Static Model')
        leg2 = plt.legend(handles=[sline,iline],ncol=1,fontsize="small",
                          frameon=False,loc="lower right")
    if leg is not None and leg2 is not None:
        plt.gca().add_artist(leg)


    # Re-align the y-axes THANK YOU VERY MUCH FOR DOING THIS YOURSELF PYPLOT
    #txts = leg.get_texts()
    #nsim = len(sims)
    #for itxt in range(0,nsim):
    #    left = txts[itxt]
    #    right = txts[itxt+nsim]
    #    right.set_y(left.get_position()[1])
    #    left.set_verticalalignment("baseline")
    #    right.set_verticalalignment("baseline")
    plt.savefig("../plots/raga/"+name+"/ragacompare_"+proftype+".pdf")
    # Now plot the fractional values
    plt.clf()
    solns = solnsspitz
    iline = 0
    icol = 0
    line = ["-","--",":","-."]
    label = sim.Label()
    '''
    for sim in sims:
        col = cols[icol]
        traga,rraga,rfrac,nfrac,vrfrac,vtfrac = solns[sim.Name()]
        plt.plot(traga,rfrac,col+line[iline],label=label+" r/r$_s$")
        icol += 1
    iline += 1
    label = ""
    icol = 0
    for sim in sims:
        col = cols[icol]
        traga,rraga,rfrac,nfrac,vrfrac,vtfrac = solns[sim.Name()]
        # Don't do rfrac, makes plot too busy, kinda already done in other plot
        #plt.plot(traga,rfrac,"-"+col,label=sim.Label()+ "r/r$_s$")
        plt.plot(traga,nfrac,col+line[iline],label=label+" n/n$_0$")
        icol += 1
    iline += 1
    '''
    label = ""
    icol = 0
    for sim in sims:
        col = cols[icol]
        traga,rraga,rfrac,nfrac,vrfrac,vtfrac, denss,vrs = solns[sim.Name()]
        # Remove negative values
        vrfrac = -vrfrac # Only show infall
        #vrfrac[vrfrac <= 0] = vrfrac[vrfrac > 0].min()
        plt.plot(traga,vrfrac,color=col,
                 linestyle=line[iline],label=label+" v$_{r}$/c$_i$")
        icol += 1
    iline += 1
    label = ""
    icol = 0
    for sim in sims:
        col = cols[icol]
        traga,rraga,rfrac,nfrac,vrfrac,vtfrac,denss,vrs = solns[sim.Name()]
        # Factor in Raga/Tremblin equation
        Ffrac = rfrac**(-3.0/4.0) * nfrac ** (-0.5)
        plt.plot(traga,Ffrac,color=col,
                 linestyle=line[iline],label=label+"F(r,t)")
        icol += 1
    iline += 1
    '''
    label = ""
    icol = 0
    for sim in sims:
        col = cols[icol]
        traga,rraga,rfrac,nfrac,vrfrac,vtfrac = solns[sim.Name()]
        vtfrac = np.sqrt(vtfrac)
        plt.plot(traga,vtfrac,color=col,
                 linestyle=line[iline],label="v$_{t}$/c$_{i}$")
        icol += 1
    '''
    iline += 1
    plt.plot([0,4],[1,1],"k--")
    plt.xlim([0,4])
    plt.ylim([1e-3,1e1])
    plt.xlabel("Time / Myr")
    plt.ylabel("Fractional Value")
    plt.yscale("log")
    leg = plt.legend(fontsize=11,loc="lower right",
                     ncol=iline,frameon=False)
    txts = leg.get_texts()
    nsim = len(sims)
    for itxt in range(0,nsim):
        left = txts[itxt]
        right = txts[itxt+nsim]
        right.set_y(left.get_position()[1])
        left.set_verticalalignment("baseline")
        right.set_verticalalignment("baseline")
    plt.savefig("../plots/raga/"+name+"/ragafractionals_"+proftype+".pdf")

    # Now plot the density and velocity at the interface
    plt.clf()
    solns = solnsspitz
    iline = 0
    icol = 0
    line = ["-","--",":","-."]
    label = sim.Label()
    # Velocity
    label = ""
    icol = 0
    for sim in sims:
        col = cols[icol]
        traga,rraga,rfrac,nfrac,vrfrac,vtfrac,denss,vrs = solns[sim.Name()]
        # Remove negative values
        vrfrac = -vrfrac # Only show infall
        #vrfrac[vrfrac <= 0] = vrfrac[vrfrac > 0].min()
        plt.plot(traga,vrfrac,color=col,
                 linestyle=line[iline],label=label+" v$_{r}$/c$_i$")
        icol += 1
    iline += 1
    # Density
    label = ""
    icol = 0
    for sim in sims:
        col = cols[icol]
        traga,rraga,rfrac,nfrac,vrfrac,vtfrac,denss,vrs = solns[sim.Name()]
        # Remove negative values
        denstoplot = np.log10(denss)
        #vrfrac[vrfrac <= 0] = vrfrac[vrfrac > 0].min()
        plt.plot(traga,denstoplot,color=col,
                 linestyle=line[iline],label=label+" $n_{ext}$")
        icol += 1
    iline += 1
    # F(r,t)
    label = ""
    icol = 0
    for sim in sims:
        col = cols[icol]
        traga,rraga,rfrac,nfrac,vrfrac,vtfrac,denss,vrs = solns[sim.Name()]
        # Factor in Raga/Tremblin equation
        Ffrac = rfrac**(-3.0/4.0) * nfrac ** (-0.5)
        plt.plot(traga,Ffrac,color=col,
                 linestyle=line[iline],label=label+"F(r,t)")
        icol += 1
    iline += 1
    # Tidy up plot
    plt.plot([0,1],[1,1],"k--")
    plt.xlim([0,1])
    plt.ylim([1e-2,1e1])
    plt.xlabel("Time / Myr")
    plt.ylabel("log($n_{H}$), $v_{r} / c_{i}$, $F(r,t)$")
    plt.yscale("log")
    leg = plt.legend(fontsize=11,loc="lower right",
                     ncol=iline,frameon=False)
    txts = leg.get_texts()
    nsim = len(sims)
    for itxt in range(0,nsim):
        left = txts[itxt]
        right = txts[itxt+nsim]
        right.set_y(left.get_position()[1])
        left.set_verticalalignment("baseline")
        right.set_verticalalignment("baseline")
    plt.savefig("../plots/raga/"+name+"/ragadensvel_"+proftype+".pdf")
    

if __name__=="__main__":
    '''
    sims = [Hamu.Simulation(s) for s in ["N47_M4_B02",
                                         "N48_M4_B02",
                                         "N49_M4_B02"]]
    PlotForSims(sims,"photons","ana")
    PlotForSims(sims,"photons","sim")
    sims = [Hamu.Simulation(s) for s in ["N48_M4_B02",
                                         "N48_M4_B02_C2",
                                         "N48_M4_B02_C"]]
    PlotForSims(sims,"compact","ana")
    PlotForSims(sims,"compact","sim")
    '''
    sims = [Hamu.Simulation(s) for s in ["noturb48",
                                         "noturb49",
                                         "noturb50"]]
    PlotForSims(sims,"freefall","ana")
    PlotForSims(sims,"freefall","sim")

    '''
    sims = [Hamu.Simulation(s) for s in ["N48_M4_B02",
                                         "N48_M4_B02_F2",
                                         "N48_M4_B02_F3"]]
    PlotForSims(sims,"tstart","ana")
    PlotForSims(sims,"tstart","sim")
    '''
