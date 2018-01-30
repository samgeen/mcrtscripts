'''
Solve the rocket equations as per Bertoldi + McKee 1990
Sam Geen, May 2015
'''

import customplot
import matplotlib.pyplot as plt

import numpy as np

import Hamu

# My code badly needs refactoring
import outflowmodel, rdotplot

msolaring = 1.9891e33
G = 6.674e-8

def rcloudforsim(simname,start=None):
    if start is None:
        rcloud = 3.0*outflowmodel.pcincm
    else:
        rcloud = start
    if "_C2" in simname:
        rcloud *= 0.75**2
    elif "_C" in simname:
        rcloud *= 0.5**2
    return rcloud

def runforsim_OLD(sim):
    simname = sim.Name()
    print "Running rocketclump for", simname
    # Cloud properties
    tff = rdotplot.Findtstart(simname)
    outflowmodel.tstart = tff
    m1 = 6e2 * msolaring # Pick the most massive clump
    mcloud = 1e4 * msolaring
    rcloud = rcloudforsim(simname) * 0.3
    vesc = np.sqrt(2.0 * G * mcloud / rcloud)
    vR = 0.8 * outflowmodel.cs
    tev = 30.0*outflowmodel.Myrins # roughly
    betaR = tev * vR / rcloud
    # Solvers
    def dmdt(m,R):
        return -2.5 * m**0.6 * R ** (-0.4) / tev
    def FindR(m):
        mpow = m ** 0.4
        R = 1.5 * betaR * ((0.4 * v1 / vR + 1) * (1 - mpow) + \
                               mpow * np.log(mpow)) + 1.0
        return R ** (5.0/3.0)
    # Set initial values and go
    ntime = 10000
    dt = 5.0*outflowmodel.Myrins / ntime
    ms = np.zeros((ntime))
    Rs = np.zeros((ntime))
    ts = np.zeros((ntime))
    ms[0] = 1.0 # this is the scaled value!
    Rs[0] = 1.0
    v1 = vesc
    print v1 / vR
    for itime in range(0,ntime-1):
        #print ms[itime], Rs[itime]
        ms[itime+1] = ms[itime]+dmdt(ms[itime],Rs[itime])*dt
        Rs[itime+1] = FindR(ms[itime])
        ts[itime+1] = ts[itime]+dt
        #import pdb; pdb.set_trace()
    ts /= outflowmodel.Myrins
    return ts, Rs, ms

def runforsim(sim):
    simname = sim.Name()
    print "Running rocketclump for", simname
    # Cloud properties
    tff = rdotplot.Findtstart(simname)
    outflowmodel.tstart = tff
    m1 = 6e2 * msolaring # Pick the most massive clump
    mcloud = 1e4 * msolaring
    radius = 0.25 * outflowmodel.pcincm
    rcloud = rcloudforsim(simname)
    R1 = rcloudforsim(simname,2.0)
    ncloud = 1e3 # I dunno
    vesc = np.sqrt(G * mcloud / rcloud)
    vR = 0.8 * outflowmodel.cs
    tev = 10*outflowmodel.Myrins # roughly
    w = 0.0#-outflowmodel.FindProfilePower(sim)
    #print "power law index w", w
    # Correct for photons
    if "N47" in simname:
        tev *= 10.0**0.2
    if "N49" in simname:
        tev /= 10.0**0.2
    betaR = tev * vR / rcloud
    mp = outflowmodel.mH / outflowmodel.X
    # Solvers
    def Menclosed(m,R):
        # Mass enclosed inside R
        return mcloud * (R/rcloud)**(3.0-w)
    def dmdt(m,R):
        return -2.5 * (m/m1)**0.6 * (R/R1) ** (-0.4) / tev * m1
    def dvdt(m,R):
        g = -G * Menclosed(m,R) / R**2
        ptherm = outflowmodel.cs**2 * ncloud * mp
        # Force = ptherm * surface
        ftherm = ptherm * FindRadius(m,R)**2.0
        atherm = ftherm / m
        # NOTE: dmtilde/dt / mttilde = dm/dt / m
        #print -vR * dmdt(m,R) / m , g, atherm
        return -vR * dmdt(m,R) / m + g# + atherm
    def FindRadius(m,R):
        # r propto m ^ (2/5)
        return radius * (m/m1)**0.4
    # Set initial values and go
    ntime = 20000
    dt = 5.0*outflowmodel.Myrins / ntime
    ms = np.zeros((ntime))
    Rs = np.zeros((ntime))
    ts = np.zeros((ntime))
    ms[0] = m1
    R1 = rcloud
    Rs[0] = R1
    v1 = 0.0#vesc
    v = v1
    for itime in range(0,ntime-1):
        #print ms[itime], Rs[itime]
        v += dvdt(ms[itime],Rs[itime])*dt
        ms[itime+1] = ms[itime] + dmdt(ms[itime],Rs[itime])*dt
        Rs[itime+1] = Rs[itime] + v*dt
        ts[itime+1] = ts[itime]+dt
        if Rs[itime+1] < 0.0:
            break
        #import pdb; pdb.set_trace()
    ms /= m1
    Rs /= R1
    ts /= outflowmodel.Myrins
    return ts, Rs, ms

def plotforsims(simnames):
    col = ["k","r","m","b","c"]
    icol = 0
    for simname in simnames:
        sim = Hamu.Simulation(simname)
        ts, Rs, ms = runforsim(sim)
        plt.plot(ts, Rs,col[icol]+"-",label=sim.Label())
        plt.plot(ts, ms,col[icol]+"--")
        icol += 1
    plt.legend(fontsize="xx-small")
    plt.yscale("log")
    plt.ylim([1e-2,10])
    plt.xlabel("Time / Myr")
    plt.ylabel("Fractional R (m = dashed)")
    plt.savefig("../plots/clumps/rocketclump.pdf")
        

if __name__=="__main__":
    simnames = ["N48_M4_B02",
                "N49_M4_B02",
                "N47_M4_B02",
                "N48_M4_B02_C2",
                "N48_M4_B02_C"]
    plotforsims(simnames)
