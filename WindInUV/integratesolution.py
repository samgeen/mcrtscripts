'''
Integrate the wind + UV photoionisation + UV radiation pressure solution through time
Sam Geen, October 2018
'''


import sys

import time
import tqdm

import customplot

import matplotlib.pyplot as plt

import numpy as np
import scipy.interpolate
import scipy.signal

from consts import *

import stars
from stars import _Star # needed to make caching work

import solvecriteria
from consts import *
from collections import OrderedDict

invsix = 1.0/6.0
def integrator(y,func,t,dt):
    # Runge-Kutta 4 integrator
    # y = system state
    # func = rate of change of y
    # t = time (useful for stuff that uses a zero point in time like the stellar evolution model)
    # dt = step size
    k1 = dt * func(t,y) 
    k2 = dt * func(t+0.5*dt,y+0.5*k1)
    k3 = dt * func(t+0.5*dt,y+0.5*k2)
    k4 = dt * func(t+dt,y+k3)
    # new t is t+dt
    return y + invsix * (k1 + 2.0*k2 + 2.0*k3 + k4)

class WindSolver(object):
    def __init__(self,star,cloud,cool=True,nowind=False):
        self._star = star
        self._cloud = cloud
        self._cool = cool
        self._nowind = nowind

    def __call__(self,t,ri):
        return solvecriteria.drdt_wind(t,ri,self._star,self._cloud,self._cool,self._nowind)

def runforwind(star,cloud,nowind,rini=0.01,dtyrs = 100.0):
    solver = WindSolver(star,cloud,cool=True,nowind=nowind)
    # Timestep every N years for 3 Myr
    dt = dtyrs * yrins
    tend = 1e6 * yrins
    times = np.arange(0.0,tend+0.1*dt,dt)
    ris = times*0.0 + rini*pc
    rws = times*0.0
    nis = times*0.0
    rpfacts = times*0.0
    momentums = times*0.0
    ri = ris[0]
    t = times[0]+dt
    print "Running!"
    rw = 0.0
    rpfact = 0.0
    for i in tqdm.tqdm(range(len(times)-1)):
        try:
            # Find ri and check for good values
            riold = ri
            ri = integrator(ri,solver,t,dt)
            if np.isnan(ri):
                raise solvecriteria.IntegratorError
            # Find ni
            ni = solvecriteria.Findni(t,star,cloud,ri,nowind)
            # Find wind radius
            if not nowind:
                rw = ri / solvecriteria.RivsRw(t,star,cloud,ri,dustfrac=0.73)
            else:
                rw = 0.0
            # Find the relative pressures
            rpfact = solvecriteria.deltaPrpoverP(t,star,cloud,ri,nowind=nowind)
            # Find shell momentum
            drdt = (ri - riold) / dt
            momentum = cloud.MassInsideR(ri) * drdt
        except solvecriteria.IntegratorError:
            ris = ris[0:i]
            rws = rws[0:i]
            nis = nis[0:i]
            rpfacts = rpfacts[0:i]  
            momentums = momentums[0:i]  
            times = times[0:i]
            break
        ris[i+1] = ri
        rws[i+1] = rw
        nis[i+1] = ni
        rpfacts[i+1] = rpfact
        momentums[i+1] = momentum
        t += dt
    print "Done!"
    # Calculate wind radius
    # Return scaled units
    tout = times/yrins/1e6
    #tout = scipy.signal.resample(tout,1000)
    riout = ris/pc
    rwout = rws/pc
    niout = nis
    rpfactout = rpfacts
    momentumout = momentums
    return tout, riout, rwout, niout, rpfactout, momentumout

def run(metal=0.014, cool=True, mcloud = 1e3, surfdens = 20.0, sfe=None, rini=0.001, seed=1):
    '''
    SFE = fraction of cloud mass (if None, use a single star)
    '''
    starmass = 30.0
    # Surface density
    surf = surfdens*Msun/pc**2
    # Stellar mass
    #starmasses = stars.starmasses
    #nstarmass = len(starmasses)
    # Cloud mass
    cloudmass = mcloud*Msun
    #surfs = np.logspace(1e1,1e4,nsurf)*Msun/pc**2
    # Other setup stuff
    #t = 1e6*yrins
    cloud = solvecriteria.Cloud(cloudmass,surf,metal,accreting=False)
    if sfe is None:
        star = stars.Star(starmass,metal)
        sfetxt = ""
    else:
        starmass = sfe*mcloud
        star = stars.ClusterOnTheFly(starmass,metal,seed=seed)
        sfetxt = "_cluster"+str(seed)
    #print "Stall radius in pc:", solvecriteria.rstall(t,star,cloud)/pc
    # Run for the wind
    tw,riw,rw,niw,rpfactw,momentumw  = runforwind(star,cloud,nowind=False,rini=rini,dtyrs = 100.0)
    tnw,rinw,dum,ninw,rpfactnw,momentumnw = runforwind(star,cloud,nowind=True,rini=rini,dtyrs = 100.0)
    #import pdb; pdb.set_trace()
    # Radius plot
    plt.clf()
    plt.plot(tw,riw,"k",label="$r_i$, wind")
    plt.plot(tw,rw,"k:",label="$r_w$, wind")
    plt.plot(tnw,rinw,"k",alpha=0.5,label="$r_i$, no wind")
    plt.legend(frameon=False,loc="upper left",fontsize="large")
    plt.xlabel("Time / Myr")
    plt.ylabel("Radius / pc")
    plt.xscale("log")
    plt.yscale("log")
    plt.savefig("plots/windradius"+sfetxt+"_M"+str(int(starmass))+".pdf")
    # Density plot
    plt.clf()
    plt.plot(tw,niw,"k",label="Wind")
    plt.plot(tnw,ninw,"k",alpha=0.5,label="No wind")
    plt.legend(frameon=False,loc="upper right",fontsize="large")
    plt.xlabel("Time / Myr")
    plt.ylabel("HII region density $n_{i}$ / cm$^{-3}$")
    plt.xscale("log")
    plt.yscale("log")
    ax = plt.gca()
    xlims = ax.get_xlim()
    ylims = ax.get_ylim()
    xpos = xlims[0]# + 0.1*(xlims[1]-xlims[0])
    ypos = ylims[0]# + 0.1*(ylims[1]-ylims[0])
    startxt = str(int(starmass))
    plt.text(0.01,0.01,startxt+" M$_{\odot}$ star in a 10$^4$ M$_{\odot}$ cloud with $\Sigma_0$=100 M$_{\odot}$/pc$^2$", 
        ha='left', va='bottom',transform=ax.transAxes)
    plt.savefig("plots/winddensity"+sfetxt+"_M"+str(int(starmass))+".pdf")
    # Radiation pressure factor
    plt.clf()
    plt.plot(tw,rpfactw,"k",label="Wind")
    plt.plot(tnw,rpfactnw,"k",alpha=0.5,label="No wind")
    plt.legend(frameon=False,loc="lower right",fontsize="large")
    plt.xlabel("Time / Myr")
    plt.ylabel("$C_{rp} \equiv \Delta P_{rp}/P_0$")
    plt.xscale("log")
    plt.yscale("linear")
    plt.ylim([0.09,0.21])
    plt.savefig("plots/dPrpoverP"+sfetxt+"_M"+str(int(starmass))+".pdf")
    # Momentum
    plt.clf()
    plt.plot(tw,momentumw,"k",label="Wind")
    plt.plot(tnw,momentumnw,"k",alpha=0.5,label="No wind")
    plt.legend(frameon=False,loc="lower right",fontsize="large")
    plt.xlabel("Time / Myr")
    plt.ylabel("Momentum / g cm/s")
    plt.xscale("log")
    plt.yscale("log")
    #plt.ylim([0.09,0.21])
    plt.savefig("plots/momentum"+sfetxt+"_M"+str(int(starmass))+".pdf")


if __name__=="__main__":
    run()
    # for i in range(1,11):
    #     run(metal=0.014, cool=True, mcloud = 1e5, surfdens = 100.0, sfe=0.1, rini=0.1,
    #         seed=i)