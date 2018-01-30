'''
Initial profiles
Sam Geen, July 2015
'''

import customplot
import matplotlib.pyplot as plt
import matplotlib.lines as mlines

import Hamu
import profilesphere
import linestyles

import outflowmodel, rdotplot, testprofilefit
from pymses.utils import constants as C

def profile(snap):
    r,p = profilesphere.profileHamu(snap,"rho",1e6)
    p[0] = p[1] # fix zero sampling problem
    p[p <= 0] = p.max()
    return r,p

def profilefit(snap):
    boxrad = snap.RawData().info["boxlen"]/2.0
    r,p = profilesphere.profileHamu(snap,"rho") 
    return testprofilefit.FindFit(r,p,
                                  rcut=boxrad*testprofilefit.rextent,
                                  params=False)
    
def run(extras=False):
    simnames = ["N48_M4_B02",
                "N48_M4_B02_C2",
                "N48_M4_B02_C"]
    labels = ["Fiducial",
              "More Compact",
              "Most Compact"]
    # Run for each sim
    for simname, label in zip(simnames, labels):
        sim = Hamu.Simulation(simname)
        snap = sim.Snapshots()[0]
        col = linestyles.col(simname)
        r, p = profile(snap)
        plt.plot(r,p,label=label,color=col)
        # Plot the profile at t=t_ff, power law fit
        if extras:
            # Get snap
            tff = rdotplot.Findtstart(simname)
            tcode = tff / snap.RawData().info["unit_time"].express(C.Myr)
            snaptff = sim.FindAtTime(tcode)
            # Sampled profile
            r, p = profile(snaptff)
            plt.plot(r,p,color=col,linestyle="--")
            # Power law fit
            r, p = profilefit(snaptff)
            plt.plot(r,p,color=col,linestyle=":")
    # Plot setup
    plt.xlabel("Radius / pc")
    plt.ylabel("Density / atoms cm$^{-3}$")
    plt.yscale("log")
    plt.xlim([0,8])
    leg1 = plt.legend(fontsize="small",frameon=False)
    if extras:
        line1 = mlines.Line2D([], [], color='k', label='$t=0$')
        line2 = mlines.Line2D([], [], color='k', linestyle=":",
                              label='$t=t_{ff}$')
        line3 = mlines.Line2D([], [], color='k', linestyle=":",
                              label='Power Law Fit')

        leg2 = plt.legend(handles=[line1,line2,line3],
                          ncol=1,fontsize="small",
                          frameon=False,loc="center right")
        plt.gca().add_artist(leg1)

    plt.savefig("../plots/compact/initdensities.pdf")

if __name__=="__main__":
    run(extras=False)
