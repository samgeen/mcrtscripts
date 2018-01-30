'''
Calculate statistics of the radius of the ionisation front
Sam Geen, November 2014
'''

import customplot
import pymses, Hamu
import numpy as np
import matzner2002, errtimeplot
from scipy.spatial import ConvexHull
import collections
import matplotlib.pyplot as plt
from pymses.utils import constants as C
from pymses.filters import CellsToPoints
import profilesphere,rayprof

sims = collections.OrderedDict()

def Samples(amr,npts=1e6):
    centre=[0.5,0.5,0.5]
    radius = 0.5
    sph = pymses.utils.regions.Sphere(centre, radius)
    return pymses.analysis.sample_points(amr, sph.random_points(npts))

def GetRanges(snap):
    amr = snap.amr_source(["rho","xHII"])
    cell_source = CellsToPoints(amr)
    cells = Samples(amr)
    rhos = cells["rho"]
    ions = cells["xHII"]
    posns = cells.points-0.5
    thresh = 0.2 # The code considers everything above 0.2 to be ionised
    try:
        pions = posns[ions > thresh,:]
    except:
        return np.zeros(5)
    # Get the points on a convex hull around the ionised gas
    if len(pions) == 0:
        return np.zeros(5)
    try:
        hull = ConvexHull(pions)
    except:
        return np.zeros(5)
    radii = np.sqrt(np.sum(hull.points[hull.vertices,:]**2,1))
    # Find ranges
    return np.percentile(radii, [1,25,50,75,100])

GetRangesHamu = Hamu.Algorithm(GetRanges)

def MedianRadiusCoarse(snap):
    # Uses the function above to get a coarse sampling of the median radius
    # median is the 50th percentile above
    return GetRangesHamu(snap.hamusnap)[2] 

def MedianRadiusProfile(snap):
    # Get the mean radius from the profiles
    # Assuming all gas is either x=0 or x=1, find point where x = 0.5
    # This gives us the mean radius of the ionisation front by angular position
    # rcut included for when the centre needs to be sampled better
    # HACK
    rcut = None
    if "_C/" in snap.hamusnap.CachePath(): 
        rcut = 0.05
    if rcut is None:
        r,x = profilesphere.profileHamu(snap.hamusnap,"xHII",1e6)
    else:
        r,x = profilesphere.profileHamu(snap.hamusnap,"xHII",1e6,rcut=rcut)
    xlim = 0.5
    if x.max() < xlim:
        # No ions? return zero
        return 0.0
    #if x.max() < xlim:
    #    # Allow low levels of ionisation to count as a radius in a pinch
    #    xlim = x.max()
    return r[x >= xlim].max()
            
def plot():
    global sims
    Myr = 3.15569e13
    plt.clf()
    for simname,col in sims.iteritems():
        sim = Hamu.Simulation(simname)
        times = sim.Times()
        utime = sim.Snapshots()[0].RawData().info["unit_time"].express(C.Myr)
        #boxlen = sim.Snapshots()[0].RawData().info["boxlen"]
        boxlen = 1.0
        times = times * utime - 1.25
        ntimes = len(times)
        rmin = np.zeros((ntimes))
        r25p = np.zeros((ntimes))
        rmed = np.zeros((ntimes))
        r75p = np.zeros((ntimes))
        rmax = np.zeros((ntimes))
        i = 0
        # Fill data
        for snap in sim.Snapshots():
            #rmin[i],r25p[i],rmed[i],r75p[i],rmax[i] = GetRangesHamu(snap)
            rmin[i],r25p[i],rmed[i],r75p[i],rmax[i] = \
                errtimeplot.radialstatsHamu(snap)
            i += 1
        # Make lines
        taba = np.concatenate((times,times[::-1])) # there and back agai 
        minmax = np.concatenate((rmin,rmax[::-1]))
        iqr = np.concatenate((r25p,r75p[::-1]))
        plt.fill(taba,minmax*boxlen,alpha=0.33,edgecolor='none',facecolor=col)
        plt.fill(taba,iqr*boxlen,   alpha=0.33,edgecolor='none',facecolor=col)
        plt.plot(times,rmed*boxlen,col,label=simname)
        # Compare with matzner
        tm2002,rm2002 = matzner2002.Findrii(sim)
        plt.plot(tm2002,rm2002,col+"--")
    plt.xlabel("Time / Myr")
    plt.ylabel("Ionisation Front Radius / pc")
    plt.legend(loc="upper left",fontsize="x-small")
    plt.savefig("../plots/radialstats.pdf")

if __name__=="__main__":
    sims["N00_M4_B02"] = "k"
    sims["N47_M4_B02"] = "b"
    sims["N48_M4_B02"] = "r"
    sims["N49_M4_B02"] = "g"
    plot()
