'''
Find the photon densities around the wind bubble
Sam Geen, March 2015
'''

from startup import *

import os, errno
from pymses.filters import CellsToPoints
from pymses.utils import constants as C

import matplotlib.patheffects as pe 

#from sklearn.cluster import DBSCAN
#import skimage.measure
#import fastcluster

import starrelations
import stellars

LSIZE = 512

class DensityFinder(object):
    # Find the densities around the wind bubble
    def __init__(self, snap):
        self._snap = snap

    def Run(self):
        # Return empty if no stars made yet to speed this up
        stellar = stellars.FindStellar(self._snap)
        if len(stellar.mass) == 0.0:
            return np.array([])
        # Find densities around the wind bubble
        dens = self._FindDensities()
        return dens
        
    def _FindDensities(self):
        # Make regular grid to find regions inside
        amr = self._snap.amr_source(["rho","xHII","P","xHeII","xHeIII","vel"])
        lmin = self._snap.info["levelmin"]
        lsize = LSIZE
        boxlen = self._snap.info["boxlen"]
        pos = [0.5,0.5,0.5]
        radius = 0.49 # Remove whatever weird edge stuff happens
        coords = np.linspace(-0.5,0.5,lsize)*2.0*radius
        grid = np.meshgrid(coords+pos[0],coords+pos[1],coords+pos[2])
        points = np.array([grid[0].flatten(),grid[1].flatten(),grid[2].flatten()])
        points = points.T
        samples = pymses.analysis.sample_points(amr,points)
        # Find cubes
        # Density cube
        scale = hydrofuncs.scale_by_units(self._snap,"rho")
        denscube = scale(samples)
        # Wind location cube
        # Here we look for two criteria, hot or fast
        # Mash the results into a cube that looks like an ionisaton state
        scale = hydrofuncs.scale_by_units(self._snap,"T")
        Tcube = scale(samples)
        scale = hydrofuncs.scale_by_units(self._snap,"spd")
        Vcube = scale(samples) 
        windcube = Tcube * 0.0 + 0.0
        windcube[Tcube > 1e5] = 1.0 # 1e5 K, find hot wind shocked gas
        windcube[Vcube > 5e2] = 1.0 # 500 km/s, find fast winds
        denscube = denscube.reshape([lsize,lsize,lsize])
        windcube = windcube.reshape([lsize,lsize,lsize])
        # Remove edge cells from mask to prevent index over/underflows
        mask = windcube * 0.0
        windmask[1:lsize-1,1:lsize-1,1:lsize-1] = windcube[1:lsize-1,1:lsize-1,1:lsize-1]
        # Indices containing winds in x, y, z
        wix, wiy, wiz = np.where(windmask == 1)
        # March values around cells neighbouring wind cells
        tosample = windcube * 0.0
        tosample[wix+1,wiy,wiz] = 1.0
        tosample[wix-1,wiy,wiz] = 1.0
        tosample[wix,wiy+1,wiz] = 1.0
        tosample[wix,wiy-1,wiz] = 1.0
        tosample[wix,wiy,wiz+1] = 1.0
        tosample[wix,wiy,wiz-1] = 1.0
        # Don't sample the wind cells
        tosample[windcube] = 0.0
        edgedens = denscube[tosample].flatten()
        # Return the densities around the wind cells
        return edgedens

def _FindHIIRegionDensities(snap,lsize=256):
    # Find the densities immediately around the wind bubbles
    global LSIZE
    LSIZE = lsize
    finder = DensityFinder(snap)
    return finder.Run()
FindHIIRegionDensities = Hamu.Algorithm(_FindHIIRegionDensities)

def plot(sims,labels):
    # TODO: CHANGE, THIS USES THE FIND HII REGIONS CODE
    plt.clf()
    numcols = len(sims)
    fig, axes = plt.subplots(1,numcols,sharex=False,sharey=True)
    first = True
    for ax, sim, label in zip(axes, sims, labels):
        lsize = 512
        t, x = timefuncs.timefunc(sim,FindHIIRegions,False,False,1,"xHII",lsize)
        t, w = timefuncs.timefunc(sim,FindHIIRegions,False,False,1,"wind",lsize)
        tcreated, sfe = starrelations.runforsim(sim.Name(),"firsttime")
        t -= tcreated 
        ax.set_xlabel("Time after 1st star formed / Myr")
        if first:
            ax.set_ylabel("Number of Regions")
        color = linestyles.Colour(sim.Name())
        ax.plot(t,x,color=color,linestyle="-",label="Photoionised regions",
                path_effects=[pe.Stroke(linewidth=5, foreground='k'), pe.Normal()]) 
        ax.plot(t,w,color=color,linestyle="--",label="Wind Bubbles",
                path_effects=[pe.Stroke(linewidth=5, foreground='k'), pe.Normal()])
        textloc = 0.1
        ax.text(0.95, textloc,label, ha='right', va="top", transform=ax.transAxes)
        if first:
            ax.legend(fontsize="x-small",loc="upper left",frameon=False)
        first = False
        ax.set_xlim([-0.5,5]) 
    fig.subplots_adjust(wspace=0)
    fig.set_size_inches(7*numcols,6)
    fig.savefig(plotfolder+"findHIIregions.pdf", dpi=80)
        
if __name__=="__main__":
    #sims = hamusims # [hamusims[x] for x in ["IMF1_04","IMF2_04",""]]
    #labels = ["IMF 1","IMF 2","Massive Cloud"]
    #plot(sims, labels)
    for name, sim in hamusims.iteritems():
        if "WIND" in name:
            for snap in sim.Snapshots():
                densities = FindHIIRegionDensities(snap)