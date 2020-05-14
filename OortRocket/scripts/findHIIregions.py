'''
Find HII regions in the simulation
Sam Geen, March 2015
'''

from startup import *

import os, errno
from pymses.filters import CellsToPoints
from pymses.utils import constants as C

import matplotlib.patheffects as pe 

#from sklearn.cluster import DBSCAN
import skimage.measure
#import fastcluster

import starrelations

LSIZE = 512

class HIIRegionFinder(object):
    def __init__(self, snap):
        self._snap = snap

    def Run(self,findertype="xHII"):
        # Make cube of xHII
        cube = self._MakeCube(findertype)
        # Run through cube
        labelcube = self._FindHIIRegions(cube)
        nregions = labelcube.max()
        return nregions
        
    def _MakeCube(self,findertype):
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
        # Sample xHII
        if findertype == "xHII":
            # Make a cube of ionisation state
            scale = hydrofuncs.scale_by_units(self._snap,"xHII")
            hydrocube = scale(samples)
        if findertype == "wind":
            # Winds
            # Here we look for two criteria, hot or fast
            # Mash the results into a cube that looks like an ionisaton state
            scale = hydrofuncs.scale_by_units(self._snap,"T")
            Tcube = scale(samples)
            scale = hydrofuncs.scale_by_units(self._snap,"spd")
            Vcube = scale(samples) 
            hydrocube = Tcube * 0.0 + 0.0
            hydrocube[Tcube > 1e5] = 1.0 # 1e5 K, find hot wind shocked gas
            hydrocube[Vcube > 5e2] = 1.0 # 500 km/s, find fast winds
        hydrocube = hydrocube.reshape([lsize,lsize,lsize])
        # Return The Hydrocube
        return hydrocube

    def _FindHIIRegions(self, xcube):
        '''
        Finds volumes of contiguous HII in the cube
        xcube - cube of xHII values
        '''
        # Numpy way of making an empty cube, don't @ me
        nx, ny, nz = xcube.shape
        intcube = np.zeros((nx,ny,nz),dtype=np.int)
        intcube[xcube > 0.1] = 1
        labelcube = skimage.measure.label(intcube)
        return labelcube



def _FindHIIRegions(snap,findertype="xHII",lsize=256):
    # Find HII regions
    # findertype - xHII or wind
    global LSIZE
    LSIZE = lsize
    finder = HIIRegionFinder(snap)
    return finder.Run(findertype)
FindHIIRegions = Hamu.Algorithm(_FindHIIRegions)

def plot(sims,labels):
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
    sims = [hamusims[x] for x in ["IMF1_04","IMF2_04","MASS_04"]]
    labels = ["IMF 1","IMF 2","Massive Cloud"]
    plot(sims, labels)
    for sim in [hamusims[x] for x in ["IMF1_04","IMF2_04","MASS_04"]]:
        for snap in sim.Snapshots():
            regions = FindHIIRegions(snap,"wind")
            regions = FindHIIRegions(snap,"xHII")
