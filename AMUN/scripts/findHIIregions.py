'''
Find HII regions in the simulation
Sam Geen, March 2015
'''

from startup import *

import os, errno
from pymses.filters import CellsToPoints
from pymses.utils import constants as C

#from sklearn.cluster import DBSCAN
import fastcluster

class HIIRegionFinder(object):
    def __init__(self, snap, xthresh=0.1):
        self._snap = snap
        self._xthresh = xthresh

    def Run(self):
        amr = self._snap.amr_source(["rho","xHII"])
        cell_source = CellsToPoints(amr)
        cells = cell_source.flatten()
        regioncells = self._FilterByxHII(cells)
        regions = self._FindHIIRegions(regioncells)
        return regions
        
    def _FilterByxHII(self, cells):
        xHII = cells["xHII"]
        return cells.filtered_by_mask(xHII > self._xthresh)

    def _FindHIIRegions(self, cells):
        '''
        Finds the nearest neighbours of the clumps
        '''
        # sklearn clustering algorithm
        # "single uses the minimum of the distances between all observations of the two sets."
        #algo = DBSCAN()
        #clustering = algo.fit(cells.points)
        #clumps = np.unique(clustering.labels_)
        clustering = fastcluster.linkage_vector(cells.points)
        import pdb; pdb.set_trace()
        return clumps

def FindHIIRegions(snap):
    finder = HIIRegionFinder(snap)
    return finder.Run()

class ClumpMaker(object):
    '''
    Identifies a single clump and returns its stats
    '''
    def __init__(self, snap, cells):
        self._snap = snap
        self._ncells = len(cells["rho"])
        self._cells = cells
        self._inclump = np.zeros((self._ncells),dtype="int32")
        self._radii = cells.get_sizes()*2.0 # Make the radii "big enough"
        self._setup = False
        self._mass = 0.0
        self._peakdens = 0.0
        self._centre = 0.0
        self._radius = 0.0
        self._distance = 0.0
        self._bfactor = 0.0

    def _Setup(self):
        if not self._setup:
            self._FindClump()
            self._FindClumpProperties()
            self._setup = True
            

    def _FindClump(self):
        newfound = np.array([0]) # Pick the first cell
        self._inclump[0] = 1
        # Loop while newly found cells exist
        icount = 0
        while len(newfound) > 0:
            icount += 1
            # Make new newfound list and search through previous list
            oldfound = newfound
            newfound = np.array([],dtype="int32")
            # Find new neighbours
            for icell in oldfound:
                found = self._FindNeighbours(icell)
                newfound = np.concatenate((newfound,found))
            # Uniquify newfound
            newfound = np.unique(newfound)
            # Identify these as being in the clump
            self._inclump[newfound] = 1
        print "Num cells in clump:", len(self._inclump[self._inclump == 1])

    def _FindClumpProperties(self):
        cells = self._cells.filtered_by_mask(self._inclump == 1)
        munit = self._snap.info["unit_mass"].express(C.Msun)
        dunit = self._snap.info["unit_density"].express(C.H_cc)
        runit = self._snap.info["unit_length"].express(C.pc)
        magunit = self._snap.info["unit_mag"].express(C.Gauss)
        rhos = cells["rho"]*dunit
        self._mass = np.sum(cells["rho"]*cells.get_sizes()**3.0*munit)
        self._peakdens = rhos.max()
        central = rhos == self._peakdens
        cells.points -= 0.5 # centre of the volume
        cells.points *= runit
        self._centre = cells.points[central,:]+0.0
        self._distance = np.sqrt(np.sum(self._centre**2))
        print "CENTRE", self._centre
        points = cells.points - self._centre
        self._radius = np.sqrt(np.sum(points**2,1).max())
        # 1e6 because Gauss -> microgauss
        #magnetic = 0.5*(cells["B-left"]+cells["B-right"])*magunit*1e6
        # b from Bertoldi & McKee 1990
        #self._bfactor = (np.sqrt(np.sum(magnetic**2,1))*rhos**(-0.5))[central]
        
    @property
    def inclump(self):
        self._Setup()
        return self._inclump

    @property
    def centre(self):
        self._Setup()
        return self._centre

    @property
    def distance(self):
        self._Setup()
        return self._distance

    @property
    def peakdensity(self):
        self._Setup()
        return self._peakdens

    @property
    def radius(self):
        self._Setup()
        return self._radius

    @property
    def mass(self):
        self._Setup()
        return self._mass

    #@property
    #def bfactor(self):
    #    self._Setup()
    #    return self._bfactor
            
    def _FindNeighbours(self,icell):
        '''
        Find neighbouring cells
        TODO: TIDY UP INTERFACE
        '''
        cells = self._cells
        radius = self._radii[icell]
        points = cells.points - cells.points[icell,:]
        # Find points that are inside the radius and not already in the clump
        nbors = (np.sum(points**2,1) < radius**2) * (self._inclump == 0)
        return np.where(nbors)[0]
    
if __name__=="__main__":
    sim = hamusims["IMF1_04"]
    snap = sim.Snapshots()[20].RawData()
    regions = FindHIIRegions(snap)
