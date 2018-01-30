'''
Find the stellar objects in a snapshot
Sam Geen, June 2016
'''

#from startup import *

import matplotlib as mpl
mpl.use('Agg')  

import os, sys, glob
import numpy as np

import pymses
from pymses.utils import constants as C

Msuning = 1.9891e33 

class StellarSnap(object):
    # A container for sink properties
    def __init__(self,i=[],x=[],y=[],z=[],mass=[],time=[],tcreated=[],lifetime=[],snapid=[]):
        self.sinkid = np.atleast_1d(i)
        self.sinkid = self.sinkid.astype(np.int64)
        self.x = np.atleast_1d(x)
        self.y = np.atleast_1d(y)
        self.z = np.atleast_1d(z)
        self.mass = np.atleast_1d(mass)
        self.time = np.atleast_1d(time)
        self.tcreated = np.atleast_1d(tcreated)
        self.lifetime = np.atleast_1d(lifetime)
        self.snapid = np.atleast_1d(snapid)
        self._dict = None
        self._Setup()

    def _Setup(self):
        self._dict = {"sinkid":self.sinkid,
                      "x": self.x,
                      "y": self.y,
                      "z": self.z,
                      "mass": self.mass,
                      "time": self.time,
                      "tcreated": self.tcreated,
                      "lifetime": self.lifetime,
                      "snapid": self.snapid}
    
    def __getitem__(self,key):
        return self._dict[key]

    def __repr__(self):
        return "StellarSnap: "+str(self._dict)

    def IsEmpty(self):
        return self.Length() == 0

    def Length(self):
        return len(np.atleast_1d(self._dict["mass"]))

    def Append(self,lstars):
        for key in self._dict.keys():
            self._dict[key] = np.append(self._dict[key],lstars._dict[key])
        return self

    def OnlyMass(self,mass):
        empty = StellarSnap()
        mask = self["mass"] == mass
        for key in self._dict.keys():
            if key != "time" and key != "snapid":
                select = self[key][mask]
            else:
                select = self[key]
            try:
                empty._dict[key] = np.atleast_1d(select)
            except:
                import pdb; pdb.set_trace()
        return empty

    def FixPosition(self,sinks):
        '''
        Make the position the host sink position and not a static value
        '''
        if sinks.Length() == 0:
            return
        sid = self["sinkid"]
        sids = sinks["id"]
        sx = sinks["x"]
        sy = sinks["y"]
        sz = sinks["z"]
        for i in range(0,len(sid)):
            mask = sids == sid[i]
            self["x"][i] = sx[mask]
            self["y"][i] = sy[mask]
            self["z"][i] = sz[mask]

def FindStellar(snap):
    try:
        ro = snap.RawData()
    except:
        ro = snap
    onum = str(ro.iout).zfill(5)
    stellarfile = ro.output_repos+"/output_"+onum+"/stellar_"+onum+".csv"
    unitt = ro.info["unit_time"].express(C.Myr)
    time = ro.info["time"]*unitt
    snapid = ro.iout
    # Particle mass units don't include boxlen (thanks RAMSES)
    unitm = ro.info["unit_mass"].express(C.g)/ro.info["boxlen"]**3.0
    unitm /= Msuning
    praw = []
    # Stellar particle CSV file contents:
    # id, mass, x, y, z, creation time, lifetime
    if os.path.exists(stellarfile):
        if os.stat(stellarfile).st_size > 0:
            praw = np.genfromtxt(stellarfile, delimiter=',')
    pid = []
    pmass = []
    px = []
    py = []
    pz = []
    tcreated = []
    lifetime = []
    if len(praw) > 0:
        if len(praw.shape) == 1:
            pid = np.array(praw[0])
            pmass = np.array(praw[1])
            px = np.array([praw[2]])
            py = np.array([praw[3]])
            pz = np.array([praw[4]])
            tcreated = np.array([praw[5]])
            lifetime = np.array([praw[6]])
        else:
            pid = praw[:,0]
            pmass = praw[:,1]
            px = praw[:,2]
            py = praw[:,3]
            pz = praw[:,4]
            tcreated = praw[:,5]
            lifetime = praw[:,6]
        stellars = StellarSnap(pid,px,py,pz,pmass*unitm,time,tcreated*unitt,lifetime*unitt,snapid)
    else:
        stellars = StellarSnap() # Make an empty object
    return stellars

if __name__=="__main__":
    outs = glob.glob("output_?????")
    outs.sort()
    for out in [outs[-1]]:
        ro = pymses.RamsesOutput("./",int(out[-5:]))
        stellars = FindStellar(ro)
        print stellars.mass
        print stellars.sinkid
