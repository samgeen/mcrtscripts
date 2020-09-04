'''
Find the sinks in a snapshot
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

class SinkSnap(object):
    # A container for sink properties
    def __init__(self,i=[],x=[],y=[],z=[],vx=[],vy=[],vz=[],mass=[],time=[]):
        self.id = np.atleast_1d(i)
        self.id = self.id.astype(np.int64)
        self.x = np.atleast_1d(x)
        self.y = np.atleast_1d(y)
        self.z = np.atleast_1d(z)
        self.vx = np.atleast_1d(vx)
        self.vy = np.atleast_1d(vy)
        self.vz = np.atleast_1d(vz)
        self.mass = np.atleast_1d(mass)
        self.time = np.atleast_1d(time)
        self._dict = None
        self._Setup()

    def _Setup(self):
        self._dict = {"id":self.id,
                      "x": self.x,
                      "y": self.y,
                      "z": self.z,
                      "vx": self.vx,
                      "vy": self.vy,
                      "vz": self.vz,
                      "mass": self.mass,
                      "time": self.time}

    def __getitem__(self,key):
        return self._dict[key]

    def IsEmpty(self):
        return self.Length() == 0

    def Length(self):
        return len(np.atleast_1d(self.mass))
        

def FindSinks(snap):
    try:
        ro = snap.RawData()
    except:
        ro = snap
    onum = str(ro.iout).zfill(5)
    sinkfile = ro.output_repos+"/output_"+onum+"/sink_"+onum+".csv"
    time = ro.info["time"]*ro.info["unit_time"].express(C.Myr)
    # Particle mass units don't include boxlen (thanks RAMSES)
    unitm = ro.info["unit_mass"].express(C.g)/ro.info["boxlen"]**3.0
    unitm /= Msuning
    praw = []
    if os.path.exists(sinkfile):
        if os.stat(sinkfile).st_size > 0:
            praw = np.genfromtxt(sinkfile, delimiter=',')
    pid = []
    pmass = []
    px = []
    py = []
    pz = []
    if len(praw) > 0:
        if len(praw.shape) == 1:
            pid = np.array(praw[0])
            pmass = np.array(praw[1])
            px = np.array([praw[3]])
            py = np.array([praw[4]])
            pz = np.array([praw[5]])
            vx = np.array([praw[6]])
            vy = np.array([praw[7]])
            vz = np.array([praw[8]])
        else:
            pid = praw[:,0]
            pmass = praw[:,1]
            px = praw[:,3]
            py = praw[:,4]
            pz = praw[:,5]
            vx = praw[:,6]
            vy = praw[:,7]
            vz = praw[:,8]
        sinks = SinkSnap(pid,px,py,pz,vx,vy,vz,pmass,time)
    else:
        sinks = SinkSnap() # Make an empty object
    return sinks
                
if __name__=="__main__":
    outs  = glob.glob("output_?????")
    outs.sort()
    for out in outs:
        print(out)
        ro = pymses.RamsesOutput("./",int(out[-5:]))
        sinks = FindSinks(ro)
        print(np.sum(sinks.mass))
        print(sinks.mass)
        print("---")
